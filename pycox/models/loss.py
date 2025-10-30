import numpy as np
import torch
from torch import Tensor

# def _reduction(loss: Tensor, reduction: str = 'mean') -> Tensor:
#     if reduction == 'none':
#         return loss
#     elif reduction == 'mean':
#         return loss.mean()
#     elif reduction == 'sum':
#         return loss.sum()
#     raise ValueError(f"`reduction` = {reduction} is not valid. Use 'none', 'mean' or 'sum'.")


def cox_ph_loss_sorted(log_h: Tensor, events: Tensor, eps: float = 1e-7) -> Tensor:
    """Requires the input to be sorted by descending duration time.
    See DatasetDurationSorted.

    We calculate the negative log of $(\frac{h_i}{\sum_{j \in R_i} h_j})^d$,
    where h = exp(log_h) are the hazards and R is the risk set, and d is event.

    We just compute a cumulative sum, and not the true Risk sets. This is a
    limitati`on, but simple and fast.
    """
    if events.dtype is torch.bool:
        events = events.float()
    events = events.view(-1)
    log_h = log_h.view(-1)
    if events.sum() == 0:
        return log_h.sum() * 0.0  # update 08/11/25: safe dummy loss
    
    gamma = log_h.max()
    log_cumsum_h = log_h.sub(gamma).exp().cumsum(0).add(eps).log().add(gamma)
    return - log_h.sub(log_cumsum_h).mul(events).sum().div(events.sum())


def cox_ph_loss(log_h: Tensor, durations: Tensor, events: Tensor, eps: float = 1e-7) -> Tensor:
    """Loss for CoxPH model. If data is sorted by descending duration, see `cox_ph_loss_sorted`.

    We calculate the negative log of $(\frac{h_i}{\sum_{j \in R_i} h_j})^d$,
    where h = exp(log_h) are the hazards and R is the risk set, and d is event.

    We just compute a cumulative sum, and not the true Risk sets. This is a
    limitation, but simple and fast.
    """
    idx = durations.sort(descending=True)[1]
    events = events[idx]
    log_h = log_h[idx]
    return cox_ph_loss_sorted(log_h, events, eps)


####### [UPDATE] 07/07/2025
def stratified_cox_ph_loss(log_h: Tensor, durations: Tensor, events: Tensor, batch_indices: Tensor, eps: float = 1e-7) -> Tensor:
    """
    Stratified CoxPH loss that computes partial likelihood across batches.

    Arguments:
        log_h {torch.Tensor} -- Log hazard predictions for each instance.
        durations {torch.Tensor} -- Duration times for each instance.
        events {torch.Tensor} -- Event indicators (1 if event, 0 if censored).
        batch_indices {numpy array} -- Batch labels for each instance.
        eps {float} -- Small epsilon for numerical stability.

    Returns:
        torch.Tensor -- The total stratified negative log partial likelihood.
    """
    device = batch_indices.device
    unique_batches = torch.unique(batch_indices)
    losses = torch.zeros(len(unique_batches), device=device)
        
    for i, batch in enumerate(unique_batches):
        # Select data for the current batch
        mask = (batch_indices == batch)
        if mask.sum() == 0 or events[mask].sum() == 0:
            continue  # skip empty batch
        
        # Sort by descending durations
        idx = torch.argsort(durations[mask], descending=True)       
        
        events_batch = events[mask][idx]
        log_h_batch = log_h[mask][idx]
        if events_batch.sum() == 0:
            continue 
        
        losses[i] = cox_ph_loss_sorted(log_h_batch, events_batch, eps)
    
    total_loss = losses.sum()
    return total_loss if total_loss.requires_grad else log_h.sum() * 0.0


class CoxPHLossSorted(torch.nn.Module):
    """Loss for CoxPH.
    Requires the input to be sorted by descending duration time.
    See DatasetDurationSorted.

    We calculate the negative log of $(\frac{h_i}{\sum_{j \in R_i} h_j})^d$,
    where h = exp(log_h) are the hazards and R is the risk set, and d is event.

    We just compute a cumulative sum, and not the true Risk sets. This is a
    limitation, but simple and fast.
    """
    def __init__(self):
        super().__init__()

    def forward(self, log_h: Tensor, events: Tensor) -> Tensor:
        return cox_ph_loss_sorted(log_h, events)


class CoxPHLoss(torch.nn.Module):
    """Loss for CoxPH model. If data is sorted by descending duration, see `cox_ph_loss_sorted`.

    We calculate the negative log of $(\frac{h_i}{\sum_{j \in R_i} h_j})^d$,
    where h = exp(log_h) are the hazards and R is the risk set, and d is event.

    We just compute a cumulative sum, and not the true Risk sets. This is a
    limitation, but simple and fast.
    """
    def forward(self, log_h: Tensor, durations: Tensor, events: Tensor) -> Tensor:
        return cox_ph_loss(log_h, durations, events)

## Update 07/07/2025
class CoxPHLossStratified(torch.nn.Module):
    """Loss for CoxPH model with batch variable.

    We calculate the batch-stratified negative log of $(\frac{h_i}{\sum_{j \in R_i} h_j})^d$,
    where h = exp(log_h) are the hazards and R is the risk set, and d is event.

    We just compute a cumulative sum, and not the true Risk sets. This is a
    limitation, but simple and fast.
    """
    # def forward(self, log_h: Tensor, durations: Tensor, events: Tensor, batch_indices: Tensor) -> Tensor:
        # return stratified_cox_ph_loss(log_h, durations, events, batch_indices)
    def forward(self, log_h: Tensor, durations: Tensor, events: Tensor, batch_indices: Tensor) -> Tensor:
        if torch.isnan(log_h).any():
            print("NaNs detected in log hazards")
        if torch.isnan(durations).any():
            print("NaNs detected in input survival time")
        if torch.isnan(events).any():
            print("NaNs detected in input events")
        if (events.sum() == 0).item():
            print("No observed events in batch (val)")
        return stratified_cox_ph_loss(log_h, durations, events, batch_indices)