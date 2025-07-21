import sys
import os
import warnings
warnings.filterwarnings("ignore")

import argparse
import json
import logging
import torch
from datetime import datetime
from runDeepSurvModels import DeepSurvPipeline
from utils import *


# class StreamToLogger:
#     """Redirects sys.stdout and sys.stderr to a logger."""
#     def __init__(self, logger, log_level=logging.INFO):
#         self.logger = logger
#         self.log_level = log_level
#         self.linebuf = ''
        
#     def write(self, buf):
#         for line in buf.rstrip().splitlines():
#             self.logger.log(self.log_level, line.rstrip())
            
#     def flush(self):
#         pass  # Needed for compatibility


# def create_logger(log_path, verbose_modules=['sklearn_pandas']):
#     """
#     Set up logger and redirect stdout/stderr to it.
    
#     Args:
#         config (dict): Config dict containing 'batchNormType' and 'dataName'.
#         log_dir (str): Directory where logs are saved.
#         verbose_modules (list): Optional list of module names to quiet.
    
#     Returns:
#         logging.Logger: Configured logger object
#     """ 
#     logging.basicConfig(
#         level=logging.INFO,
#         format="%(asctime)s - %(levelname)s - %(message)s",
#         handlers=[
#             logging.FileHandler(log_path),
#             logging.StreamHandler()
#         ]
#     )
    
#     logger = logging.getLogger(__name__)
    
#     # Quiet specified submodules
#     if verbose_modules is None:
#         verbose_modules = []
#     for mod in verbose_modules:
#         logging.getLogger(mod).setLevel(logging.WARNING)
        
#     # Redirect print() output
#     sys.stdout = StreamToLogger(logger, logging.INFO)
#     sys.stderr = StreamToLogger(logger, logging.ERROR)
    
#     return logger


def main(config_path):
    """
    Example config file
    ---------------------------------------------------------------
    {
        "batchNormType": "BE00Asso00_normNone",
        "dataName": "linear-moderate",
        "keywords": [""],
        "time_col": "time",
        "status_col": "status",
        "batch_col": "batch.id",
        "keep_batch": true,
        "hyperparameters": {
            "num_nodes": {"type": "categorical", "choices": [[64, 64], [32, 32], [16,16]]},
            "dropout": {"type": "float", "low": 0.1, "high": 0.7},
            "learning_rate": {"type": "float", "low": 1e-4, "high": 1e-2, "log": true},
            "weight_decay": {"type": "float", "low": 1e-5, "high": 1e-2, "log": true},
            "batch_size": {"type": "categorical", "choices": [128, 64, 32, 16]},
            "activation": {"type": "categorical", "choices": ["ReLU", "LeakyReLU", "SELU"]}
        },
        "subset_sizes": [100, 500, 1000, 2000, 5000],
        "runs_per_size": [20, 20, 20, 20, 20],
        "splits_per_size": [3, 5, 5, 10, 10],
        "trials_per_size": [30, 30, 30, 30, 30],
        "trial_threshold": 30,
        "is_tune": true,
        "is_save": true,
        "file_name": "model_results.csv"
    }
    """
    # ============= Load Config ==============  
    with open(config_path) as f:
        config = json.load(f)
        
    # ============ Setup Logging =============
    today = datetime.now().strftime("%Y%m%d")
    log_filename = f"{config['batchNormType']}-{config['dataName']}_{today}.log"
    os.makedirs("logs", exist_ok=True)
    log_path = os.path.join("logs", log_filename)
    
    logger = create_logger(log_path=log_path)
    
    # ======== Check GPU Availability ========
    if torch.cuda.is_available():
        logger.info("GPU is available")
    else:
        logger.info("GPU is not available, using CPU instead")
    
    # =============== Load Data =============== 
    logger.info("Loading miRNAseq dataset %s-%s...", config['batchNormType'], config['dataName'])
    train_df, test_df = load_simulate_survival_data(
        batchNormType=config["batchNormType"],
        dataName=config["dataName"],
        keywords=config["keywords"], 
        keep_batch=True
    )
    logger.info("Train shape: %s | Test shape: %s", train_df.shape, test_df.shape)
    
    # ============ Training Settings ============ 
    # Initialize pipeline
    pipeline = DeepSurvPipeline(
        train_df=train_df,
        test_df=test_df,
        time_col=config['time_col'],
        status_col=config['status_col'],
        batch_col=config['batch_col'],
        batchNormType=config['batchNormType'],
        dataName=config['dataName'],
        hyperparameters=config['hyperparameters'],
        is_stratified=config.get('stratified', True),
        storage_url=config['storage_url']
    )
    if "param_override" in config:
        pipeline.param_override = config['param_override']
        logger.info("Applying parameter overrides: %s", pipeline.param_override)
        
    # Optional override for early stopping settings
    if "early_stop_per_size" in config:
        pipeline.early_stop_per_size = config["early_stop_per_size"]
        
    # ============ Run Training Pipeline ============ 
    logger.info("Launching training pipeline with subset sizes %s...", config["subset_sizes"])
    logger.info("Training with %s CoxPH neural network" % ['non-stratified', 'stratified'][pipeline.is_stratified])
    
    pipeline.train_over_subsets(
        subset_sizes=config['subset_sizes'],
        runs_per_size=config['runs_per_size'],
        splits_per_size=config['splits_per_size'],
        trials_per_size=config['trials_per_size'],
        trial_threshold=config['trial_threshold'],
        n_jobs=config.get('n_jobs', 1),
        is_tune=config.get('is_tune', True),
        is_save=config.get('is_save', True),
        fileName=config.get('fileName', None)
    )
    
    logger.info("Congrats! Training completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PRECISION.seq.survival.DL - Run DeepSurv training pipeline.')
    parser.add_argument('-c', '--config', type=str, required=True, help='Path to JSON config file')
    
    args = parser.parse_args()
    
    main(args.config)
    
    logging.shutdown()
    sys.exit(0)

