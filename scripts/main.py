import os
import sys
import json
import logging
import argparse
from datetime import datetime
import torch
from pss.utils import load_simulate_survival_data, create_logger
from pss.run_models import DeepSurvPipeline, train_over_subsets

import warnings
warnings.filterwarnings("ignore")


def main(config_path):
    """
    Example config file:
    ----------------------------------------------
    {
        "batchNormType": "BE00Asso00_normNone",
        "dataName": "linear-moderate",
        "keywords": [""],
        "time_col": "time",
        "status_col": "status",
        "batch_col": "batch.id",
        "keep_batch": true,
        "hyperparameters": {
            "num_nodes": {"type": "categorical", "choices": [[64, 32], [32, 16], [32]]},
            "dropout": {"type": "float", "low": 0.1, "high": 0.7},
            "learning_rate": {"type": "float", "low": 1e-4, "high": 1e-2, "log": true},
            "weight_decay": {"type": "float", "low": 1e-5, "high": 1e-2, "log": true},
            "batch_size": {"type": "categorical", "choices": [128, 64, 32]},
            "activation": {"type": "categorical", "choices": ["ReLU", "LeakyReLU", "SELU"]}
        },
        "subset_sizes": [100, 500, 1000, 2000, 5000, 10000],
        "runs_per_size": [20, 20, 20, 20, 20, 20],
        "splits_per_size": [3, 5, 5, 10, 10, 10],
        "trials_per_size": [30, 30, 30, 30, 30, 30],
        "trial_threshold": 30,
        "is_tune": true,
        "is_save": true,
        "fileName": "model_results.csv"
    }
    """
    # =============== Load Config ================  
    with open(config_path) as f:
        config = json.load(f)
        
    # ============== Setup Logging ===============
    # today = datetime.now().strftime("%m%d%y")
    log_filename = f"{config['batchNormType']}-{config['dataName']}.log"
    os.makedirs("logs", exist_ok=True)
    log_path = os.path.join("logs", log_filename)
    logger = create_logger(log_path=log_path)
    
    # ========== Check GPU Availability ==========
    if torch.cuda.is_available():
        logger.info("GPU is available")
    else:
        logger.info("GPU is not available, using CPU instead")
    
    # ================= Load Data ================ 
    logger.info("Loading miRNAseq dataset %s-%s...", config['batchNormType'], config['dataName'])
    train_df, test_df = load_simulate_survival_data(
        batchNormType=config["batchNormType"],
        dataName=config["dataName"],
        keywords=config["keywords"], 
        keep_batch=True
    )
    logger.info("Train shape: %s | Test shape: %s", train_df.shape, test_df.shape)
    
    # ============ Training Settings ============= 
    # Initialize pipeline
    dl = DeepSurvPipeline(
        train_df=train_df,
        test_df=test_df,
        batchNormType=config['batchNormType'],
        dataName=config['dataName'],
        is_stratified=config.get('stratified', False),
        time_col=config['time_col'],
        status_col=config['status_col'],
        batch_col=config['batch_col'],
        hyperparameters=config['hyperparameters'],
        storage_url=config['storage_url']
    )
    if "param_override" in config:
        dl.param_override = config['param_override']
        logger.info("Applying parameter overrides: %s", dl.param_override)
        
    # Optional override for early stopping settings
    if "early_stop_per_size" in config:
        dl.early_stop_per_size = config["early_stop_per_size"]
        
    # ============ Run Training Pipeline ============ 
    logger.info("Launching training pipeline with subset sizes %s...", config["subset_sizes"])
    logger.info("Training with %s CoxPH neural network" % ['non-stratified', 'stratified'][dl.is_stratified])
    
    _ = train_over_subsets(
        dl,
        subset_sizes=config['subset_sizes'],
        runs_per_size=config['runs_per_size'],
        splits_per_size=config['splits_per_size'],
        trials_per_size=config['trials_per_size'],
        trial_threshold=config['trial_threshold'],
        n_jobs=config.get('n_jobs', 1),
        is_tune=config.get('is_tune', True),
        is_save=config.get('is_save', True),
        file_path=config.get('filePath', None),
        file_name=config.get('fileName', None)
    )
    
    logger.info("Congrats! Training completed.")
    
    logging.shutdown()
    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PRECISION.seq.survival.DL - Run DeepSurv training pipeline.')
    parser.add_argument('-c', '--config', type=str, required=True, help='Path to JSON config file')
    
    args = parser.parse_args()
    main(args.config)


