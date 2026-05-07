from __future__ import annotations
import os
import sys
import logging
import argparse
from ..utils.loading import create_logger
from ..run_models import MLSurvivalPipeline, train_over_subsets

import warnings
warnings.filterwarnings("ignore")


def run_with_args(args):
    
    # ============ Setup Logging =============
    os.makedirs("logs", exist_ok=True)
    log_filename = f"{args.modelType}-{args.batchNormType}-{args.dataName}.log"
    log_path = os.path.join("logs", log_filename)
    logger = create_logger(log_path=log_path)
    

    # ============ Run Training Pipeline ============ 
    
    # Initialize survival model class
    ml = MLSurvivalPipeline(
        model_type=args.modelType,
        train_df=None, test_df=None,
        batchNormType=args.batchNormType,
        dataName=args.dataName,
        is_stratified=args.is_stratified
    )
    
    subset_sizes = [int(x) for x in args.subsets.split(",")] if args.subsets else [100, 200, 500, 1000, 2000, 5000, 10000]
    runs_per_size = [int(x) for x in args.runs.split(",")] if args.runs else [20]*len(subset_sizes)
    run_seeds = [int(x) for x in args.seeds.split(",")] if args.seeds else None
    splits_per_size = [int(x) for x in args.splits.split(",")] if args.splits else [3,5,5,5,10,10,10]
    trials_per_size = [int(x) for x in args.trials.split(",")] if args.trials else [20]*len(subset_sizes)
    
    logger.info("Launching training pipeline with subset sizes %s..." % subset_sizes)
    logger.info("Training with %s%s model" % (['', 'stratified '][ml.is_stratified], args.modelType))
    
    _ = train_over_subsets(
        ml,
        subset_sizes=subset_sizes,
        runs_per_size=runs_per_size,
        run_seeds=run_seeds,
        splits_per_size=splits_per_size,
        trials_per_size=trials_per_size,
        trial_threshold=args.trial_threshold,
        n_jobs=args.n_jobs,
        is_tune=args.is_tune,
        is_save=args.is_save,
        file_path=args.filePath,
        file_name=args.fileName
    )
    
    logger.info("Congrats! Training completed.")
    
    logging.shutdown()
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser(description="PRECISION.seq.survival - Run ML survival model pipelines.")
    parser.add_argument("--batchNormType", type=str, required=True, help="Batch effects and normalization type (e.g., BE00Asso00_normNone)")
    parser.add_argument("--dataName", type=str, required=True, help="Survival simulation type (e.g., linear-moderate)")
    parser.add_argument("--modelType", type=str, required=True, choices=["coxnet", "ssvm", "rsf", "sgb"], help="Model type")
    parser.add_argument("--keywords", type=str, default="", help="Keywords for input file search, comma-separated (optional)")
    parser.add_argument("--keep_batch", action="store_true", help="Keep batch column when loading input data", default=True)
    parser.add_argument("--subsets", type=str, help="Comma-separated training subset sizes (e.g., 100,500,1000,2000,5000,10000)")
    parser.add_argument("--runs", type=str, help="Comma-separated runs per size (default 20 each)")
    parser.add_argument("--seeds", type=str, default=None, help="Comma-separated list of integer seeds for each run (e.g., ``5,6,7,8,9``)")
    parser.add_argument("--splits", type=str, help="Comma-separated k-fold splits per size (e.g., 3,5,5,10,10,10)")
    parser.add_argument("--trials", type=str, help="Comma-separated trials per size (default 30 each)")
    parser.add_argument("--trial_threshold", type=int, default=30, help="Trial threshold for tuning reuse (default 30)")
    parser.add_argument("--n_jobs", type=int, help="Number of parallel jobs for Optuna (default 1)", default=1)
    parser.add_argument("--is_tune", action="store_true", help="Enable hyperparameter tuning with optuna (default: False)", default=False)
    parser.add_argument("--is_save", action="store_true", help="Save results to file (default: False)", default=False)
    parser.add_argument("--is_stratified", action="store_true", help="Fit batch-stratified models by inputting batch labels", default=False)
    parser.add_argument("--filePath", type=str, default=None, help="Output file pathway")
    parser.add_argument("--fileName", type=str, default=None, help="Output file name")

    args = parser.parse_args()
    run_with_args(args)


if __name__ == "__main__":
    main()
