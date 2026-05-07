import os
import re
import glob
import sys
import pandas as pd
import logging


# --------------------- Load Data ---------------------------

class StreamToLogger:
    """Redirects sys.stdout and sys.stderr to a logger."""
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''
        
    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())
            
    def flush(self):
        pass  # Needed for compatibility


def create_logger(log_path, verbose_modules=['sklearn_pandas']):
    """
    Set up logger and redirect stdout/stderr to it.
    
    Args:
        config (dict): Config dict containing 'batchNormType' and 'dataName'.
        log_dir (str): Directory where logs are saved.
        verbose_modules (list): Optional list of module names to quiet.
    
    Returns:
        logging.Logger: Configured logger object
    """ 
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger(__name__)
    
    # Quiet specified submodules
    if verbose_modules is None:
        verbose_modules = []
    for mod in verbose_modules:
        logging.getLogger(mod).setLevel(logging.WARNING)
        
    # Redirect print() output
    sys.stdout = StreamToLogger(logger, logging.INFO)
    sys.stderr = StreamToLogger(logger, logging.ERROR)
    
    return logger


# ------- Updated data loading function -------

def _find_one(pattern: str) -> str:
    hits = glob.glob(pattern)
    if len(hits) == 0:
        raise FileNotFoundError(f"No file matches: {pattern}")
    if len(hits) > 1:
        # pick the most recent if multiple
        hits.sort(key=os.path.getmtime, reverse=True)
    return hits[0]

def _parse_be_flags(batchNormType: str) -> tuple[int, int]:
    """
    Parse 'BEXX...' where XX are train/test flags (0=clean, 1=dirty).
    Example: 'BE10...' => train=1 (dirty), test=0 (clean)
    """
    m = re.match(r"^BE([01])([01])", batchNormType)
    if not m:
        raise ValueError(f"batchNormType must start with 'BE00/01/10/11'. Got: {batchNormType}")
    return int(m.group(1)), int(m.group(2))

def _col_slice(flag: int) -> slice:
    return list(range(538)) if flag == 0 else list(range(538, 1076))

def load_prepared_split(
    batchNormType: str,
    dataName: str,
    train_size: int,
    iter_i: int,
    *,
    keep_batch: bool = True,
    time_col: str = "time",
    status_col: str = "status",
    batch_col: str = "batch_id"
    ):    

    # [NOTE] Survival outcome saved per exp condition per iteration
    # --> reason being conditions with surv-batch correlation will reorganize sample order 
    surv_folder = os.path.join('data', batchNormType, dataName, f"train{train_size}")
    file_pat = "_".join([batchNormType, dataName, f"train{train_size}", f"iter{iter_i}"])
    try:
        train_pat = os.path.join(surv_folder, f"simSurv_train_{file_pat}.csv")
        test_pat  = os.path.join(surv_folder, f"simSurv_test_{file_pat}.csv")
        train_path = _find_one(train_pat)
        test_path  = _find_one(test_pat)
        surv_train = pd.read_csv(train_path)[[time_col, status_col]].reset_index(drop=True)
        surv_test = pd.read_csv(test_path)[[time_col, status_col]].reset_index(drop=True)
    except FileExistsError:
        raise FileExistsError(f"No train / test index data found under: {surv_folder}")
   
    # Load partitioned log count data (w batch id!)
    mirna_tr_folder = os.path.join('data', "gene_exp", batchNormType, f"train{train_size}")
    mirna_te_folder = os.path.join('data', "gene_exp", batchNormType)
    mirna_tr_file = "_".join([batchNormType, f"train{train_size}", f"iter{iter_i}"])
    mirna_te_file = "_".join([batchNormType, f"iter{iter_i}"])
    try: 
        mirna_tr_pat = os.path.join(mirna_tr_folder, f"simGeneExp_train_{mirna_tr_file}.parquet")
        mirna_te_pat = os.path.join(mirna_te_folder, f"simGeneExp_test_{mirna_te_file}.parquet")
        mirna_tr_path = _find_one(mirna_tr_pat)
        mirna_te_path = _find_one(mirna_te_pat)
        x_train = pd.read_parquet(mirna_tr_path)
        x_test  = pd.read_parquet(mirna_te_path)
    except FileExistsError:
        raise FileExistsError("No gene expression data file found")
    
    train_df = pd.concat([surv_train, x_train], axis=1)
    test_df = pd.concat([surv_test, x_test], axis=1)
    
    if not keep_batch:
        train_df = train_df.drop(columns=[batch_col])
        test_df = test_df.drop(columns=[batch_col])
    
    return train_df, test_df

# def load_simulate_survival_data(batchNormType,
#                                 dataName,
#                                 keywords='',
#                                 keep_batch=False,
#                                 batch_col='batch.id'):    
#     # prepare simulated survival data
#     keywords = [keywords] if not isinstance(keywords, list) else keywords
    
#     surv_folder = os.path.join('data', batchNormType, dataName)
#     mirna_folder = os.path.join('data', batchNormType)
    
#     found_surv = False
#     for file in os.listdir(surv_folder):
#         if all([kw in file for kw in keywords+['simSurvival', 'train']]):
#             surv_train = pd.read_csv(os.path.join(surv_folder, file))
#             found_surv = True      
#         if all([kw in file for kw in keywords+['simSurvival', 'test']]):
#             surv_test = pd.read_csv(os.path.join(surv_folder, file))
#             found_surv = True    
#     if not found_surv:
#         raise FileExistsError("No survival data file found with the given query keywords")
    
#     found_mirna = False
#     for file in os.listdir(mirna_folder):
#         if all([kw in file for kw in keywords+['simGeneExp', 'train']]):
#             x_train = pd.read_csv(os.path.join(mirna_folder, file))
#             found_mirna = True      
#         if all([kw in file for kw in keywords+['simGeneExp', 'test']]):
#             x_test = pd.read_csv(os.path.join(mirna_folder, file))
#             found_mirna = True    
#     if not found_mirna:
#         raise FileExistsError("No gene expression data file found with the given query keywords")
    
#     train_df = pd.concat([surv_train, x_train], axis=1)
#     test_df = pd.concat([surv_test, x_test], axis=1)
    
#     if not keep_batch:
#         train_df = train_df.drop(columns=[batch_col])
#         test_df = test_df.drop(columns=[batch_col])
    
#     return train_df, test_df


modelname_dict = {'coxnet':'CoxPH',
                'svm': 'SSVM',
                'rsf': 'RSF',
                'gb':  "SGB",
                'deepsurv-torch': "DeepSurv"}

def load_simulate_results(dataFolderName, 
                        modelnames=['coxnet','svm','rsf','gb','deepsurv-torch'],
                        fileName = 'model.results.10runs.txt'):
    
    results = pd.DataFrame({'n train': []})
    for mdl in modelnames:
        file_dir = os.path.join('models', dataFolderName, mdl, fileName)
        try:
            result_df = pd.read_table(file_dir, index_col=0)
        except:
            result_df = pd.read_table(os.path.join('models', dataFolderName, mdl, 'model.results.txt'), index_col=0)
        result_df['model'] = modelname_dict[mdl]
        results = pd.concat([results, result_df], axis=0)

    return results
