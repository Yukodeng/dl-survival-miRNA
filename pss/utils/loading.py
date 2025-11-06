import os
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



def load_simulate_survival_data(batchNormType,
                                dataName,
                                keywords='',
                                keep_batch=False,
                                batch_col='batch.id'):    
    # prepare simulated survival data
    keywords = [keywords] if not isinstance(keywords, list) else keywords
    
    surv_folder = os.path.join('data', batchNormType, dataName)
    mirna_folder = os.path.join('data', batchNormType)
    
    found_surv = False
    for file in os.listdir(surv_folder):
        if all([kw in file for kw in keywords+['simSurvival', 'train']]):
            surv_train = pd.read_csv(os.path.join(surv_folder, file))
            found_surv = True      
        if all([kw in file for kw in keywords+['simSurvival', 'test']]):
            surv_test = pd.read_csv(os.path.join(surv_folder, file))
            found_surv = True    
    if not found_surv:
        raise FileExistsError("No survival data file found with the given query keywords")
    
    found_mirna = False
    for file in os.listdir(mirna_folder):
        if all([kw in file for kw in keywords+['simGeneExp', 'train']]):
            x_train = pd.read_csv(os.path.join(mirna_folder, file))
            found_mirna = True      
        if all([kw in file for kw in keywords+['simGeneExp', 'test']]):
            x_test = pd.read_csv(os.path.join(mirna_folder, file))
            found_mirna = True    
    if not found_mirna:
        raise FileExistsError("No gene expression data file found with the given query keywords")
    
    train_df = pd.concat([surv_train, x_train], axis=1)
    test_df = pd.concat([surv_test, x_test], axis=1)
    
    if not keep_batch:
        train_df = train_df.drop(columns=[batch_col])
        test_df = test_df.drop(columns=[batch_col])
    
    return train_df, test_df


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
