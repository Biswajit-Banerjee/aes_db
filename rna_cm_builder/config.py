import os
from pathlib import Path

class Config:
    BASE_DIR = Path(os.getenv('RNA_CM_BUILDER_BASE', './cm_builder_storage'))
    DATA_DIR = BASE_DIR / 'data'
    RESULTS_DIR = BASE_DIR / 'results'
    LOG_DIR = BASE_DIR / 'logs'

    @classmethod
    def create_directories(cls):
        for dir_path in [cls.BASE_DIR, cls.DATA_DIR, cls.RESULTS_DIR, cls.LOG_DIR]:
            dir_path.mkdir(parents=True, exist_ok=True)

# Create directories when the module is imported
Config.create_directories()

class Binaries:
    BASE_DIR   = Path("../executables")
    CM_COMPARE = BASE_DIR / "cmcompare"
    IP_KNOT    = BASE_DIR / "ipknot"
    
    
    '''The below is getting used only for CMCOMPARE as it requires older Covariance model version'''
    INF_102 = Path("/home/sumon/repos/infernal-1.0.2")
    
    '''The below is  '''
    INF_115 = Path("/home/sumon/repos/rna/infernal-1.1.5")
    CM_CONVERT = INF_115 / "src/cmconvert"
    
    HMMER = Path("/home/sumon/repos/rna/hmmer-3.4/src")
    HMM_BUILD  = HMMER / "hmmbuild"
    