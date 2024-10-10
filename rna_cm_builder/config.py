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