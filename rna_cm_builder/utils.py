import logging
from datetime import datetime
import os
from pathlib import Path
from .config import Config

def setup_logger(name, log_file, level=logging.INFO):
    formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
    
    # Ensure the directory exists
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger

logger = setup_logger('rna_cm_builder', Config.LOG_DIR / 'rna_cm_builder.log')

def log(message, level="info"):
    getattr(logger, level.lower())(message)

def modify_extension(path, ext, dir_path=""):
    base_name = os.path.splitext(os.path.basename(path))[0]
    if not dir_path:
        new_path = os.path.join(os.path.dirname(path), f"{base_name}.{ext}")
    else:
        new_path = os.path.join(dir_path, f"{base_name}.{ext}")
    
    os.makedirs(os.path.dirname(new_path), exist_ok=True)
    return new_path

class PathManager:
    def __init__(self, base_path=Config.DATA_DIR):
        super().__setattr__("base_path", Path(base_path))
        self.base_path.mkdir(parents=True, exist_ok=True)

    def __setattr__(self, name: str, value: str):
        path = self.base_path / value

        if "." not in path.name:
            path.mkdir(parents=True, exist_ok=True)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
        super().__setattr__(name, path)

path_manager = PathManager()