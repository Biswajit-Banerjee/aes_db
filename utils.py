from datetime import datetime as dt
import os
from pathlib import Path

def log(message, level="info"):
    print(f"[{dt.now().strftime('%x %X')}] | {level:^8}| {message}")


def modify_extension(path, ext, dir_path=""):
    base_name = os.path.splitext(os.path.basename(path))[0]
    if not dir_path:
        new_path = os.path.join(os.path.dirname(path), f"{base_name}.{ext}")
    else:
        new_path = os.path.join(dir_path, f"{base_name}.{ext}")
    
    os.makedirs(os.path.dirname(new_path), exist_ok=True)
    return new_path

class PathManager:
    def __init__(self):
        super().__setattr__("base_path", Path("./Data"))
        self.base_path.mkdir(parents=True, exist_ok=True)

    def __setattr__(self, name: str, value: str):
        path = self.base_path / value

        if "." not in path.name:
            path.mkdir(parents=True, exist_ok=True)

        else:
            path.parent.mkdir(parents=True, exist_ok=True)
        super().__setattr__(name, path)

path_manager = PathManager()