#!/usr/bin/env python3
import os
import sys
import json
import logging
import hail as hl

def load_config(config_path="config/config.json"):
    """Loads configuration from a JSON file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found at: {config_path}")
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Inject Bucket into output paths dynamically if needed
    bucket = os.getenv("WORKSPACE_BUCKET")
    if bucket:
        # Helper to prepend bucket to relative output paths
        for key in config['outputs']:
            if not config['outputs'][key].startswith("gs://"):
                config['outputs'][key] = f"{bucket}/{config['outputs'][key]}"
    
    return config

def get_env_var(name: str) -> str:
    val = os.getenv(name)
    if not val:
        raise EnvironmentError(f"Required environment variable {name} is not set.")
    return val

def init_hail(log_prefix: str):
    hl.init(default_reference='GRCh38', log=f'/tmp/hail_{log_prefix}.log')

def setup_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(name)
