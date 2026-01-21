#!/usr/bin/env python3
import os
import subprocess
import logging
from utils import setup_logger

logger = setup_logger("install_tools")

def install_gcta():
    tools_dir = "tools"
    os.makedirs(tools_dir, exist_ok=True)
    
    gcta_path = os.path.join(tools_dir, "gcta64")
    if os.path.exists(gcta_path):
        logger.info("GCTA is already installed.")
        return

    logger.info("Installing GCTA...")
    url = "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.4-linux-kernel-3-x86_64.zip"
    zip_name = "gcta.zip"
    
    try:
        subprocess.run(["wget", "-O", zip_name, url], check=True)
        subprocess.run(["unzip", "-o", zip_name], check=True)
        
        # Find binary and move it
        extracted_folder = "gcta-1.94.4-linux-kernel-3-x86_64"
        src = os.path.join(extracted_folder, "gcta64")
        
        if os.path.exists(src):
            subprocess.run(["mv", src, gcta_path], check=True)
            subprocess.run(["chmod", "+x", gcta_path], check=True)
            logger.info("GCTA installed successfully.")
        else:
            logger.error(f"Could not find gcta64 in {extracted_folder}")
            
        # Cleanup
        if os.path.exists(zip_name):
            os.remove(zip_name)
        if os.path.exists(extracted_folder):
            subprocess.run(["rm", "-rf", extracted_folder], check=True)
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Installation failed: {e}")

if __name__ == "__main__":
    install_gcta()
