"""
Utility functions for the Molecular AOP Analyser application.
"""
import os
import time
import logging
from config import Config

logger = logging.getLogger(__name__)

def cleanup_file(filepath):
    """
    Safely remove a file if it exists.
    
    Args:
        filepath (str): Path to the file to remove
    
    Returns:
        bool: True if file was removed or didn't exist, False if error occurred
    """
    try:
        if os.path.exists(filepath):
            os.remove(filepath)
            logger.info(f"Cleaned up file: {filepath}")
            return True
        return True
    except OSError as e:
        logger.error(f"Failed to cleanup file {filepath}: {e}")
        return False

def cleanup_old_uploads(max_age_hours=24):
    """
    Remove old files from the uploads directory.
    
    Args:
        max_age_hours (int): Maximum age of files to keep in hours
    
    Returns:
        int: Number of files cleaned up
    """
    uploads_dir = Config.UPLOAD_FOLDER
    if not os.path.exists(uploads_dir):
        return 0
    
    cleaned_count = 0
    current_time = time.time()
    cutoff_time = current_time - (max_age_hours * 3600)
    
    try:
        for filename in os.listdir(uploads_dir):
            filepath = os.path.join(uploads_dir, filename)
            if os.path.isfile(filepath):
                file_age = os.path.getmtime(filepath)
                if file_age < cutoff_time:
                    if cleanup_file(filepath):
                        cleaned_count += 1
                        
    except OSError as e:
        logger.error(f"Error during cleanup of uploads directory: {e}")
    
    if cleaned_count > 0:
        logger.info(f"Cleaned up {cleaned_count} old upload files")
    
    return cleaned_count

def validate_file_path(filepath):
    """
    Validate that a file path is safe and within expected directories.
    
    Args:
        filepath (str): File path to validate
    
    Returns:
        bool: True if path is valid and safe
    """
    try:
        # Resolve the path and check if it's within uploads directory
        abs_filepath = os.path.abspath(filepath)
        abs_uploads = os.path.abspath(Config.UPLOAD_FOLDER)
        
        # Check if file is within uploads directory
        if not abs_filepath.startswith(abs_uploads):
            logger.warning(f"File path outside uploads directory: {filepath}")
            return False
            
        # Check if file exists
        if not os.path.exists(abs_filepath):
            logger.warning(f"File does not exist: {filepath}")
            return False
            
        return True
        
    except Exception as e:
        logger.error(f"Error validating file path {filepath}: {e}")
        return False