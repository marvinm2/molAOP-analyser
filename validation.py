"""
Input validation utilities for the Molecular AOP Analyser application.
"""
import re
import logging
from typing import Dict, Any, List, Tuple, Union
from werkzeug.datastructures import FileStorage
from config import Config

logger = logging.getLogger(__name__)

class ValidationError(Exception):
    """Custom exception for validation errors."""
    pass

def validate_file_upload(file: FileStorage) -> Tuple[bool, str]:
    """
    Validate uploaded file.
    
    Args:
        file: Werkzeug FileStorage object
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not file or not file.filename:
        return False, "No file was uploaded"
    
    # Check file extension
    if '.' not in file.filename:
        return False, "File must have an extension"
    
    extension = file.filename.rsplit('.', 1)[1].lower()
    if extension not in Config.ALLOWED_EXTENSIONS:
        return False, f"File type '{extension}' not allowed. Use: {', '.join(Config.ALLOWED_EXTENSIONS)}"
    
    # Check file size (approximate, as we can't get exact size from FileStorage easily)
    file.seek(0, 2)  # Seek to end
    size = file.tell()
    file.seek(0)  # Reset to beginning
    
    if size > Config.MAX_FILE_SIZE:
        return False, f"File too large. Maximum size is {Config.MAX_FILE_SIZE // (1024*1024)} MB"
    
    if size == 0:
        return False, "File is empty"
    
    return True, ""

def validate_column_name(column_name: str) -> Tuple[bool, str]:
    """
    Validate column name for safety.
    
    Args:
        column_name: Column name to validate
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not column_name or not isinstance(column_name, str):
        return False, "Column name is required"
    
    # Check length
    if len(column_name) > 100:
        return False, "Column name too long"
    
    # Check for potentially dangerous characters
    if re.search(r'[<>;"\'\\]', column_name):
        return False, "Column name contains invalid characters"
    
    return True, ""

def validate_threshold(threshold_str: str) -> Tuple[bool, Union[float, str]]:
    """
    Validate and parse log2FC threshold.
    
    Args:
        threshold_str: String representation of threshold
    
    Returns:
        Tuple of (is_valid, threshold_value_or_error_message)
    """
    if not threshold_str:
        return True, 0.0  # Default value
    
    try:
        threshold = float(threshold_str)
        
        # Reasonable bounds check
        if threshold < 0:
            return False, "Threshold must be non-negative"
        
        if threshold > 10:
            return False, "Threshold seems unreasonably high (>10)"
        
        return True, threshold
        
    except ValueError:
        return False, "Threshold must be a valid number"

def validate_aop_selection(aop_id: str) -> Tuple[bool, str]:
    """
    Validate AOP selection.
    
    Args:
        aop_id: AOP identifier
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not aop_id:
        return False, "AOP selection is required"
    
    # Check if AOP exists in our configuration
    valid_aops = [aop_data.get('id') for aop_data in Config.CASE_STUDY_AOPS.values() 
                  if aop_data.get('id') and aop_data.get('enabled', True)]
    
    if aop_id not in valid_aops:
        return False, f"Invalid AOP selection: {aop_id}"
    
    return True, ""

def validate_form_data(form_data: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Validate complete form data for analysis.
    
    Args:
        form_data: Dictionary of form data
    
    Returns:
        Tuple of (is_valid, list_of_error_messages)
    """
    errors = []
    
    # Required fields
    required_fields = ['filename', 'id_column', 'fc_column', 'pval_column', 'aop_selection']
    for field in required_fields:
        if not form_data.get(field):
            errors.append(f"Required field '{field}' is missing")
    
    # Validate individual fields if present
    if form_data.get('id_column'):
        is_valid, error = validate_column_name(form_data['id_column'])
        if not is_valid:
            errors.append(f"Gene ID column: {error}")
    
    if form_data.get('fc_column'):
        is_valid, error = validate_column_name(form_data['fc_column'])
        if not is_valid:
            errors.append(f"Log2FC column: {error}")
    
    if form_data.get('pval_column'):
        is_valid, error = validate_column_name(form_data['pval_column'])
        if not is_valid:
            errors.append(f"P-value column: {error}")
    
    if form_data.get('logfc_threshold'):
        is_valid, result = validate_threshold(form_data['logfc_threshold'])
        if not is_valid:
            errors.append(f"Log2FC threshold: {result}")
    
    if form_data.get('aop_selection'):
        is_valid, error = validate_aop_selection(form_data['aop_selection'])
        if not is_valid:
            errors.append(error)
    
    return len(errors) == 0, errors

def sanitize_filename(filename: str) -> str:
    """
    Sanitize filename for safe storage.
    
    Args:
        filename: Original filename
    
    Returns:
        str: Sanitized filename
    """
    # Remove path traversal attempts
    filename = filename.replace('..', '').replace('/', '').replace('\\', '')
    
    # Remove potentially dangerous characters
    filename = re.sub(r'[<>:"|?*]', '_', filename)
    
    # Ensure reasonable length
    if len(filename) > 100:
        name, ext = filename.rsplit('.', 1) if '.' in filename else (filename, '')
        filename = f"{name[:90]}.{ext}" if ext else name[:100]
    
    return filename

def log_validation_error(error_type: str, details: str, user_data: Dict[str, Any] = None):
    """
    Log validation errors for security monitoring.
    
    Args:
        error_type: Type of validation error
        details: Detailed error message
        user_data: Sanitized user data (no sensitive info)
    """
    log_data = {
        'error_type': error_type,
        'details': details,
        'user_data': {k: str(v)[:100] for k, v in (user_data or {}).items()}  # Truncate for safety
    }
    
    logger.warning(f"Validation error: {log_data}")