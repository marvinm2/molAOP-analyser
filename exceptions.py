"""
Custom exception classes for the Molecular AOP Analyser application.
"""

class AOPAnalysisError(Exception):
    """Base exception for AOP analysis errors."""
    
    def __init__(self, message: str, error_code: str = None, details: dict = None):
        self.message = message
        self.error_code = error_code or "GENERAL_ERROR"
        self.details = details or {}
        super().__init__(self.message)

class DataValidationError(AOPAnalysisError):
    """Raised when input data validation fails."""
    
    def __init__(self, message: str, field: str = None, value: str = None):
        self.field = field
        self.value = value
        super().__init__(
            message, 
            error_code="VALIDATION_ERROR",
            details={"field": field, "value": value}
        )

class FileProcessingError(AOPAnalysisError):
    """Raised when file processing fails."""
    
    def __init__(self, message: str, filename: str = None, line_number: int = None):
        self.filename = filename
        self.line_number = line_number
        super().__init__(
            message,
            error_code="FILE_ERROR", 
            details={"filename": filename, "line_number": line_number}
        )

class EnrichmentAnalysisError(AOPAnalysisError):
    """Raised when enrichment analysis fails."""
    
    def __init__(self, message: str, ke_id: str = None, gene_count: int = None):
        self.ke_id = ke_id
        self.gene_count = gene_count
        super().__init__(
            message,
            error_code="ENRICHMENT_ERROR",
            details={"ke_id": ke_id, "gene_count": gene_count}
        )

class AOPDataError(AOPAnalysisError):
    """Raised when AOP reference data is missing or invalid."""
    
    def __init__(self, message: str, aop_id: str = None, missing_data: str = None):
        self.aop_id = aop_id
        self.missing_data = missing_data
        super().__init__(
            message,
            error_code="AOP_DATA_ERROR", 
            details={"aop_id": aop_id, "missing_data": missing_data}
        )

class NetworkBuildError(AOPAnalysisError):
    """Raised when network visualization building fails."""
    
    def __init__(self, message: str, node_count: int = None, edge_count: int = None):
        self.node_count = node_count
        self.edge_count = edge_count
        super().__init__(
            message,
            error_code="NETWORK_ERROR",
            details={"node_count": node_count, "edge_count": edge_count}
        )

class ConfigurationError(AOPAnalysisError):
    """Raised when application configuration is invalid."""
    
    def __init__(self, message: str, config_key: str = None):
        self.config_key = config_key
        super().__init__(
            message,
            error_code="CONFIG_ERROR",
            details={"config_key": config_key}
        )

# User-friendly error messages
USER_ERROR_MESSAGES = {
    "VALIDATION_ERROR": "Please check your input data and try again.",
    "FILE_ERROR": "There was a problem processing your file. Please ensure it's a valid CSV/TSV file with the correct format.",
    "ENRICHMENT_ERROR": "The enrichment analysis could not be completed. This may be due to insufficient gene overlap or data quality issues.",
    "AOP_DATA_ERROR": "The selected AOP pathway data is not available. Please try a different pathway.",
    "NETWORK_ERROR": "The network visualization could not be generated. The analysis results are still available in the table.",
    "CONFIG_ERROR": "A system configuration issue occurred. Please contact support.",
    "GENERAL_ERROR": "An unexpected error occurred. Please try again or contact support if the problem persists."
}

def get_user_friendly_message(error_code: str) -> str:
    """
    Get user-friendly error message for an error code.
    
    Args:
        error_code: Error code from exception
        
    Returns:
        User-friendly error message
    """
    return USER_ERROR_MESSAGES.get(error_code, USER_ERROR_MESSAGES["GENERAL_ERROR"])

def format_error_response(exception: AOPAnalysisError, include_details: bool = False) -> dict:
    """
    Format exception for API response.
    
    Args:
        exception: AOPAnalysisError instance
        include_details: Whether to include technical details
        
    Returns:
        Dictionary suitable for JSON response
    """
    response = {
        "error": True,
        "error_code": exception.error_code,
        "message": get_user_friendly_message(exception.error_code),
        "technical_message": str(exception)
    }
    
    if include_details and exception.details:
        response["details"] = exception.details
    
    return response