"""
Gene ID validation service supporting multiple identifier types.
Supports HGNC symbols, Ensembl Gene IDs, and NCBI Gene IDs.
"""
import re
import logging
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from collections import Counter
import pandas as pd

logger = logging.getLogger(__name__)

@dataclass
class ValidationResult:
    """Result of gene ID validation."""
    id_type: str
    total_ids: int
    valid_ids: int
    invalid_ids: List[str]
    confidence: float
    mixed_types: bool
    type_distribution: Dict[str, int]

@dataclass
class GeneIDAnalysis:
    """Comprehensive analysis of a gene ID column."""
    primary_type: str
    confidence: float
    validation_result: ValidationResult
    recommendations: List[str]
    warnings: List[str]

class GeneIDValidator:
    """Enhanced gene ID validator supporting multiple ID types."""
    
    # Gene ID patterns for different identifier types
    GENE_ID_PATTERNS = {
        'HGNC': re.compile(r'^[A-Z][A-Z0-9-]*[A-Z0-9]$|^[A-Z]$'),  # HGNC symbols
        'Ensembl': re.compile(r'^ENSG\d{11}(\.\d+)?$'),  # Ensembl with optional version
        'NCBI': re.compile(r'^\d+$'),  # NCBI Gene IDs (Entrez)
        'Probe': re.compile(r'^[A-Za-z0-9_-]+(_at|_s_at|_x_at)$'),  # Microarray probes
        'RefSeq': re.compile(r'^(NM_|XM_|NR_|XR_)\d+(\.\d+)?$'),  # RefSeq IDs
    }
    
    # Confidence thresholds
    HIGH_CONFIDENCE_THRESHOLD = 0.85
    MEDIUM_CONFIDENCE_THRESHOLD = 0.60
    MIXED_TYPE_THRESHOLD = 0.20  # If secondary type >20%, consider mixed
    
    def __init__(self):
        """Initialize the validator."""
        self.supported_types = {'HGNC', 'Ensembl', 'NCBI'}  # Main supported types
        
    def validate_gene_ids(self, gene_ids: List[str], sample_size: int = 200) -> ValidationResult:
        """
        Validate a list of gene identifiers.
        
        Args:
            gene_ids: List of gene identifiers to validate
            sample_size: Maximum number of IDs to analyze for performance
            
        Returns:
            ValidationResult with detailed analysis
        """
        if not gene_ids:
            return ValidationResult(
                id_type="Unknown",
                total_ids=0,
                valid_ids=0,
                invalid_ids=[],
                confidence=0.0,
                mixed_types=False,
                type_distribution={}
            )
        
        # Sample IDs for performance (analyze first N non-empty entries)
        clean_ids = [str(id_).strip().upper() for id_ in gene_ids if pd.notna(id_) and str(id_).strip()]
        sample_ids = clean_ids[:sample_size]
        
        if not sample_ids:
            return ValidationResult(
                id_type="Unknown",
                total_ids=len(gene_ids),
                valid_ids=0,
                invalid_ids=gene_ids,
                confidence=0.0,
                mixed_types=False,
                type_distribution={}
            )
        
        # Count matches for each ID type
        type_counts = Counter()
        valid_ids_set = set()
        invalid_ids = []
        
        for gene_id in sample_ids:
            matched = False
            for id_type, pattern in self.GENE_ID_PATTERNS.items():
                if pattern.match(gene_id):
                    type_counts[id_type] += 1
                    valid_ids_set.add(gene_id)
                    matched = True
                    break  # First match wins
            
            if not matched:
                invalid_ids.append(gene_id)
        
        # Determine primary ID type and confidence
        total_analyzed = len(sample_ids)
        
        if not type_counts:
            primary_type = "Unknown"
            confidence = 0.0
            mixed_types = False
        else:
            primary_type = type_counts.most_common(1)[0][0]
            primary_count = type_counts[primary_type]
            confidence = primary_count / total_analyzed
            
            # Check for mixed types
            if len(type_counts) > 1:
                secondary_count = type_counts.most_common(2)[1][1] if len(type_counts.most_common(2)) > 1 else 0
                secondary_ratio = secondary_count / total_analyzed
                mixed_types = secondary_ratio > self.MIXED_TYPE_THRESHOLD
            else:
                mixed_types = False
        
        return ValidationResult(
            id_type=primary_type,
            total_ids=total_analyzed,
            valid_ids=len(valid_ids_set),
            invalid_ids=invalid_ids[:20],  # Limit invalid examples
            confidence=confidence,
            mixed_types=mixed_types,
            type_distribution=dict(type_counts)
        )
    
    def analyze_gene_column(self, gene_series: pd.Series) -> GeneIDAnalysis:
        """
        Perform comprehensive analysis of a gene ID column.
        
        Args:
            gene_series: Pandas Series containing gene identifiers
            
        Returns:
            GeneIDAnalysis with recommendations and warnings
        """
        # Validate the gene IDs
        gene_ids = gene_series.dropna().astype(str).tolist()
        validation_result = self.validate_gene_ids(gene_ids)
        
        # Generate recommendations and warnings
        recommendations = []
        warnings = []
        
        primary_type = validation_result.id_type
        confidence = validation_result.confidence
        
        # Confidence-based recommendations
        if confidence >= self.HIGH_CONFIDENCE_THRESHOLD:
            recommendations.append(f"High confidence {primary_type} identifiers detected")
        elif confidence >= self.MEDIUM_CONFIDENCE_THRESHOLD:
            recommendations.append(f"Medium confidence {primary_type} identifiers detected")
            if primary_type not in self.supported_types:
                warnings.append(f"ID type '{primary_type}' may require conversion for analysis")
        else:
            warnings.append("Low confidence in ID type detection")
            recommendations.append("Consider checking ID format or using different column")
        
        # Mixed type warnings
        if validation_result.mixed_types:
            warnings.append("Mixed ID types detected - may cause analysis issues")
            recommendations.append("Consider standardizing to single ID type")
        
        # Invalid IDs warnings
        if validation_result.invalid_ids:
            invalid_count = len(validation_result.invalid_ids)
            invalid_percentage = (invalid_count / validation_result.total_ids) * 100
            warnings.append(f"{invalid_count} invalid/unrecognized IDs found ({invalid_percentage:.1f}%)")
            
            if invalid_percentage > 30:
                warnings.append("High percentage of invalid IDs may indicate wrong column or format")
        
        # Unsupported type warnings
        if primary_type not in self.supported_types and confidence > 0.5:
            warnings.append(f"Primary ID type '{primary_type}' requires conversion for enrichment analysis")
            recommendations.append("Consider using HGNC symbols, Ensembl IDs, or NCBI Gene IDs")
        
        return GeneIDAnalysis(
            primary_type=primary_type,
            confidence=confidence,
            validation_result=validation_result,
            recommendations=recommendations,
            warnings=warnings
        )
    
    def get_confidence_level(self, confidence: float) -> str:
        """Get human-readable confidence level."""
        if confidence >= self.HIGH_CONFIDENCE_THRESHOLD:
            return "High"
        elif confidence >= self.MEDIUM_CONFIDENCE_THRESHOLD:
            return "Medium"
        else:
            return "Low"
    
    def is_supported_type(self, id_type: str) -> bool:
        """Check if ID type is supported for analysis."""
        return id_type in self.supported_types
    
    def get_type_description(self, id_type: str) -> str:
        """Get human-readable description of ID type."""
        descriptions = {
            'HGNC': 'HGNC Gene Symbols (e.g., TP53, BRCA1)',
            'Ensembl': 'Ensembl Gene IDs (e.g., ENSG00000141510)',
            'NCBI': 'NCBI Gene IDs (e.g., 7157, 672)',
            'Probe': 'Microarray Probe IDs (e.g., 1007_s_at)',
            'RefSeq': 'RefSeq Transcript IDs (e.g., NM_000546)',
            'Unknown': 'Unrecognized ID format'
        }
        return descriptions.get(id_type, f"{id_type} identifiers")
    
    def suggest_conversion(self, id_type: str) -> Optional[str]:
        """Suggest conversion strategy for unsupported ID types."""
        if id_type in self.supported_types:
            return None
            
        suggestions = {
            'Probe': 'Convert microarray probe IDs to gene symbols using platform annotation',
            'RefSeq': 'Map RefSeq IDs to HGNC symbols or Ensembl Gene IDs',
            'Unknown': 'Verify ID format and consider using standard gene identifiers'
        }
        
        return suggestions.get(id_type, 'Consider converting to HGNC symbols or Ensembl Gene IDs')

# Global validator instance
gene_id_validator = GeneIDValidator()