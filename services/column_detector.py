"""
Intelligent column detection service for gene expression datasets.
Automatically identifies gene ID, log2FC, and p-value columns.
"""
import re
import logging
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy import stats
from services.gene_id_validator import gene_id_validator, GeneIDAnalysis

logger = logging.getLogger(__name__)

@dataclass
class ColumnMatch:
    """Represents a potential column match with confidence score."""
    column_name: str
    column_type: str  # 'gene_id', 'log2fc', 'pvalue'
    confidence: float
    reasons: List[str]
    data_analysis: Optional[Dict] = None

@dataclass
class ColumnSuggestions:
    """Complete set of column suggestions for a dataset."""
    gene_id_suggestions: List[ColumnMatch]
    log2fc_suggestions: List[ColumnMatch]
    pvalue_suggestions: List[ColumnMatch]
    
    # Best matches (highest confidence)
    best_gene_id: Optional[ColumnMatch] = None
    best_log2fc: Optional[ColumnMatch] = None 
    best_pvalue: Optional[ColumnMatch] = None

class ColumnDetector:
    """Intelligent column detector for gene expression data."""
    
    # Column name patterns for different data types
    GENE_ID_PATTERNS = [
        r'gene.*symbol', r'gene.*name', r'gene.*id', r'symbol', r'hgnc',
        r'ensembl', r'entrez', r'ncbi', r'gene', r'identifier', r'id',
        r'probe.*id', r'probe', r'feature'
    ]
    
    LOG2FC_PATTERNS = [
        r'log2.*fc', r'log2.*fold', r'logfc', r'log.*fold.*change',
        r'fold.*change', r'fc', r'lfc', r'log2.*ratio', r'ratio'
    ]
    
    PVALUE_PATTERNS = [
        r'p.*val', r'pval', r'p\.val', r'p_val', r'p.*adj', r'padj',
        r'fdr', r'q.*val', r'qval', r'significance', r'p.*value'
    ]
    
    def __init__(self):
        """Initialize the column detector."""
        self.min_confidence = 0.3
        self.high_confidence = 0.8
        
    def detect_columns(self, df: pd.DataFrame) -> ColumnSuggestions:
        """
        Detect all relevant columns in a DataFrame.
        
        Args:
            df: Input DataFrame to analyze
            
        Returns:
            ColumnSuggestions with ranked matches for each column type
        """
        gene_id_matches = self._detect_gene_id_columns(df)
        log2fc_matches = self._detect_log2fc_columns(df)
        pvalue_matches = self._detect_pvalue_columns(df)
        
        # Sort by confidence (highest first)
        gene_id_matches.sort(key=lambda x: x.confidence, reverse=True)
        log2fc_matches.sort(key=lambda x: x.confidence, reverse=True)
        pvalue_matches.sort(key=lambda x: x.confidence, reverse=True)
        
        # Determine best matches
        best_gene_id = gene_id_matches[0] if gene_id_matches and gene_id_matches[0].confidence >= self.min_confidence else None
        best_log2fc = log2fc_matches[0] if log2fc_matches and log2fc_matches[0].confidence >= self.min_confidence else None
        best_pvalue = pvalue_matches[0] if pvalue_matches and pvalue_matches[0].confidence >= self.min_confidence else None
        
        return ColumnSuggestions(
            gene_id_suggestions=gene_id_matches,
            log2fc_suggestions=log2fc_matches,
            pvalue_suggestions=pvalue_matches,
            best_gene_id=best_gene_id,
            best_log2fc=best_log2fc,
            best_pvalue=best_pvalue
        )
    
    def _detect_gene_id_columns(self, df: pd.DataFrame) -> List[ColumnMatch]:
        """Detect potential gene ID columns."""
        matches = []
        
        for col_name in df.columns:
            match = self._analyze_gene_id_column(df, col_name)
            if match and match.confidence >= self.min_confidence:
                matches.append(match)
                
        return matches
    
    def _detect_log2fc_columns(self, df: pd.DataFrame) -> List[ColumnMatch]:
        """Detect potential log2 fold change columns."""
        matches = []
        
        for col_name in df.columns:
            match = self._analyze_log2fc_column(df, col_name)
            if match and match.confidence >= self.min_confidence:
                matches.append(match)
                
        return matches
    
    def _detect_pvalue_columns(self, df: pd.DataFrame) -> List[ColumnMatch]:
        """Detect potential p-value columns."""
        matches = []
        
        for col_name in df.columns:
            match = self._analyze_pvalue_column(df, col_name)
            if match and match.confidence >= self.min_confidence:
                matches.append(match)
                
        return matches
    
    def _analyze_gene_id_column(self, df: pd.DataFrame, col_name: str) -> Optional[ColumnMatch]:
        """Analyze a potential gene ID column."""
        column_data = df[col_name].dropna()
        
        if len(column_data) == 0:
            return None
        
        reasons = []
        confidence_scores = []
        
        # Name-based scoring
        name_score = self._score_column_name(col_name, self.GENE_ID_PATTERNS)
        if name_score > 0:
            reasons.append(f"Column name matches gene ID patterns (score: {name_score:.2f})")
            confidence_scores.append(name_score * 0.4)  # 40% weight for name
        
        # Content-based analysis using gene ID validator
        try:
            gene_analysis = gene_id_validator.analyze_gene_column(column_data)
            content_confidence = gene_analysis.confidence
            
            if content_confidence > 0.5:
                reasons.append(f"Contains {gene_analysis.primary_type} identifiers ({content_confidence:.1%} valid)")
                confidence_scores.append(content_confidence * 0.6)  # 60% weight for content
            
            # Bonus for supported ID types
            if gene_id_validator.is_supported_type(gene_analysis.primary_type):
                reasons.append("ID type supported for analysis")
                confidence_scores.append(0.1)  # 10% bonus
            
            data_analysis = {
                'gene_id_analysis': gene_analysis,
                'sample_data': column_data.head(5).tolist()
            }
            
        except Exception as e:
            logger.warning(f"Failed to analyze gene column {col_name}: {e}")
            content_confidence = 0
            data_analysis = None
        
        # Calculate final confidence
        final_confidence = sum(confidence_scores) if confidence_scores else 0
        final_confidence = min(final_confidence, 1.0)  # Cap at 100%
        
        if final_confidence < self.min_confidence:
            return None
            
        return ColumnMatch(
            column_name=col_name,
            column_type='gene_id',
            confidence=final_confidence,
            reasons=reasons,
            data_analysis=data_analysis
        )
    
    def _analyze_log2fc_column(self, df: pd.DataFrame, col_name: str) -> Optional[ColumnMatch]:
        """Analyze a potential log2FC column."""
        column_data = df[col_name].dropna()
        
        if len(column_data) == 0:
            return None
            
        reasons = []
        confidence_scores = []
        
        # Name-based scoring
        name_score = self._score_column_name(col_name, self.LOG2FC_PATTERNS)
        if name_score > 0:
            reasons.append(f"Column name matches fold change patterns (score: {name_score:.2f})")
            confidence_scores.append(name_score * 0.4)
        
        # Content-based analysis
        try:
            # Convert to numeric, skip non-numeric values
            numeric_data = pd.to_numeric(column_data, errors='coerce').dropna()
            
            if len(numeric_data) < len(column_data) * 0.8:
                reasons.append("Contains significant non-numeric values")
                return None  # Too many non-numeric values
            
            content_score = 0
            
            # Check value range (typical log2FC range: -10 to +10)
            data_range = numeric_data.max() - numeric_data.min()
            mean_val = numeric_data.mean()
            std_val = numeric_data.std()
            
            if -15 <= numeric_data.min() <= 15 and -15 <= numeric_data.max() <= 15:
                content_score += 0.3
                reasons.append("Values in typical log2FC range")
            
            # Check for distribution around zero (expected for fold changes)
            if abs(mean_val) < 2:  # Mean close to zero
                content_score += 0.2
                reasons.append("Mean close to zero (expected for fold changes)")
            
            # Check for reasonable spread
            if 0.5 <= std_val <= 5:
                content_score += 0.2
                reasons.append("Standard deviation in expected range")
            
            # Check for both positive and negative values
            if (numeric_data > 0).any() and (numeric_data < 0).any():
                content_score += 0.2
                reasons.append("Contains both up and down regulation")
            
            confidence_scores.append(content_score * 0.6)  # 60% weight for content
            
            data_analysis = {
                'numeric_stats': {
                    'mean': float(mean_val),
                    'std': float(std_val),
                    'min': float(numeric_data.min()),
                    'max': float(numeric_data.max()),
                    'range': float(data_range)
                },
                'sample_data': numeric_data.head(5).tolist()
            }
            
        except Exception as e:
            logger.warning(f"Failed to analyze log2FC column {col_name}: {e}")
            data_analysis = None
        
        final_confidence = sum(confidence_scores) if confidence_scores else 0
        final_confidence = min(final_confidence, 1.0)
        
        if final_confidence < self.min_confidence:
            return None
            
        return ColumnMatch(
            column_name=col_name,
            column_type='log2fc',
            confidence=final_confidence,
            reasons=reasons,
            data_analysis=data_analysis
        )
    
    def _analyze_pvalue_column(self, df: pd.DataFrame, col_name: str) -> Optional[ColumnMatch]:
        """Analyze a potential p-value column."""
        column_data = df[col_name].dropna()
        
        if len(column_data) == 0:
            return None
            
        reasons = []
        confidence_scores = []
        
        # Name-based scoring
        name_score = self._score_column_name(col_name, self.PVALUE_PATTERNS)
        if name_score > 0:
            reasons.append(f"Column name matches p-value patterns (score: {name_score:.2f})")
            confidence_scores.append(name_score * 0.4)
        
        # Content-based analysis
        try:
            numeric_data = pd.to_numeric(column_data, errors='coerce').dropna()
            
            if len(numeric_data) < len(column_data) * 0.8:
                reasons.append("Contains significant non-numeric values")
                return None
            
            content_score = 0
            
            # Check if values are in [0, 1] range
            if (numeric_data >= 0).all() and (numeric_data <= 1).all():
                content_score += 0.4
                reasons.append("All values in [0, 1] range")
            
            # Check for typical p-value distribution (more small values)
            small_pvals = (numeric_data < 0.05).sum()
            medium_pvals = ((numeric_data >= 0.05) & (numeric_data < 0.5)).sum()
            large_pvals = (numeric_data >= 0.5).sum()
            
            if small_pvals > 0 and large_pvals > 0:
                content_score += 0.2
                reasons.append("Contains mix of significant and non-significant values")
            
            # Check for absence of exactly 0 or 1 (which are suspicious for p-values)
            exact_zeros = (numeric_data == 0).sum()
            exact_ones = (numeric_data == 1).sum()
            
            if exact_zeros == 0 and exact_ones == 0:
                content_score += 0.1
                reasons.append("No exact 0 or 1 values (good for p-values)")
            
            confidence_scores.append(content_score * 0.6)
            
            data_analysis = {
                'pvalue_stats': {
                    'min': float(numeric_data.min()),
                    'max': float(numeric_data.max()),
                    'mean': float(numeric_data.mean()),
                    'significant_count': int(small_pvals),
                    'total_count': len(numeric_data)
                },
                'sample_data': numeric_data.head(5).tolist()
            }
            
        except Exception as e:
            logger.warning(f"Failed to analyze p-value column {col_name}: {e}")
            data_analysis = None
        
        final_confidence = sum(confidence_scores) if confidence_scores else 0
        final_confidence = min(final_confidence, 1.0)
        
        if final_confidence < self.min_confidence:
            return None
            
        return ColumnMatch(
            column_name=col_name,
            column_type='pvalue',
            confidence=final_confidence,
            reasons=reasons,
            data_analysis=data_analysis
        )
    
    def _score_column_name(self, col_name: str, patterns: List[str]) -> float:
        """Score column name against pattern list."""
        col_lower = col_name.lower().strip()
        max_score = 0
        
        for pattern in patterns:
            if re.search(pattern, col_lower):
                # Exact match gets higher score
                if re.fullmatch(pattern, col_lower):
                    max_score = max(max_score, 1.0)
                else:
                    max_score = max(max_score, 0.8)
        
        return max_score
    
    def get_confidence_description(self, confidence: float) -> str:
        """Get human-readable confidence description."""
        if confidence >= self.high_confidence:
            return "High"
        elif confidence >= 0.6:
            return "Medium"
        elif confidence >= self.min_confidence:
            return "Low"
        else:
            return "Very Low"

# Global detector instance
column_detector = ColumnDetector()