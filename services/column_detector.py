"""
Intelligent column detection service for gene expression datasets.
Automatically identifies gene ID, log2FC, and p-value columns.
"""
import re
import logging
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
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

    # Adjusted (FDR / q-value) p-value suggestions — D-05 (Plan 11-01).
    # Independent of `pvalue_suggestions` so callers can carry both raw and
    # adjusted columns through to the gene-by-KE export. Default factory keeps
    # backward compat for any code constructing ColumnSuggestions without it.
    pvalue_adj_suggestions: List[ColumnMatch] = field(default_factory=list)

    # Best matches (highest confidence)
    best_gene_id: Optional[ColumnMatch] = None
    best_log2fc: Optional[ColumnMatch] = None
    best_pvalue: Optional[ColumnMatch] = None
    # Best adjusted p-value match (D-05). None when no adjusted column detected.
    best_pvalue_adj: Optional[ColumnMatch] = None

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

    # Raw (unadjusted) p-value patterns — D-05 (Plan 11-01).
    # Used by `_detect_pvalue_adj_columns` companion `_detect_pvalue_raw_columns`
    # to identify the *raw* column when both raw and adjusted are present.
    # Anchored exact matches keep DESeq2 ('pvalue') and limma ('P.Value' / 'p.value')
    # conventions disambiguated from adjusted variants. Case-insensitive matching
    # is handled by `_score_column_name` via `.lower()`.
    PVALUE_RAW_PATTERNS = [
        r'^pvalue$', r'^p_value$', r'^p\.value$',
        r'^pval$', r'^p_val$', r'^p\.val$',
        r'^p$', r'^pvalues$', r'^p_values$', r'^p\.values$',
        r'^raw_?p(value)?$', r'^raw\.p(value)?$',
        r'^unadjusted_?p(value)?$', r'^nominal_?p(value)?$',
    ]

    # Adjusted (FDR / corrected / q-value) p-value patterns — D-05 (Plan 11-01).
    # By convention 'fdr', 'qvalue', 'padj', and limma's 'adj.P.Val' all denote
    # an adjusted column; we treat them as a single category. Patterns are
    # anchored to avoid spurious matches against raw 'pvalue' columns.
    PVALUE_ADJ_PATTERNS = [
        r'^padj$', r'^p_adj$', r'^p\.adj$',
        r'^fdr$', r'^fdr_p(value)?$', r'^fdr\.p(value)?$',
        r'^qvalue$', r'^q_value$', r'^q\.value$',
        r'^qval$', r'^q_val$', r'^q\.val$',
        r'^adj_?p(value|val)?$', r'^adj\.p(value|val|\.val)?$',
        r'^adjusted_?p(value)?$', r'^adjusted\.p(value)?$',
        r'^corrected_?p(value)?$', r'^corrected\.p(value)?$',
        r'^bh_?p(value)?$', r'^bh\.p(value)?$',
        r'^benjamini[._-]?hochberg$',
        # Generic catch-alls for "adj.*p" / "p.*adj" style headers.
        r'adj.*p', r'p.*adj',
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
        # D-05: detect adjusted p-value column independently. The legacy
        # PVALUE_PATTERNS umbrella already covers both raw and adjusted, so
        # `best_pvalue` semantics are unchanged. `best_pvalue_adj` is purely
        # additive and consumed by the gene-by-KE export writer (Plan 11-03).
        pvalue_adj_matches = self._detect_pvalue_adj_columns(df)

        # Sort by confidence (highest first)
        gene_id_matches.sort(key=lambda x: x.confidence, reverse=True)
        log2fc_matches.sort(key=lambda x: x.confidence, reverse=True)
        pvalue_matches.sort(key=lambda x: x.confidence, reverse=True)
        pvalue_adj_matches.sort(key=lambda x: x.confidence, reverse=True)

        # Determine best matches
        best_gene_id = gene_id_matches[0] if gene_id_matches and gene_id_matches[0].confidence >= self.min_confidence else None
        best_log2fc = log2fc_matches[0] if log2fc_matches and log2fc_matches[0].confidence >= self.min_confidence else None
        best_pvalue = pvalue_matches[0] if pvalue_matches and pvalue_matches[0].confidence >= self.min_confidence else None
        best_pvalue_adj = pvalue_adj_matches[0] if pvalue_adj_matches and pvalue_adj_matches[0].confidence >= self.min_confidence else None

        return ColumnSuggestions(
            gene_id_suggestions=gene_id_matches,
            log2fc_suggestions=log2fc_matches,
            pvalue_suggestions=pvalue_matches,
            pvalue_adj_suggestions=pvalue_adj_matches,
            best_gene_id=best_gene_id,
            best_log2fc=best_log2fc,
            best_pvalue=best_pvalue,
            best_pvalue_adj=best_pvalue_adj,
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

    def _detect_pvalue_adj_columns(self, df: pd.DataFrame) -> List[ColumnMatch]:
        """Detect potential adjusted p-value (FDR / q-value / padj) columns.

        D-05 (Plan 11-01): Independent of `_detect_pvalue_columns` so the gene-by-KE
        export can echo both raw and adjusted columns from the upload as separate
        passthrough fields. Uses `PVALUE_ADJ_PATTERNS` for name matching and a
        content-scoring path that tolerates the BH ceiling (many legitimate 1.0
        values) without losing confidence.

        The companion raw-only pattern set `PVALUE_RAW_PATTERNS` is intentionally
        not consumed by this detector (which is adjusted-only by design); it is
        defined alongside as the canonical raw-only pattern catalogue for any
        downstream caller that needs to disambiguate raw from adjusted (e.g. the
        gene-by-KE export writer in Plan 11-03).
        """
        # Sanity guard: keep both pattern sets in sync structurally so that any
        # future tightening of `PVALUE_RAW_PATTERNS` does not silently bypass
        # the adjusted matcher when raw and adjusted regex sets need to remain
        # mutually exclusive on common DE-tool conventions.
        assert isinstance(self.PVALUE_RAW_PATTERNS, list) and self.PVALUE_RAW_PATTERNS, \
            "PVALUE_RAW_PATTERNS must be a non-empty list (D-05 invariant)"
        matches = []

        for col_name in df.columns:
            match = self._analyze_pvalue_adj_column(df, col_name)
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

            # Issue #69: rank by what enrichment can actually match, not by
            # whether the type is recognised. A DESeq2 export carries both
            # 'GeneID' (Ensembl) and 'GeneSymbol' (HGNC); both used to score
            # identically, so column order picked the Ensembl one and the run
            # silently intersected the gene sets at zero. Only symbol columns
            # get the bonus; other recognised types are penalised so a symbol
            # column beside them always wins.
            if gene_analysis.primary_type in gene_id_validator.ENRICHMENT_MATCHABLE_TYPES:
                reasons.append("Gene symbols — directly matchable against KE gene sets")
                confidence_scores.append(0.1)
            elif gene_id_validator.is_supported_type(gene_analysis.primary_type):
                reasons.append(
                    f"{gene_analysis.primary_type} identifiers cannot be matched against "
                    f"KE gene sets without translation — prefer a gene symbol column"
                )
                confidence_scores.append(-0.25)
            
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

    def _analyze_pvalue_adj_column(self, df: pd.DataFrame, col_name: str) -> Optional[ColumnMatch]:
        """Analyze a potential adjusted p-value column (D-05, Plan 11-01).

        Mirrors `_analyze_pvalue_column` with two material differences:
        1. Matches against `PVALUE_ADJ_PATTERNS` (FDR / qvalue / padj / adj.P.Val).
        2. Drops the "no exact 0 or 1 values" bonus — adjusted columns
           legitimately contain many 1.0 values from the BH ceiling, so
           penalising them would suppress valid padj column detection.
        """
        column_data = df[col_name].dropna()

        if len(column_data) == 0:
            return None

        reasons = []
        confidence_scores = []

        # Name-based scoring against the adjusted-only pattern set.
        # Adjusted detection REQUIRES a name match (D-05): content alone
        # (values in [0,1] + mix of significant/non-significant) describes any
        # p-value column, raw or adjusted. Without a name signal we cannot
        # claim a column is specifically "adjusted" — that would flag every
        # raw `pvalue` upload as also adjusted, breaking the only-raw case.
        name_score = self._score_column_name(col_name, self.PVALUE_ADJ_PATTERNS)
        if name_score == 0:
            return None
        reasons.append(f"Column name matches adjusted p-value patterns (score: {name_score:.2f})")
        confidence_scores.append(name_score * 0.4)

        # Content-based analysis (same as raw analyser, minus the no-1.0 bonus)
        try:
            numeric_data = pd.to_numeric(column_data, errors='coerce').dropna()

            if len(numeric_data) < len(column_data) * 0.8:
                reasons.append("Contains significant non-numeric values")
                return None

            content_score = 0

            # Adjusted p-values must be in [0, 1] just like raw p-values
            if (numeric_data >= 0).all() and (numeric_data <= 1).all():
                content_score += 0.4
                reasons.append("All values in [0, 1] range")

            # Mix of significant and non-significant indicates a real p/q column.
            small_pvals = (numeric_data < 0.05).sum()
            large_pvals = (numeric_data >= 0.5).sum()

            if small_pvals > 0 and large_pvals > 0:
                content_score += 0.2
                reasons.append("Contains mix of significant and non-significant values")

            # NB: the "no exact 1.0 values" bonus from `_analyze_pvalue_column`
            # is intentionally OMITTED here. BH adjustment routinely produces
            # padj=1.0 for unranked genes; rewarding columns lacking 1.0 would
            # bias against legitimate adjusted-p columns.

            confidence_scores.append(content_score * 0.6)

            data_analysis = {
                'pvalue_adj_stats': {
                    'min': float(numeric_data.min()),
                    'max': float(numeric_data.max()),
                    'mean': float(numeric_data.mean()),
                    'significant_count': int(small_pvals),
                    'total_count': len(numeric_data),
                    'ones_count': int((numeric_data == 1).sum()),
                },
                'sample_data': numeric_data.head(5).tolist()
            }

        except Exception as e:
            logger.warning(f"Failed to analyze adjusted p-value column {col_name}: {e}")
            data_analysis = None

        final_confidence = sum(confidence_scores) if confidence_scores else 0
        final_confidence = min(final_confidence, 1.0)

        if final_confidence < self.min_confidence:
            return None

        return ColumnMatch(
            column_name=col_name,
            column_type='pvalue_adj',
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