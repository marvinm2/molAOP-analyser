"""
Unit tests for column detection service.
"""

import pytest
import pandas as pd
from services.column_detector import ColumnDetector, ColumnMatch, ColumnSuggestions


@pytest.mark.unit
class TestColumnDetector:
    """Test column detection functionality."""
    
    def test_column_detector_initialization(self, column_detector):
        """Test ColumnDetector initialization."""
        assert column_detector.min_confidence == 0.3
        assert column_detector.high_confidence == 0.8
        assert len(column_detector.GENE_ID_PATTERNS) > 0
        assert len(column_detector.LOG2FC_PATTERNS) > 0
        assert len(column_detector.PVALUE_PATTERNS) > 0
    
    def test_gene_id_column_detection(self, column_detector):
        """Test gene ID column detection."""
        # Create test DataFrame with clear gene ID column
        df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
            'expression': [1.5, -2.1, 3.2, -1.8, 2.4],
            'p_value': [0.001, 0.005, 0.0001, 0.01, 0.03]
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Should detect Gene_Symbol as gene ID column
        assert suggestions.best_gene_id is not None
        assert suggestions.best_gene_id.column_name == 'Gene_Symbol'
        assert suggestions.best_gene_id.confidence > 0.5
        assert 'gene' in ' '.join(suggestions.best_gene_id.reasons).lower()
    
    def test_log2fc_column_detection(self, column_detector):
        """Test log2 fold change column detection."""
        df = pd.DataFrame({
            'genes': ['GENE1', 'GENE2', 'GENE3'],
            'log2FoldChange': [2.5, -1.8, 3.2],
            'pval': [0.001, 0.005, 0.01]
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Should detect log2FoldChange column
        assert suggestions.best_log2fc is not None
        assert suggestions.best_log2fc.column_name == 'log2FoldChange'
        assert suggestions.best_log2fc.confidence > 0.5
    
    def test_pvalue_column_detection(self, column_detector):
        """Test p-value column detection."""
        df = pd.DataFrame({
            'genes': ['GENE1', 'GENE2', 'GENE3'],
            'fold_change': [2.0, -1.5, 1.8],
            'padj': [0.001, 0.05, 0.1]
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Should detect padj as p-value column
        assert suggestions.best_pvalue is not None
        assert suggestions.best_pvalue.column_name == 'padj'
        assert suggestions.best_pvalue.confidence > 0.5
    
    def test_complete_column_detection(self, column_detector, sample_gene_data):
        """Test detection of all column types together."""
        suggestions = column_detector.detect_columns(sample_gene_data)
        
        # Should detect all three column types
        assert suggestions.best_gene_id is not None
        assert suggestions.best_log2fc is not None
        assert suggestions.best_pvalue is not None
        
        # Check column names are correct
        assert suggestions.best_gene_id.column_name == 'Gene_Symbol'
        assert suggestions.best_log2fc.column_name == 'log2FoldChange'
        assert suggestions.best_pvalue.column_name == 'padj'
    
    def test_column_name_scoring(self, column_detector):
        """Test column name pattern scoring."""
        # Test exact matches
        score = column_detector._score_column_name('gene', column_detector.GENE_ID_PATTERNS)
        assert score == 1.0  # Exact match
        
        # Test partial matches
        score = column_detector._score_column_name('gene_symbol', column_detector.GENE_ID_PATTERNS)
        assert score > 0.5  # Should match gene pattern
        
        # Test no match
        score = column_detector._score_column_name('random_column', column_detector.GENE_ID_PATTERNS)
        assert score == 0.0
    
    def test_gene_id_content_analysis(self, column_detector):
        """Test gene ID content analysis."""
        # Test with clear HGNC symbols
        df = pd.DataFrame({
            'good_genes': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
            'bad_column': ['random', 'text', 'values', 'here', 'test']
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Should detect good_genes despite non-ideal name
        assert suggestions.best_gene_id is not None
        assert suggestions.best_gene_id.column_name == 'good_genes'
        
        # Should have content-based confidence
        assert any('identifier' in reason.lower() or 'hgnc' in reason.lower() 
                  for reason in suggestions.best_gene_id.reasons)
    
    def test_log2fc_value_range_analysis(self, column_detector):
        """Test log2FC value range analysis."""
        # Test with typical log2FC values
        df = pd.DataFrame({
            'genes': ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'],
            'good_fc': [2.5, -1.8, 3.2, -2.1, 1.9],  # Good log2FC range
            'bad_fc': [100, 200, 300, 400, 500]       # Bad range
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Should prefer good_fc over bad_fc
        assert suggestions.best_log2fc is not None
        assert suggestions.best_log2fc.column_name == 'good_fc'
        
        # Should mention range in reasons
        assert any('range' in reason.lower() for reason in suggestions.best_log2fc.reasons)
    
    def test_pvalue_range_analysis(self, column_detector):
        """Test p-value range analysis."""
        df = pd.DataFrame({
            'genes': ['GENE1', 'GENE2', 'GENE3'],
            'good_pval': [0.001, 0.05, 0.1],    # Good p-value range
            'bad_pval': [1.5, 2.3, -0.1]       # Bad range (outside [0,1])
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Should prefer good_pval
        assert suggestions.best_pvalue is not None
        assert suggestions.best_pvalue.column_name == 'good_pval'
    
    def test_low_confidence_detection(self, column_detector):
        """Test handling of low confidence detections."""
        # Create DataFrame with ambiguous column names and content
        df = pd.DataFrame({
            'col1': ['A', 'B', 'C'],
            'col2': [1, 2, 3],
            'col3': [0.5, 0.6, 0.7]
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # With ambiguous data, confidence should be low or None
        if suggestions.best_gene_id:
            assert suggestions.best_gene_id.confidence < 0.8
        if suggestions.best_log2fc:
            assert suggestions.best_log2fc.confidence < 0.8
        if suggestions.best_pvalue:
            assert suggestions.best_pvalue.confidence < 0.8
    
    def test_empty_dataframe(self, column_detector):
        """Test behavior with empty DataFrame."""
        empty_df = pd.DataFrame()
        
        suggestions = column_detector.detect_columns(empty_df)
        
        assert suggestions.best_gene_id is None
        assert suggestions.best_log2fc is None
        assert suggestions.best_pvalue is None
        assert len(suggestions.gene_id_suggestions) == 0
        assert len(suggestions.log2fc_suggestions) == 0
        assert len(suggestions.pvalue_suggestions) == 0
    
    def test_missing_values_handling(self, column_detector):
        """Test handling of columns with missing values."""
        df = pd.DataFrame({
            'Gene_Symbol': ['BRCA1', 'TP53', None, 'MYC', 'KRAS'],
            'log2FoldChange': [2.5, None, 3.2, -2.1, 1.9],
            'padj': [0.001, 0.005, 0.01, None, 0.03]
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Should still detect columns despite missing values
        assert suggestions.best_gene_id is not None
        assert suggestions.best_log2fc is not None
        assert suggestions.best_pvalue is not None
    
    def test_column_match_dataclass(self):
        """Test ColumnMatch dataclass functionality."""
        match = ColumnMatch(
            column_name='test_column',
            column_type='gene_id',
            confidence=0.85,
            reasons=['Pattern match', 'Content analysis'],
            data_analysis={'sample_data': ['A', 'B', 'C']}
        )
        
        assert match.column_name == 'test_column'
        assert match.column_type == 'gene_id'
        assert match.confidence == 0.85
        assert len(match.reasons) == 2
        assert match.data_analysis is not None
    
    def test_column_suggestions_sorting(self, column_detector):
        """Test that suggestions are sorted by confidence."""
        # Create DataFrame with multiple potential gene ID columns
        df = pd.DataFrame({
            'genes': ['GENE1', 'GENE2', 'GENE3'],        # Ambiguous name
            'Gene_Symbol': ['BRCA1', 'TP53', 'EGFR'],    # Clear name + content
            'identifier': ['ID1', 'ID2', 'ID3']          # Good name, poor content
        })
        
        suggestions = column_detector.detect_columns(df)
        
        # Suggestions should be sorted by confidence (highest first)
        if len(suggestions.gene_id_suggestions) > 1:
            confidences = [s.confidence for s in suggestions.gene_id_suggestions]
            assert confidences == sorted(confidences, reverse=True)
    
    def test_confidence_description(self, column_detector):
        """Test confidence level descriptions."""
        assert column_detector.get_confidence_description(0.9) == "High"
        assert column_detector.get_confidence_description(0.7) == "Medium"
        assert column_detector.get_confidence_description(0.4) == "Low"
        assert column_detector.get_confidence_description(0.1) == "Very Low"