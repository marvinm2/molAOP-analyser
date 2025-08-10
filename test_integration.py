"""
Integration test script for new metadata and report generation features.
Tests the complete workflow including database operations.
"""

import os
import sys
import requests
import json
from datetime import datetime

# Test configuration
BASE_URL = "http://localhost:5000"
TEST_DATA = {
    'demo_file': 'GSE90122_TO90137.tsv',
    'dataset_id': 'TEST001',
    'stressor': 'Test Chemical',
    'dosing': '10 µM for 24h', 
    'owner': 'Test User',
    'description': 'Integration test experiment'
}

def test_metadata_workflow():
    """Test the complete metadata collection and storage workflow."""
    print("🧪 Testing metadata collection workflow...")
    
    # Start session
    session = requests.Session()
    
    try:
        # Step 1: Load demo dataset with metadata
        print("  1. Loading demo dataset with metadata...")
        response = session.post(f"{BASE_URL}/preview", data=TEST_DATA)
        
        if response.status_code != 200:
            print(f"  ❌ Preview failed: {response.status_code} - {response.text}")
            return False
        
        print(f"  ✅ Dataset loaded successfully")
        
        # Step 2: Extract column suggestions from response (basic check)
        if 'column-suggestion' not in response.text.lower():
            print("  ⚠️  Column auto-detection may not be working")
        else:
            print("  ✅ Column auto-detection is active")
        
        # Step 3: Confirm column selection (simulate user selecting columns)
        column_data = dict(TEST_DATA)
        column_data.update({
            'filename': 'GSE90122_TO90137.tsv',
            'id_column': 'Gene_Symbol',  # Common column name in demo files
            'fc_column': 'log2FoldChange',  # Common column name
            'pval_column': 'padj'  # Common column name
        })
        
        print("  2. Confirming column selection...")
        response = session.post(f"{BASE_URL}/preview#volcano-section", data=column_data)
        
        if response.status_code != 200:
            print(f"  ❌ Column confirmation failed: {response.status_code}")
            return False
        
        print("  ✅ Column selection confirmed")
        
        # Step 4: Run analysis (simplified)
        analysis_data = dict(column_data)
        analysis_data.update({
            'aop_selection': 'AOP:1',  # Use the PXR steatosis AOP
            'logfc_threshold': '1.0'
        })
        
        print("  3. Running enrichment analysis...")
        response = session.post(f"{BASE_URL}/analyze", data=analysis_data)
        
        if response.status_code != 200:
            print(f"  ❌ Analysis failed: {response.status_code} - {response.text[:200]}")
            return False
        
        print("  ✅ Analysis completed successfully")
        
        # Step 5: Test report generation
        print("  4. Testing report generation...")
        report_data = {
            'format': 'html',
            'filename': 'GSE90122_TO90137.tsv',
            'gene_count': '1000',
            'significant_genes': '100',
            'aop_id': 'AOP:1',
            'aop_label': 'PXR activation leading to liver steatosis',
            'logfc_threshold': '1.0',
            'pval_cutoff': '0.05',
            'id_column': 'Gene_Symbol',
            'fc_column': 'log2FoldChange',
            'pval_column': 'padj',
            'id_type': 'HGNC',
            'enrichment_results': '[]'
        }
        
        response = session.post(f"{BASE_URL}/generate_report", data=report_data)
        
        if response.status_code != 200:
            print(f"  ❌ Report generation failed: {response.status_code} - {response.text[:200]}")
            return False
        
        # Check if response is HTML
        if 'Content-Type' in response.headers and 'html' in response.headers['Content-Type']:
            print("  ✅ HTML report generated successfully")
        else:
            print("  ⚠️  Report generated but content type unexpected")
        
        return True
        
    except requests.exceptions.ConnectionError:
        print("  ❌ Cannot connect to Flask application. Is it running on localhost:5000?")
        return False
    except Exception as e:
        print(f"  ❌ Unexpected error: {e}")
        return False


def test_database_operations():
    """Test database operations directly."""
    print("🗄️  Testing database operations...")
    
    try:
        # Import database components
        sys.path.append('.')
        from database import db_manager
        from config import ExperimentMetadata
        
        # Test metadata creation
        print("  1. Creating test metadata...")
        test_metadata = ExperimentMetadata(
            dataset_id="DB_TEST001",
            stressor="Database Test Chemical", 
            dosing="Test dosing",
            owner="DB Test User"
        )
        
        metadata_dict = test_metadata.to_dict()
        print(f"  ✅ Metadata created: {metadata_dict['dataset_id']}")
        
        # Test database save
        print("  2. Saving to database...")
        experiment_id = db_manager.save_experiment_metadata(
            metadata=metadata_dict,
            analysis_params={'aop_id': 'AOP:1', 'logfc_threshold': 1.0},
            results={'gene_count': 1000, 'significant_genes': 50}
        )
        
        if experiment_id:
            print(f"  ✅ Saved to database with ID: {experiment_id}")
        else:
            print("  ❌ Database save failed")
            return False
        
        # Test retrieval
        print("  3. Retrieving from database...")
        retrieved = db_manager.get_experiment(experiment_id)
        
        if retrieved and retrieved['dataset_id'] == "DB_TEST001":
            print("  ✅ Successfully retrieved experiment data")
        else:
            print("  ❌ Database retrieval failed")
            return False
        
        # Test search
        print("  4. Testing database search...")
        search_results = db_manager.search_experiments(dataset_id="DB_TEST")
        
        if len(search_results) > 0:
            print(f"  ✅ Search found {len(search_results)} experiments")
        else:
            print("  ⚠️  Search returned no results")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Database test failed: {e}")
        return False


def test_report_generation():
    """Test report generation service directly."""
    print("📊 Testing report generation service...")
    
    try:
        sys.path.append('.')
        from services.report_service import report_generator, ReportData, get_software_versions
        
        # Create test report data
        print("  1. Creating test report data...")
        test_report_data = ReportData(
            metadata={
                'dataset_id': 'REPORT_TEST001',
                'stressor': 'Test Chemical',
                'dosing': '10 µM for 24h',
                'owner': 'Report Test User',
                'upload_timestamp': datetime.now().isoformat()
            },
            filename='test_data.csv',
            gene_count=1000,
            significant_genes=150,
            aop_id='AOP:1',
            aop_label='Test AOP',
            logfc_threshold=1.0,
            pval_cutoff=0.05,
            id_column='Gene_Symbol',
            fc_column='log2FC',
            pval_column='pvalue',
            id_type='HGNC',
            enrichment_results=[
                {'KE_ID': 'KE:115', 'KE_title': 'Test KE', 'pvalue': 0.001, 'pvalue_adjusted': 0.01}
            ],
            software_versions=get_software_versions()
        )
        
        print("  ✅ Test data created")
        
        # Test HTML report generation
        print("  2. Generating HTML report...")
        html_report = report_generator.generate_html_report(test_report_data)
        
        if len(html_report) > 1000 and '<html' in html_report and 'REPORT_TEST001' in html_report:
            print("  ✅ HTML report generated successfully")
        else:
            print("  ❌ HTML report generation failed or incomplete")
            return False
        
        # Test software versions
        print("  3. Checking software versions...")
        versions = get_software_versions()
        if len(versions) > 0:
            print(f"  ✅ Software versions collected: {list(versions.keys())}")
        else:
            print("  ⚠️  No software versions found")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Report generation test failed: {e}")
        return False


def main():
    """Run all integration tests."""
    print("🚀 Starting integration tests for new functionalities...\n")
    
    tests = [
        ("Database Operations", test_database_operations),
        ("Report Generation Service", test_report_generation),
        ("Complete Metadata Workflow", test_metadata_workflow),
    ]
    
    results = {}
    
    for test_name, test_func in tests:
        print(f"\n{'='*50}")
        print(f"Running: {test_name}")
        print('='*50)
        
        try:
            results[test_name] = test_func()
        except Exception as e:
            print(f"❌ Test {test_name} crashed: {e}")
            results[test_name] = False
    
    # Summary
    print(f"\n{'='*50}")
    print("TEST RESULTS SUMMARY")
    print('='*50)
    
    passed = sum(results.values())
    total = len(results)
    
    for test_name, result in results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test_name}: {status}")
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! New functionalities are working correctly.")
    else:
        print("⚠️  Some tests failed. Check the output above for details.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)