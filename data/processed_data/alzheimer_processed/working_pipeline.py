"""
ALZHEIMER PREPROCESSING PIPELINE
It handles metadata-only files (GSE125050, GSE67333), multiple platforms from GSE84422 (GPL570, GPL96, GPL97)
"""
import pandas as pd
import numpy as np
import gzip 
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from scipy import stats
import warnings
import json
warnings.filterwarnings('ignore')

print("ALZHEIMER DATA PREPROCESSING PIPELINE")

class AlzheimerPreprocessor:
    """
    Working preprocessor that handles the data structure
    """
    def __init__(self):
        self.datasets_info = {
            "GSE125050": {
                "description": "ROSMAP - Alzheimer's gene expression by cell type",
                "platform": "metadata_only",
                "expected_size": "small",
                "priority": 4,  # Low priority (lack of expression data)
            },
            "GSE84422-GPL570": {
                "description": "MSBB - Affymetrix Human Genome U133 Plus 2.0",
                "platform": "microarray_GPL570", 
                "expected_size": "large",
                "priority": 1,  # High priority
            },
            "GSE84422-GPL96": {
                "description": "MSBB - Affymetrix GeneChip Human Gene 1.0 ST",
                "platform": "microarray_GPL96",
                "expected_size": "large", 
                "priority": 2, 
            },
            "GSE84422-GPL97": {
                "description": "MSBB - Affymetrix GeneChip Human Gene 1.0 ST",
                "platform": "microarray_GPL97",
                "expected_size": "large",
                "priority": 3,  
            },
            "GSE67333": {
                "description": "Transcriptomics profiling of Alzheimer disease",
                "platform": "metadata_only",
                "expected_size": "small", 
                "priority": 5,  # Low priority (lack of expression data)
            }
        }
        
        self.processed_data = {}
        print("Working Alzheimer Preprocessor initialized")
        print(f"Targeting {len(self.datasets_info)} known datasets")

    def detect_files(self, data_directory):
        """
        Detect files based on your actual file structure
        """
        print("\n DETECTING FILES IN YOUR DATA STRUCTURE")
        
        data_dir = Path(data_directory)
        detected_files = {}
        
        # Find all GEO series matrix files
        all_geo_files = list(data_dir.rglob("*_series_matrix.txt.gz"))
        
        print(f"Found {len(all_geo_files)} GEO files:")
        for f in all_geo_files:
            size_mb = f.stat().st_size / (1024*1024)
            print(f"{f.name}: {size_mb:.3f} MB")
        
        # Match files to our known datasets
        for dataset_id, info in self.datasets_info.items():
            pattern = f"{dataset_id}_series_matrix.txt.gz"
            matches = [f for f in all_geo_files if f.name == pattern]
            
            if matches:
                file_path = matches[0]
                size_mb = file_path.stat().st_size / (1024*1024)
                
                detected_files[dataset_id] = {
                    'filepath': file_path,
                    'info': info,
                    'size_mb': size_mb,
                    'status': 'found'
                }
                print(f"{dataset_id}: {size_mb:.3f} MB")
            else:
                detected_files[dataset_id] = {
                    'filepath': None,
                    'info': info,
                    'status': 'missing'
                }
                print(f"{dataset_id}: Not found")
        
        return detected_files

    def parse_geo_file(self, file_path, dataset_id):
        """
        Parse GEO file with proper handling of metadata-only files
        """
        print(f"\n PARSING {dataset_id}")
        
        try:
            # Read the file
            with gzip.open(file_path, 'rt', encoding='utf-8') as f:
                lines = f.readlines()
            
            print(f" File contains {len(lines)} lines")
            
            # Parse metadata
            metadata = {}
            sample_info = {}
            matrix_start = None
            matrix_end = None
            
            for i, line in enumerate(lines):
                line = line.strip()
                
                # Sample information
                if line.startswith('!Sample_geo_accession'):
                    sample_info['geo_accession'] = line.split('\t')[1:]
                elif line.startswith('!Sample_title'):
                    sample_info['title'] = line.split('\t')[1:]
                elif line.startswith('!Sample_source_name'):
                    sample_info['source_name'] = line.split('\t')[1:]
                elif line.startswith('!Sample_characteristics'):
                    parts = line.split('\t')
                    if len(parts) > 1:
                        char_name = parts[0].replace('!Sample_characteristics_ch1', '').strip('_')
                        if not char_name:
                            char_name = 'characteristics'
                        metadata[char_name] = parts[1:]
                
                # Matrix boundaries
                elif line.startswith('!series_matrix_table_begin'):
                    matrix_start = i + 1
                elif line.startswith('!series_matrix_table_end'):
                    matrix_end = i
                    break
            
            # Create metadata DataFrame
            metadata_df = None
            if sample_info.get('geo_accession'):
                n_samples = len(sample_info['geo_accession'])
                records = []
                
                for i in range(n_samples):
                    record = {}
                    
                    # Basic sample info
                    for key, values in sample_info.items():
                        if i < len(values):
                            record[key] = values[i].strip('"')
                    
                    # Characteristics
                    for key, values in metadata.items():
                        if i < len(values):
                            record[key] = values[i].strip('"')
                    
                    records.append(record)
                
                metadata_df = pd.DataFrame(records)
                print(f"Metadata: {len(records)} samples, {len(metadata_df.columns)} fields")
            
            # Parse expression matrix if it exists
            expression_df = None
            if matrix_start is not None and matrix_end is not None and matrix_end > matrix_start:
                matrix_lines = lines[matrix_start:matrix_end]
                
                if matrix_lines:
                    print(f"Parsing expression matrix: {len(matrix_lines)} lines")
                    
                    # Parse header
                    header = matrix_lines[0].strip().split('\t')
                    if len(header) > 1:
                        sample_names = [name.strip('"') for name in header[1:]]
                        
                        # Parse data
                        gene_ids = []
                        expression_values = []
                        
                        for line in matrix_lines[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) >= len(sample_names) + 1:
                                gene_id = parts[0].strip('"')
                                values = []
                                
                                for val_str in parts[1:len(sample_names)+1]:
                                    val_str = val_str.strip().strip('"')
                                    try:
                                        if val_str.lower() in ['null', 'na', 'nan', '', '--']:
                                            values.append(np.nan)
                                        else:
                                            values.append(float(val_str))
                                    except ValueError:
                                        values.append(np.nan)
                                
                                gene_ids.append(gene_id)
                                expression_values.append(values)
                        
                        if expression_values:
                            expression_df = pd.DataFrame(
                                expression_values,
                                index=gene_ids,
                                columns=sample_names
                            )
                            print(f"Expression data: {expression_df.shape[0]} genes × {expression_df.shape[1]} samples")
                        else:
                            print(f"No valid expression data in matrix")
                    else:
                        print(f"Invalid matrix header")
                else:
                    print(f"Empty matrix section")
            else:
                print(f"No expression matrix found (metadata-only file)")
            
            return {
                'expression': expression_df,
                'metadata': metadata_df,
                'dataset_id': dataset_id,
                'has_expression': expression_df is not None,
                'sample_count': len(metadata_df) if metadata_df is not None else 0
            }
            
        except Exception as e:
            print(f"Error parsing {dataset_id}: {e}")
            return None

    def extract_labels(self, metadata_df, dataset_id):
        """
        Extract labels from metadata
        """
        if metadata_df is None:
            return None
        
        print(f"\n EXTRACTING LABELS for {dataset_id}")
        
        labels_df = pd.DataFrame({
            'sample_id': metadata_df.get('geo_accession', metadata_df.index),
            'dataset_id': [dataset_id] * len(metadata_df)
        })
        
        # Look for diagnosis information
        diag_columns = [col for col in metadata_df.columns 
                       if any(term in col.lower() for term in 
                             ['diagnosis', 'disease', 'condition', 'status', 'group'])]
        
        print(f"Potential diagnosis columns: {diag_columns}")
        
        if diag_columns:
            primary_col = diag_columns[0]
            print(f" Using: {primary_col}")
            
            # Show unique values
            unique_vals = metadata_df[primary_col].value_counts()
            print(f"   Values found: {list(unique_vals.index)}")
            
            # Simple mapping
            labels_df['diagnosis'] = metadata_df[primary_col]
        else:
            labels_df['diagnosis'] = 'unknown'
            print("No clear diagnosis column found")
        
        return labels_df

    def process_expression_data(self, expression_df, dataset_id):
        """
        Process and normalize expression data
        """
        if expression_df is None:
            return None
        
        print(f"\n PROCESSING EXPRESSION DATA for {dataset_id}")
        
        # Basic QC
        missing_pct = (expression_df.isna().sum().sum() / (expression_df.shape[0] * expression_df.shape[1])) * 100
        print(f"Missing values: {missing_pct:.2f}%")
        
        # Handle missing values
        if missing_pct > 0:
            imputer = SimpleImputer(strategy='median')
            expression_imputed = pd.DataFrame(
                imputer.fit_transform(expression_df.T).T,
                index=expression_df.index,
                columns=expression_df.columns
            )
            print(f"Imputed missing values")
        else:
            expression_imputed = expression_df.copy()
        
        # Log transform for microarray data (if not already log-transformed)
        platform = self.datasets_info[dataset_id]['platform']
        if 'microarray' in platform:
            # Check if data needs log transformation
            max_val = expression_imputed.max().max()
            if max_val > 50:  # Likely not log-transformed
                expression_log = np.log2(expression_imputed + 1)
                print(f"Applied log2 transformation")
            else:
                expression_log = expression_imputed
                print(f"Data appears already log-transformed")
        else:
            expression_log = expression_imputed
        
        # Z-score normalization
        scaler = StandardScaler()
        expression_normalized = pd.DataFrame(
            scaler.fit_transform(expression_log.T).T,
            index=expression_log.index,
            columns=expression_log.columns
        )
        
        print(f"Normalization complete: {expression_normalized.shape}")
        
        return expression_normalized

    def save_results(self, dataset_id, expression_df, metadata_df, labels_df):
        """
        Save processed results
        """
        print(f"\n SAVING RESULTS for {dataset_id}")
        
        output_dir = Path("results") / dataset_id
        output_dir.mkdir(parents=True, exist_ok=True)
        
        saved_files = []
        
        if expression_df is not None:
            expr_file = output_dir / "expression_normalized.csv"
            expression_df.to_csv(expr_file)
            saved_files.append(f"Expression: {expr_file.name}")
            print(f"Expression saved: {expression_df.shape}")
        
        if metadata_df is not None:
            meta_file = output_dir / "metadata.csv"
            metadata_df.to_csv(meta_file, index=False)
            saved_files.append(f"Metadata: {meta_file.name}")
            print(f"Metadata saved: {metadata_df.shape}")
        
        if labels_df is not None:
            labels_file = output_dir / "labels.csv"
            labels_df.to_csv(labels_file, index=False)
            saved_files.append(f"Labels: {labels_file.name}")
            print(f"Labels saved: {labels_df.shape}")
        
        # Save summary
        summary = {
            'dataset_id': dataset_id,
            'has_expression': expression_df is not None,
            'n_genes': expression_df.shape[0] if expression_df is not None else 0,
            'n_samples': len(metadata_df) if metadata_df is not None else 0,
            'platform': self.datasets_info[dataset_id]['platform'],
            'files_saved': saved_files
        }
        
        summary_file = output_dir / "summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"Summary saved: {summary_file}")
        
        return summary

def run_working_pipeline(data_directory):
    """
    Run the complete working pipeline
    """
    print("RUNNING WORKING ALZHEIMER PIPELINE")
    
    # Initialize
    preprocessor = AlzheimerPreprocessor()
    
    # Detect files
    print("\n STEP 1: FILE DETECTION")
    detected_files = preprocessor.detect_files(data_directory)
    
    # Process each dataset
    print("\n STEP 2: PROCESSING DATASETS")
    results = {}
    
    # Sort by priority (highest first)
    sorted_datasets = sorted(detected_files.items(), 
                           key=lambda x: x[1]['info'].get('priority', 999))
    
    for dataset_id, file_info in sorted_datasets:
        if file_info['status'] == 'found':
            print(f"PROCESSING {dataset_id}")
            print(f"Size: {file_info['size_mb']:.3f} MB")
            print(f"Platform: {file_info['info']['platform']}")
            
            try:
                # Parse file
                parsed_data = preprocessor.parse_geo_file(
                    file_info['filepath'], dataset_id
                )
                
                if parsed_data:
                    # Extract labels
                    labels_df = preprocessor.extract_labels(
                        parsed_data['metadata'], dataset_id
                    )
                    
                    # Process expression data (if exists)
                    if parsed_data['has_expression']:
                        normalized_expr = preprocessor.process_expression_data(
                            parsed_data['expression'], dataset_id
                        )
                    else:
                        normalized_expr = None
                        print("   No expression data to process")
                    
                    # Save results
                    summary = preprocessor.save_results(
                        dataset_id, normalized_expr, parsed_data['metadata'], labels_df
                    )
                    
                    results[dataset_id] = {
                        'status': 'success',
                        'summary': summary,
                        'has_expression': parsed_data['has_expression']
                    }
                    
                    print(f" {dataset_id} COMPLETED SUCCESSFULLY")
                    
                else:
                    results[dataset_id] = {'status': 'parsing_failed'}
                    print(f" {dataset_id} parsing failed")
                    
            except Exception as e:
                results[dataset_id] = {'status': 'error', 'error': str(e)}
                print(f" {dataset_id} error: {e}")
        
        else:
            results[dataset_id] = {'status': 'file_not_found'}
            print(f" {dataset_id}: File not found")
    
    # Final summary
    print("PIPELINE SUMMARY")
    
    successful = [k for k, v in results.items() if v['status'] == 'success']
    with_expression = [k for k, v in results.items() 
                      if v.get('has_expression', False)]
    
    print(f"\n RESULTS:")
    print(f"Total datasets processed: {len(results)}")
    print(f"Successful: {len(successful)}")
    print(f"With expression data: {len(with_expression)}")
    
    print(f"\n DATASET STATUS:")
    for dataset_id, result in results.items():
        status = result['status']
        if status == 'success':
            has_expr = " Expression" if result['has_expression'] else " Metadata only"
            print(f"    {dataset_id}: {has_expr}")
        else:
            print(f"    {dataset_id}: {status}")
    
    if with_expression:
        print(f"\n SUCCESS! {len(with_expression)} datasets with expression data ready for analysis")
        print(f"Results saved in 'results/' directory")
        print(f"Ready for feature engineering and machine learning")
    else:
        print(f"\n No datasets with expression data found")
        print(f"Consider downloading full expression datasets")
    
    return results

# Main execution
if __name__ == "__main__":
    # Your data directory
    data_dir = "C:/Users/diego/Projects2025/neurocardiac_project/data/raw/alzheimer_data"
    
    print("Starting pipeline with your actual data...")
    
    # Run the working pipeline
    results = run_working_pipeline(data_dir)
    
    print(f"\n Pipeline completed!")
    print(f"Check the 'results/' directory for processed data.")