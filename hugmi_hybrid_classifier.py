
""" q2-BLAST-sklearn-hybrid-classifier

Note: sklearn classifier needs to be provided separately.
 
"""

# import packages
import os
import pandas as pd
from qiime2 import Artifact, Metadata
from qiime2.plugins.feature_classifier.methods import makeblastdb
from qiime2.plugins.feature_classifier.pipelines import classify_consensus_blast
from qiime2.plugins.feature_classifier.methods import classify_sklearn
from qiime2.plugins.feature_table.methods import filter_seqs

def create_hybrid_classifier(
    rep_seqs_path,
    database_taxonomy_path,
    database_sequences_path,
    classifier_path,
    max_accepts,
    perc_identity,
    query_cov,
    confidence,
    num_threads,
    output_dir="hybrid_classification_results"   
    
):
    """
    Create a hybrid classifier using BLAST and pre-trained sklearn classifier.
    
    Parameters:
    -----------
    rep_seqs_path : str
        Path to representative sequences artifact (.qza) from DADA2
    database_taxonomy_path : str
        Path to database taxonomy file (.qza)
    database_sequences_path : str
        Path to database sequences file (.qza)
    classifier_path : str
        Path to pre-trained sklearn classifier (.qza)
    output_dir : str
        Directory to save outputs
    """
    print("Loading input files...")
    rep_seqs = Artifact.load(rep_seqs_path)
    database_taxonomy = Artifact.load(database_taxonomy_path)
    database_sequences = Artifact.load(database_sequences_path)
    classifier = Artifact.load(classifier_path)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Create BLAST database
    print("Creating BLAST database...")
    blast_db, = makeblastdb(
        sequences=database_sequences
    )
    
    # Step 2: Run BLAST classification
    print("Running BLAST classification...")
    blast_results = classify_consensus_blast(
        query=rep_seqs,
        blastdb=blast_db,
        reference_taxonomy=database_taxonomy,
        maxaccepts=max_accepts,
        perc_identity=perc_identity,
        query_cov=query_cov,
        num_threads=num_threads
    )
    
    # BLAST results
    blast_taxonomy_metadata = blast_results.classification.view(Metadata)
    blast_taxonomy_df = blast_taxonomy_metadata.to_dataframe()
    blast_taxonomy_df = blast_taxonomy_df.reset_index()
    
    # Step 3: Filter ASVs with insufficient taxonomic classification
    def is_low_resolution(taxon):
        """Check if taxonomy is less than genus level or Unassigned"""
        if taxon == 'Unassigned' or pd.isna(taxon):
            return True
        
        # Count the number of levels in the taxonomy
        levels = taxon.split(';')
        # Check if we have at least 6 levels (kingdom to genus)
        if len(levels) < 6:
            return True
        
        # Check if the genus level exists and is not empty
        genus_level = levels[5] if len(levels) > 5 else ""
        if not genus_level or genus_level.endswith('__'):
            return True
            
        return False
    
    # Mark ASVs with low-resolution taxonomy
    blast_taxonomy_df['LowResolution'] = blast_taxonomy_df['Taxon'].apply(is_low_resolution)
    
    # Create a dataframe with results that have sufficient resolution
    high_res_df = blast_taxonomy_df[~blast_taxonomy_df['LowResolution']].copy()
    high_res_df['Method'] = 'BLAST'
    
    # Create list of ASVs that need better classification
    low_res_asv_ids = blast_taxonomy_df[blast_taxonomy_df['LowResolution']]['Feature ID'].tolist()
    
    if low_res_asv_ids:
        print(f"Filtering {len(low_res_asv_ids)} ASVs for sklearn classification...")
        
        low_res_metadata = pd.DataFrame(index=low_res_asv_ids)
        low_res_metadata.index.name = 'feature id'
        low_res_metadata_qiime = Metadata(low_res_metadata)
        
        print("Extracting sequences for low-resolution ASVs...")
        filtered_seqs_result = filter_seqs(
            data=rep_seqs,
            metadata=low_res_metadata_qiime
        )
        
        # Step 5: Run sklearn classifier on the filtered sequences
        print("Running sklearn classification on low-resolution ASVs...")
        sklearn_results = classify_sklearn(
            reads=filtered_seqs_result.filtered_data,
            classifier=classifier,
            confidence=confidence
        )
        
        # Convert sklearn result to dataframe
        sklearn_taxonomy_metadata = sklearn_results.classification.view(Metadata)
        sklearn_taxonomy_df = sklearn_taxonomy_metadata.to_dataframe()
        
        # Reset the index to make the feature IDs a column
        sklearn_taxonomy_df = sklearn_taxonomy_df.reset_index()
        
        # Rename columns to ensure consistency
        if 'index' in sklearn_taxonomy_df.columns:
            sklearn_taxonomy_df = sklearn_taxonomy_df.rename(columns={'index': 'Feature ID'})
        
        if 'Taxon' not in sklearn_taxonomy_df.columns and 'Taxonomy' in sklearn_taxonomy_df.columns:
            sklearn_taxonomy_df = sklearn_taxonomy_df.rename(columns={'Taxonomy': 'Taxon'})
               
        # Add method column
        sklearn_taxonomy_df['Method'] = 'sklearn'        
        
        # Rename the Confidence column to match our desired output format
        if 'Confidence' in sklearn_taxonomy_df.columns:
            sklearn_taxonomy_df = sklearn_taxonomy_df.rename(columns={'Confidence': 'Consensus/Confidence'})
        else:
            print("No low-resolution ASVs to process with sklearn")
            sklearn_taxonomy_df = pd.DataFrame(columns=['Feature ID', 'Taxon', 'Method', 'Consensus/Confidence'])

    # Step 6: Combine results from both methods
    
    # Renaming high_res_df column
    high_res_df.rename({"Consensus":"Consensus/Confidence"}, axis=1, inplace=True)
    
    print("Combining results from both methods...")
    
    # Ensure all required columns exist in both dataframes
    required_columns = ['Feature ID', 'Taxon', 'Method', 'Consensus/Confidence']
    
    for df in [high_res_df, sklearn_taxonomy_df]:
        for col in required_columns:
            if col not in df.columns:
                df[col] = None
    
    combined_df = pd.concat([high_res_df[required_columns], 
                           sklearn_taxonomy_df[required_columns]])
    
    # Saving combined results
    combined_output_path = os.path.join(output_dir, "hybrid_taxonomy.tsv")
    combined_df.to_csv(combined_output_path, sep='\t', index=False)
    
    print(f"Hybrid classification complete. Results saved to: {combined_output_path}")
    return combined_df

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Hybrid taxonomy classifier using QIIME 2")
    parser.add_argument("--rep-seqs", required=True, help="Path to representative sequences artifact (.qza)")
    parser.add_argument("--database-taxonomy", required=True, help="Path to database taxonomy file (.qza)")
    parser.add_argument("--database-sequences", required=True, help="Path to database sequences file (.qza)")
    parser.add_argument("--classifier", required=True, help="Path to pre-trained sklearn classifier (.qza)")
    parser.add_argument("--max_accepts", type=int, default=10, help="Maximum number of hits to keep for each query (default 10")
    parser.add_argument("--perc_identity", type=float, default=1.0, help="Rejects if percent identity to query is lower (default 1.0)")
    parser.add_argument("--query_cov", type=float, default=0.95, help="Rejects if alignment coverage is lower (default 0.95")
    parser.add_argument("--classifier_confidence", type=float, default=0.7, help="Confidence threshold for sklearn (default 0.7)")
    parser.add_argument("--num_threads", type=int, default=5, help="Number of threads (CPUs) used for both classifiers (default 5)")
    parser.add_argument("--output-dir", default="hybrid_classification_results", help="Directory to save outputs")
    
    
    args = parser.parse_args()
    
    create_hybrid_classifier(
        args.rep_seqs,
        args.database_taxonomy,
        args.database_sequences,
        args.classifier,
        args.max_accepts,
        args.perc_identity,
        args.query_cov,
        args.classifier_confidence,
        args.num_threads,
        args.output_dir
    )
