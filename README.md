![Title](https://github.com/debaleena-bhowmik/HUGMi/blob/main/hugmi_logo_nb.png)

Human Uro-Genital Microbiome database (HUGMi) and hybrid classifier based on q2-plugins, for 16S rRNA amplicon sequencing data.

## Summary:
The HUGMi workflow incorporates a niche-specific 16S rRNA database with an innovative Hybrid Classification algorithm, enabling enhanced bacterial species resolution from amplicon sequencing data.

HUGMi represents the sole reference database specifically optimized for bacterial taxa inhabiting human urogenital environments. Such ecological niche-focused resources minimize false positives while providing phylogenetically coherent reference structures, thereby enhancing taxonomic precision. The database maintains standardized prokaryotic nomenclature across all seven taxonomic hierarchical levels in accordance with International Committee on Systematics of Prokaryotes (ICSP) protocols.

The Hybrid Classifier algorithm unites the QIIME2 BLAST and sklearn-based classification methods to facilitate superior species-level taxonomic assignment. This classification system offers parameterization flexibility through adjustable confidence thresholds and is compatible with any 16S rRNA databases.


### Requirements:

To run the program create a QIIME 2 (q2) conda environment and execute the program from within that environment. The program is compatible for QIIME 2 amplicon distribution version 2023.9 and above. The given (.qza) files were created in q2 version 2024.10.

Note - Make sure all the .qza files are generated in the same q2 version to prevent version conflict.



## Syntax:

### General usage (with default parameters):
```
python hugmi_hybrid_classifier.py --rep-seqs test-rep-seqs.qza --database-taxonomy HUGMi_taxa.qza --database-sequences HUGMi_seq.qza --classifier HUGMi_v4_classifier.qza
```


### For help on how to run the program:
```
python hugmi_hybrid_classifier.py --help

usage: hugmi_hybrid_classifier.py [-h] --rep-seqs REP_SEQS --database-taxonomy DATABASE_TAXONOMY --database-sequences DATABASE_SEQUENCES --classifier
                                  CLASSIFIER [--max_accepts MAX_ACCEPTS] [--perc_identity PERC_IDENTITY] [--query_cov QUERY_COV]
                                  [--classifier_confidence CLASSIFIER_CONFIDENCE] [--num_threads NUM_THREADS] [--output-dir OUTPUT_DIR]

Hybrid taxonomy classifier using QIIME 2

-h, --help                      show this help message and exit
--rep-seqs REP_SEQS             Path to representative sequences artifact (.qza)
--database-taxonomy             DATABASE_TAXONOMY
                                Path to database taxonomy file (.qza)
--database-sequences            DATABASE_SEQUENCES
                                Path to database sequences file (.qza)
--classifier CLASSIFIER         Path to pre-trained sklearn classifier (.qza)
--max_accepts MAX_ACCEPTS       Maximum number of hits to keep for each query (default 10)
--perc_identity PERC_IDENTITY   Rejects if percent identity to query is lower (default 1.0)
--query_cov QUERY_COV           Rejects if alignment coverage is lower (default 0.95)
--classifier_confidence         CLASSIFIER_CONFIDENCE
                                Confidence threshold for sklearn (default 0.7)
--num_threads NUM_THREADS       Number of threads (CPUs) used for both classifiers (default 5)
--output-dir OUTPUT_DIR         Directory to save outputs
```


### Usage of the hybrid classifier with other 16S rRNA databases:
```
python hugmi_hybrid_classifier.py --rep-seqs rep-seqs.qza --database-taxonomy ~/path/to/database-taxonomy-qza --database-sequences ~/path/to/database-sequence-qza --classifier ~/path/to/16S-region-soecific-classifier-qza
```


Information regarding how to train 16S region specific classier can be found [here](https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494#heading--sixth-header).

### File description:

hugmi_hybrid_classifier.py:: The hybrid classfier program

HUGMi.zip:: Contains HUGMi database sequence file (HUGMi_seq.fasta) & taxonomy file (HUGMi_taxa.txt)

qza_files.zip:: Contains the following q2 compatible files:

i) HUGMi_seq.qza - HUGMi sequence (.qza) file
                
ii) HUGMi_taxa.qza - HUGMi taxonomy (.qza) file
                
iii) HUGMi_v3v4_classifier.qza - HUGMi sklearn trained classifier for V3V4 region
                
iv) HUGMi_v4_classifier.qza - HUGMi sklearn trained classifier for V4 region
                
v) test-rep-seq.qza - Sample representative sequence file
                
                

Citation [pre-print](https://doi.org/10.1101/2025.05.01.651608) :

Bhowmik, D., & Paul, S. (2025). HUGMi: Human Uro-Genital Microbiome database and hybrid classifier for improved species level annotation of 16S rRNA amplicon sequences. bioRxiv, 2025-05.
