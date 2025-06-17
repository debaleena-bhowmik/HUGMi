![Title](https://github.com/debaleena-bhowmik/HUGMi/blob/main/hugmi_logo_nb.png)

Human Uro-Genital Microbiome database (HUGMi) and hybrid classifier based on q2-plugins, for 16S rRNA amplicon sequencing data.

## Summary:
The HUGMi workflow incorporates a niche-specific 16S rRNA database with an innovative Hybrid Classification algorithm, enabling enhanced bacterial species resolution from amplicon sequencing data. < br />
HUGMi represents the sole reference database specifically optimized for bacterial taxa inhabiting human urogenital environments. Such ecological niche-focused resources minimize false positives while providing phylogenetically coherent reference structures, thereby enhancing taxonomic precision. The database maintains standardized prokaryotic nomenclature across all seven taxonomic hierarchical levels in accordance with International Committee on Systematics of Prokaryotes (ICSP) protocols.
The Hybrid Classifier algorithm unites the QIIME2 BLAST and sklearn-based classification methods to facilitate superior species-level taxonomic assignment. This classification system offers parameterization flexibility through adjustable confidence thresholds and is compatible with any 16S rRNA databases.



## Syntax:
```
python hugmi_hybrid_classifier.py --rep-seqs representative-sequence.qza --database-taxonomy hugmi_taxonomy.qza --database-sequences hugmi_sequences.qza --classifier hugmi_classifier.qza
```



Citation [pre-print](https://doi.org/10.1101/2025.05.01.651608) :

Bhowmik, D., & Paul, S. (2025). HUGMi: Human Uro-Genital Microbiome database and hybrid classifier for improved species level annotation of 16S rRNA amplicon sequences. bioRxiv, 2025-05.
