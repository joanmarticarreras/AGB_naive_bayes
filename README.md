# AGB_naive_bayes
Project of AGB class: Naive Bayes model to classify tumor types gene expression patterns

INFO AT: http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/index.html

## TO DO
- [ ] Are there prior probabilities to take into account?
- [ ] Get the sets and parse them.
- [ ] Define and create the training sets and the test set (equal n).
- [ ] Z-scores to descrete values.
- [ ] Measure Likelyhoods.
- [ ] Do we need pseudocounts?

### Get the sets
```bash

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/brca_gene_zscore_full-filtered.txt -O brca.tbl

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/coad_gene_zscore_full-filtered.txt -O coad.tbl

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/hnsc_gene_zscore_full-filtered.txt -O hnsc.tbl

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/kirc_gene_zscore_full-filtered.txt -O kric.tbl

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/luad_gene_zscore_full-filtered.txt -O luad.tbl

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/lusc_gene_zscore_full-filtered.txt -O lusc.tbl

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/prad_gene_zscore_full-filtered.txt -O prad.tbl

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/thca_gene_zscore_full-filtered.txt -O thca.tbl

```	