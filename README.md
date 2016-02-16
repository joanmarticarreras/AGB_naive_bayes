# AGB_naive_bayes
Project of AGB class: Naive Bayes model to classify tumor types gene expression patterns

INFO AT: http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/index.html

## TO DO
- [x] Are there prior probabilities to take into account? No, all 1/8.
- [x] Get the sets and parse them.
- [ ] Define and create the training sets and the test set (equal n).
- [x] Z-scores to descrete values.
- [ ] Measure Likelyhoods.
- [ ] Measure Mutual information.
- [x] Do we need pseudocounts? YES

### Get the sets
```bash
mkdir sources

wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/brca_gene_zscore_full-filtered.txt -O ./sources/brca.tbl && \
wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/coad_gene_zscore_full-filtered.txt -O ./sources/coad.tbl && \
wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/hnsc_gene_zscore_full-filtered.txt -O ./sources/hnsc.tbl && \
wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/kirc_gene_zscore_full-filtered.txt -O ./sources/kric.tbl && \
wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/luad_gene_zscore_full-filtered.txt -O ./sources/luad.tbl && \
wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/lusc_gene_zscore_full-filtered.txt -O ./sources/lusc.tbl && \
wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/prad_gene_zscore_full-filtered.txt -O ./sources/prad.tbl && \
wget http://regulatorygenomics.upf.edu/courses/Master_AGB/Exercise_NaiveBayes/thca_gene_zscore_full-filtered.txt -O ./sources/thca.tbl

```
### Look for funy values with AWK
```bash
gawk '{print $3}' coad.tbl | sort | uniq -c | more

gawk '{print $3}' coad.tbl | sort | uniq -c | tail

etc


```
### Filter genes with NA or Inf in at least 1 samples
```bash
cat sources/* | egrep '\bNA\b|\bInf\b' | gawk '{print $1 }' | sort | uniq > not_valid_genes.txt

```
Load them into a hash, and just check in NoSoNaive.

We need to know the minumum number of samples in each cancer type, because it will be the number of samples used in our study.

```bash
gawk '{print FILENAME, NF}' f_* | sort | uniq | more

```
The number is 288.


### Do parser.pl
- [x] Filter train & testing (144 each), where 288 is the minumum number of samples from the cancer type files.
    - WARNING: 144 --> 147 when applying pseudocounts.

- [x] The filter was pseudo-randomized (odds -> test, even -> train). Better against stratification than 0..143 vs 144..end.
- [x] Create a bunch of test and train files.

### Do NotSoNaive.pl
- [x] Filter which genes have at leat 1 value as NA or Inf, and store it within a hash.
    - Filter from the test & training set those genes.
- [x] Calcule Liklihoods.
- [x] Calcule Conditional entropies.
- [x] Calcule entropies.
- [x] Information gain.

### One command to rule them all

```bash
perl NotSoNaive.pl -train sources/train_brca.tbl,sources/train_thca.tbl,sources/train_coad.tbl,sources/train_hnsc.tbl,sources/train_kric.tbl,sources/train_luad.tbl,sources/train_lusc.tbl,sources/train_prad.tbl -test sources/test_brca.tbl,sources/test_thca.tbl,sources/test_coad.tbl,sources/test_hnsc.tbl,sources/test_kric.tbl,sources/test_luad.tbl,sources/test_lusc.tbl,sources/test_prad.tbl -not not_valid_genes.txt | egrep "IG:" > information_gain.tbl
```

### Information gain
* 3rd quantile = 0.241287
* Number of genes with IG > 3rd quantile = 4755
