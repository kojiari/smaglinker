MetaSAG
------
MetaSAG is bioinformatic tool that can integrate single amplified genome (SAG) and metagenome to reconstruct qualified microbial genomes.

Table of Contents
------

  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Example](#example)
  * [License](#license)

Requirements
------
MetaSAG is written by python3. It requires the following to run:

  * python (>= 3.7.6)
  * BWA (>= 0.7.17)
  * Samtools (>= 1.9)
  * Prokka (>= 1.14.6)
  * CheckM (>= 1.1.2)
  * HaploMerger2 (>= 20180603)

In addition, MetaSAG requires the following Python libraries:

  * pysam (>= 0.15.4)

Installation
------

* Download the packages from https://github.com/kojiari/metasag and grant permission to run MetaSAG.
```
git clone https://github.com/kojiari/metasag
cd metasag
chmod a+x metasag.py

```

* Add path of HaploMerger2 to your environment.
```
export HAPLOMERGER2_PATH=/path/to/haplomerger2
```

Usage
------
### 1. Binning
```
// example
metasag.py binning -o sgbin/ example/mock-ma.fasta example/mock-sag*.fasta
```

Metagenome assembly (e.g. example/mock-ma.fasta) is classified to single-cell genome-guided bin (sgBin) by mapping on reference SAGs (e.g. example/mock-sag*.fasta).

### 2. Merge
```
// example
metasag.py merge -o merged_genome/mock-sag01 example/mock-sag01.fasta sgbin/mock-sag01.fasta
metasag.py merge -o merged_genome/mock-sag02 example/mock-sag02.fasta sgbin/mock-sag02.fasta
metasag.py merge -o merged_genome/mock-sag02 example/mock-sag03.fasta sgbin/mock-sag03.fasta
```

Paired SAGs (e.g. example/mock-sag01.fasta) and sgBins (e.g. sgbin/mock-sag01.fasta) are merged to single-cell genome-guided MAG (sgMAG) or metagenome-guided SAG (mgSAG).

Output
------
### 1. Binning

* classification
* *.fasta
* unclassified.fasta

### 2. Merge

* sgmag.fasta or mgsag.fasta
+ checkm.result

Example data
------
Examples are the mock microbial community (Cell-Mock-001) obtained from the National Institute of Technology and Evaluation Biological Resource Center, Japan.  This mock microbial community was composed of 15 bacterial strains detected in various environments (intestinal, oral, skin, and natural environment).

|Accession|contigs|total_length|chrom|plasmid|asm_level|NBRC|Name|
|----|----|----|----|----|----|----|----|
|GCF_006828885.1|5|2,277,780|5|0|Contig|113353|Bifidobacterium pseudocatenulatum Scardovi et al. 1979|
|GCF_006538485.1|40|3,113,921|40|0|Contig|15291|Corynebacterium striatum (Chester 1901) Eberson 1918|
|GCF_006739385.1|1|2,494,738|1|0|Complete|107605|Cutibacterium acnes subsp. acnes (Gilchrist 1900) Nouioui et al. 2018|
|GCF_006742345.1|5|4,989,532|1|4|Complete|113350|Bacteroides uniformis Eggerth and Gagnon 1933|
|GCF_006739545.1|1|5,179,960|1|0|Complete|113806|Parabacteroides distasonis Sakamoto and Benno 2006|
|GCF_006741845.1|2|4,295,305|1|1|Complete|13719|Bacillus subtilis subsp. subtilis (Ehrenberg 1835) Cohn 1872|
|GCF_006742205.1|2|2,427,041|1|1|Complete|100911|Staphylococcus epidermidis (Winslow and Winslow 1908) Evans 1916|
|GCF_006740305.1|1|1,910,306|1|0|Complete|3202|Lactobacillus delbrueckii subsp. delbrueckii (Leichmann 1896) Beijerinck 1901|
|GCF_006739205.1|1|2,018,796|1|0|Complete|13955|Streptococcus mutans Clarke 1924|
|GCF_006742065.1|4|4,705,096|1|3|Complete|13949|Clostridium butyricum Prazmowski 1880|
|GCF_006538465.1|2|5,687,315|2|0|Contig|113352|Clostridium clostridioforme corrig. (Burri and Ankersmit 1906) Kaneuchi et al. 1976|
|GCF_001515545.1|23|4,629,566|23|0|Contig|13299|Comamonas terrigena (ex Hugh 1962) De Vos et al. 1985 emend. Wauters et al. 2003|
|GCF_010509415.1|2|4,755,096|1|1|Complete|3301|Escherichia coli (Migula 1895) Castellani and Chalmers 1919|
|GCF_006757745.1|9|3,433,938|1|8|Complete|102413|Acinetobacter radioresistens Nishimura et al. 1988|
|GCF_000412675.1|1|6,156,701|1|0|Complete|14164|Pseudomonas putida (Trevisan 1889) Migula 1895|

* example/mock-ma.fasta : metagenome assembly of cell mock community
* example/mock-sag01~15.fasta : SAGs of cell mock community

License
------
MetaSAG is licensed under the  LGPL-2.1 License.
