# cnnLSV
cnnLSV: detecting structural variants by encoding long-read alignment information and convolutional neural network

<h1>
<p align="center">
  <img src="logo.png" alt="Logo" width="600" height="150">
</h1>

---
## Installation
### Clone
	git clone https://github.com/mhuidong/cnnLSV.git
### Requirements
python3, cv2, numpy, torch, torchvision, pysam, cigar

---
## Usage
### Detecting SVs
	python cnnLSV.py <sorted.bam> <initial.vcf> <filtered.vcf>
### [OPTIONAL] Removing Redundant information
#### CnnLSV outputs the callset merged of several callers. This means that one SV may be detected by several callers. You can use the foolowing command to obtain unique result.
	python cnnLSV_rmrd.py <multi.vcf> <unique.vcf>
---

## Datasets
HG00512, HG00513, HG00514, HG00731, HG00732, HG00733, NA19238, NA19239 and NA19240 datasets can be downloaded from:

http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/

HG002 CLR and HG002 CCS datasets can be downloaded from:

https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/

## Contact
huidongma@st.gxu.edu.cn

