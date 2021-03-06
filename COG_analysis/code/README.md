
## Details

eggNOG based data and analysis for the RosR project. 

Created: 2018-08-16
Author: Rylee K. Hackley	rkh11@duke.edu

## Description
determine orthology between called RosR binding regions (identified via ChIP-Seq) for Haloferax med and volcanii. 

Annotation data comes from emapper run numbers MM_X2qKDu (H medd.) and (H vol.). on peak files provided by CD (file location). Annotations were computed using eggnog-mapper [1] based on eggNOG 4.5 orthology data [2].

eggnog details:
Seed orthologs are the most similar ortholog, and then eggNOGs are different ortholog groups based on taxa.
@halNOG = Halobacteria
@eurNOG = Euryarchaeota
@arNOG = Archaea
@NOG = all Domains

H medd. run info: http://eggnogdb.embl.de/#/app/emapper?jobname=MM_X2qKDu
H vol. run info: http://eggnogdb.embl.de/#/app/emapper?jobname=MM_4N0HNn

## References
[1]  Huerta-Cepas J, et. al. Fast genome-wide functional annotation through orthology assignment by eggNOG-mapper.Mol Biol Evol (2017).doi: 10.1093/molbev/msx148
[2] Huerta-Cepas J, et. al. eggNOG 4.5: a hierarchical orthology framework with improved functional annotations for eukaryotic, prokaryotic and viral sequences. Nucl. Acids Res. (04 January 2016) 44 (D1): D286-D293. doi: 10.1093/nar/gkv1248

## Edit log
2018-08-16
	add downloadable .tsv output to eggNOG directory
	add headers to raw data and save as 'preprocessed' .csv files
2018-08-17
	created data/ code/ output/ directories and sub directories 
	copied raw directories from online results page into data/raw directory 
	create compressed backup of raw data files
2018-08-21
	downloaded and added Halobacterium salinarum results to raw and preprocessed folders. created another instance of backup, to be handled via git at a later date.
	
