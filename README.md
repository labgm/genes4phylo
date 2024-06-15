# genes4phylo
Script to generate files for genes of interest in separate directories for later phylogenetic analysis performance

## Usage

Before running the script it is necessary to download the genbank files from your organisms of interest

python gene4phylo.py --input files_list.txt --genes genes_of_interest.tsv --output result

  -h, --help       show this help message and exit
  
  --input INPUT    list of input genbank files.
  
  --genes GENES    List of genes of interest.
  
  --output OUTPUT  Output directory.
