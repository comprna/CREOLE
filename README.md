# CREOLE
Reference-free identification of open reading frames and encoded proteins from nanopore transcriptomic long reads. 
# Requirements
- Minimap2
- GeneMark
- DIAMOND
# Installation
- Clone the repository

`git clone --recurse-submodules https://github.com/comprna/CREOLE`

# Usage example 
**Warning** All the commands and parameters are still highly experimental and subject to changes in future versions

- Sample commands
```
python creole.py -gmst gmst.pl -reads read.fa -o file.out -ct prok -st direct -diamond diamond -diamond_bd nr.dmnd -f 2 -minimap2 minimap2 -mr gencode.v31.transcripts.fa -prot gencode_v31/gencode.v31.pc_translations.fa

```

- multiple files in cluster (Slurm)

```
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH -J creoleT4_RNA_mecat
#SBATCH --mem=91G
#SBATCH -c 24
#SBATCH -e error/creoleT4_RNA_mecat.err
#SBATCH -o error/creoleT4_RNA_mecat.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joelamaque@gmail.com

module load Python

export PATH=/users/genomics/joel/creole:$PATH
export PATH="$PATH:/genomics/users/joel/soft/ST_alg/gmst_linux_64"
export PATH="$PATH:/projects_rg/joel/tools/diamond"

dir="creole/RNA/mecat"
mkdir -p $dir
cd $dir

array=(`ls /projects_rg/joel/charlote_et_all/nanopore/corr_reads/RNA/mecat/*_cut_150_cor.fa_mix | cut -d"/" -f9 | cut -d"." -f1`)
array1=(`ls /projects_rg/joel/charlote_et_all/nanopore/corr_reads/RNA/mecat/*_cut_150_cor.fa_mix`)

for index in ${!array[@]}; do
        echo ${array[$index]} ${array1[$index]}
        python /users/genomics/joel/creoleT3/creole.py\
        -gmst /genomics/users/joel/soft/ST_alg/gmst_linux_64/gmst.pl\ 
        -reads ${array1[$index]}\ 
        -o ${array[$index]}\  
        -ct prok\ 
        -st direct\ 
        -diamond /projects_rg/joel/tools/diamond\ 
        -diamond_bd /users/genomics/joel/charlote_et_all/orfs/prot_db/nr.dmnd\ 
        -f 2\ 
        -minimap2 /projects_rg/joel/tools/minimap2-2.17_x64-linux/minimap2\ 
        -mr /projects_rg/joel/gencode_v31/gencode.v31.transcripts_noPseudogene.fa\ 
        -prot /projects_rg/joel/gencode_v31/gencode.v31.pc_translations.fa
done
```

- Commands and options

```
usage: creole.py [-h] [-gmst GMST_PATH] -reads INPUT_READ -o OUTPUT_FILE
                 [-ct {prok,euka}] [-st {direct,reverse,both}]
                 [-diamond DIAMOND_PATH] [-diamond_bd DIAMOND_BD] [-t THREADS]
                 [-f {1,2,3,4}] [-minimap2 MINIMAP2_PATH]
                 [-txtRef TXT_REFERENCE] [-protRef PROTEIN_REFERENCE] [-v]

Parses command CREOLE

optional arguments:
  -h, --help            show this help message and exit
  -gmst GMST_PATH, --gmst_path GMST_PATH
                        points the direction of the gmst.pl file
  -reads INPUT_READ, --input_read INPUT_READ
                        input file (fasta or fastq format)
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        output file
  -ct {prok,euka}, --cell_type {prok,euka}
                        type of cell(prok/euka)
  -st {direct,reverse,both}, --strand {direct,reverse,both}
                        <string> sequence strand to predict genes in (default:
                        'both'; supported: direct, reverse and both )#
  -diamond DIAMOND_PATH, --diamond_path DIAMOND_PATH
                        diamond path
  -diamond_bd DIAMOND_BD, --diamond_bd DIAMOND_BD
                        diamond BD
  -t THREADS, --threads THREADS
                        number of threads
  -f {1,2,3,4}, --filter_type {1,2,3,4}
                        choose a filter type: 1 - my score, 2 - mapping
                        quality score, 3 - max mapping quality (60), 4 - DP
                        alignment score
  -minimap2 MINIMAP2_PATH, --minimap2_path MINIMAP2_PATH
                        minimap2_path
  -txtRef TXT_REFERENCE, --txt_reference TXT_REFERENCE
                        reference to map the read (fasta or fastq format)
  -protRef PROTEIN_REFERENCE, --protein_reference PROTEIN_REFERENCE
                        annoted protein
  -v, --verbose         Verbose mode.
```
