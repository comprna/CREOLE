import nanopore_orfs as nano
import illumina as ill

import concat
import csv
import sys

# transcriptome ids ambigous
idAmabigous = "/home/joel/NetBeansProjects/LongReads/longreads/datas/id_annotation_ambigous"
# Paf nanopore
nn = sys.argv[1]
# TPMs illumina
ii = sys.argv[2]

# Txt level
txt = concat.concat_illimina_nanopore(ill.normalize_TPM_illumina(ii),\
    nano.normalize_TPM_nanopore(nn, idAmabigous)[1])
concat.ploting(txt)

# Gene level
gen = concat.concat_illimina_nanopore(ill.txtPerGene_illumina(ii),\
    nano.txtPerGene_nanopore(nn))
concat.ploting(gen)
