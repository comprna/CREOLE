# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

if __name__ == "__main__":
    print ("Hello World")


mnm = "/genomics/users/joel/CEPH1463/nanopore/cDna1Dpass/Minimap2/run_1/--secondary=no_Bham_Run1.paf"



import numpy

def ambigous_nanopore(idAmabigous, dic_count_txt):
    dic_ambigous = {}

    with open(idAmabigous) as id_ambig:
        for i in id_ambig:
            i = i.rstrip().split(',')
            dic_ambigous[i[0]] = i[1:]

        
    for key, value in dic_ambigous.items():
        
        if key in dic_count_txt:
            dic_count_txt[key] /=(len(value)+1)
            for k in value:
                dic_count_txt[k] = dic_count_txt[key]
    return dic_count_txt


def normalize_TPM_nanopore(nanopore, idAmabigous):

    j = open(nanopore)
    dic_ids = {}
    dic_ambigous = {}

    for line in j:
#        sc = float(line.rstrip().split(" ")[2])
        line = line.rstrip().split("\t")[5]
        
        if not line in dic_ids:
            dic_ids[line] = 1
        else:
            dic_ids[line] += 1

    dic_ids = ambigous_nanopore(idAmabigous, dic_ids)

    # dic_write(dic_ids, "K562_rattle_nanoTPM_txt")
    no_nor = dic_ids
    N = sum(dic_ids.values())
    nor = 10**6
    dic_ids = {k:numpy.log10(nor * (n / N)) for (k,n) in dic_ids.items()}
    # dic_write(dic_ids, "count_transcripts_nor")
    return (no_nor, dic_ids)

idAmabigous = "../datas/id_annotation_ambigous"