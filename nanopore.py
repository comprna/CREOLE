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
    
    if nanopore.split(".")[1] == "paf":
        
        for line in j:
#        sc = float(line.rstrip().split(" ")[2])
            line = line.rstrip().split("\t")[5]

            if not line in dic_ids:
                dic_ids[line] = 1
            else:
                dic_ids[line] += 1
    else:
        
        j.readline()
        
        for line in j:
#        sc = float(line.rstrip().split(" ")[2])
            line = line.rstrip().split(" ")[1]

            if not line in dic_ids:
                dic_ids[line] = 1
            else:
                dic_ids[line] += 1

#    j = open(nanopore)
#    j.readline()
#    dic_ids = {}
#    dic_ambigous = {}
#
#    for line in j:
##        sc = float(line.rstrip().split(" ")[2])
#        line = line.rstrip().split(" ")[1]
#        
#        if not line in dic_ids:
#            dic_ids[line] = 1
#        else:
#            dic_ids[line] += 1

    dic_ids = ambigous_nanopore(idAmabigous, dic_ids)

    # dic_write(dic_ids, "K562_rattle_nanoTPM_txt")
    no_nor = dic_ids
    N = sum(dic_ids.values())
    nor = 10**6
    dic_ids = {k:numpy.log10(nor * (n / N)) for (k,n) in dic_ids.items()}
    # dic_write(dic_ids, "count_transcripts_nor")
    return (no_nor, dic_ids)

idAmabigous = "../datas/id_annotation_ambigous"
def txtPerGene_nanopore(nanopore):
    dic = normalize_TPM_nanopore(nanopore, idAmabigous)[0]
    dic_gene = {}

    for key, val in dic.items():
        gen = key.split("|")[1]
        if not gen in dic_gene:
            dic_gene[gen] = val
        else:
            dic_gene[gen] += val

    N = sum(dic_gene.values())
    nor = 10**6
    dic_gene = {k:numpy.log10(nor * (n / N)) for (k,n) in dic_gene.items()}
       # dic_write(dic_ids, "count_genes")
    return dic_gene

def count_txtPerGene_nanopore(nanopore):
    dic = normalize_TPM_nanopore(nanopore, idAmabigous)[0]
    dic_gene = {}

    for key, val in dic.items():
        txt, gen = key.split("|")[:2]
        if not gen in dic_gene:
            dic_gene[gen] = 1
        else:
            dic_gene[gen] += 1

#    N = sum(dic_gene.values())
#    nor = 10**6
#    dic_gene = {k:numpy.log10(nor * (n / N)) for (k,n) in dic_gene.items()}
       # dic_write(dic_ids, "count_genes")
    return dic_gene

def selectGeneWith2txt_nanopore(nanopore):

    dic_gene_2 = {}
    dic_gene = count_txtPerGene_nanopore(nanopore)
    dic = normalize_TPM_nanopore(nanopore, idAmabigous)[0]
    print(dic)

    for key, val in dic_gene.items():
        gen = key.split("|")[1]
        print(gen)
        if not gen in dic_gene_2 and dic_gene[gen]== 2:
            dic_gene_2[gen] = val
        elif gen in dic_gene_2:
            dic_gene_2[gen] += val
    # dic_write(dic_gene, "count_genes_two_txt")
    print(dic_gene_2)
    N = sum(dic_gene_2.values())
    nor = 10**6
    dic_gene_2 = {k:numpy.log10(nor * (n / N)) for (k,n) in dic_gene_2.items()}
    
    return dic_gene_2
