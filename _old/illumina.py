#!/usr/bin/env python3
#encoding: UTF-8

import numpy

def normalize_TPM_illumina(illumina):

    with open(illumina) as illum:
        illum.readline()
        dic_illumina = {}

        for line in illum:
            line = line.rstrip().split("\t")

#            rep1, rep2 = numpy.log10(float(line[1])), numpy.log10(float(line[2])) Normal
#            rep1, rep2 = numpy.log10(float(line[1])), numpy.log10(float(line[2])) # Contol_fetal
            #modificado para 1rep
            rep1, rep2 = numpy.log10(float(line[1])), numpy.log10(float(line[1])) # Modificado

            if rep1 >= 0 and rep2 >= 0:
                dic_illumina[line[0]] = [rep1, rep2]
            elif rep1 >= 0 and rep2 < 0:
                dic_illumina[line[0]] = [rep1, numpy.nan]
            elif rep2 >= 0 and rep1 < 0:
                dic_illumina[line[0]] = [numpy.nan, rep2]

    return dic_illumina


def txtPerGene_illumina(illumina):

    with open(illumina) as illum:
        illum.readline()
        dic_aux = {}
        dic_illumina_gene = {}

        for line in illum:

            line = line.rstrip().split("\t")
#            v_rep1, v_rep2 = float(line[1]), float(line[2])
            v_rep1, v_rep2 = float(line[1]), float(line[2])
            gen = line[0].split("|")[1]

            if not gen in dic_aux:
                dic_aux[gen] = [v_rep1, v_rep2]
            else:
                dic_aux[gen]= [dic_aux[gen][i] + (v_rep1, v_rep2)[i] for i in range(len(dic_aux[gen]))]

        for key, values in dic_aux.items():

            rep1, rep2 = numpy.log10(values[0]), numpy.log10(values[1])

            if rep1 >= 0 and rep2 >= 0:
                dic_illumina_gene[key] = [rep1, rep2]
            elif rep1 >= 0 and rep2 < 0:
                dic_illumina_gene[key] = [rep1, numpy.nan]
            elif rep2 >= 0 and rep1 < 0:
                dic_illumina_gene[key] = [numpy.nan, rep2]

    return dic_illumina_gene

def count_txtPerGene_illumina(illumina):

    with open(illumina) as illum:
        illum.readline()
        dic_aux = {}

        for line in illum:

            line = line.rstrip().split("\t")
            key = line[0].split("|")[1]
#            value = [float(line[1]), float(line[2])]
            value = [float(line[1]), float(line[2])]

            if key in dic_aux:
                if value[0] > 0 < value[1]:
                    dic_aux[key][0] +=1
                    dic_aux[key][1] +=1
                elif value[0] > 0:
                    dic_aux[key][0] +=1
                elif value[1] > 0:
                    dic_aux[key][1] +=1
            else:
                if value[0] > 0 < value[1]:
                    dic_aux[key] =[1, 1]
                elif value[0] > 0:
                    dic_aux[key] =[1, 0]
                elif value[1] > 0:
                    dic_aux[key] = [0, 1]
    # dic_write(dic_aux, "dic_aux")
    return dic_aux

def selectGeneWith2txt_illumina(illumina):
    dic_count = count_txtPerGene_illumina(illumina)
    with open(illumina) as illum:
        illum.readline()
        dic_aux1 = {}

        for line in illum:

            line = line.rstrip().split("\t")
            key = line[0].split("|")[1]

            # value = [float(line[1]), float(line[2])]
#            value = [numpy.log10(float(line[1])), numpy.log10(float(line[2]))]
            value = [numpy.log10(float(line[1])), numpy.log10(float(line[2]))]
            # print(dic_aux[key])

            if key in dic_count:
                if key in dic_aux1:
                    if dic_count[key][0] == 2 == dic_count[key][1]:
                        dic_aux1[key] = [dic_aux1[key][i] + value[i] for i in range(len(dic_aux1[key]))]
                    elif dic_count[key][0] == 2:
                        dic_aux1[key][0] += value[0]
                    elif dic_count[key][1] == 2:
                        dic_aux1[key][1] += value[1]
                else:
                    # print(dic_aux[key])
                    if dic_count[key][0] == 2 == dic_count[key][1]:
                        dic_aux1[key] = value
                    elif dic_count[key][0] == 2:
                        dic_aux1[key] = [value[0], 0]
                    elif dic_count[key][1] == 2:
                        dic_aux1[key] = [0,value[1]]

    return dic_aux1
