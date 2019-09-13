import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy

def dic_write(dicionary, output_file):

    with open(output_file, "w") as write_dic:
        w = csv.writer(write_dic)
        for key, val in dicionary.items():
            w.writerow([key, val])


def concat_illimina_nanopore(dic_illumina, dic_nanopore):
    dic_mix ={}

    for key, val in dic_nanopore.items():
        if key in dic_illumina:
            dic_illumina[key].append(val)

    for key in dic_illumina:
        value = dic_illumina[key]
        if len(dic_illumina[key]) == 3:
            dic_mix[key] = value

    df = pd.DataFrame.from_dict(data=dic_mix, orient="index")
    df.columns = ["ILLUMINA_1", "ILLUMINA_2", "NANOPORE"]

    return df

def concat_illimina_nanopore_rep1(dic_illumina, dic_nanopore):
    dic_mix ={}

    for key, val in dic_nanopore.items():
        if key in dic_illumina:# and val >= 1:
            dic_illumina[key].append(val)

    for key in dic_illumina:
        value = dic_illumina[key]
        if len(dic_illumina[key]) == 2 :
            dic_mix[key] = value
            
    df = pd.DataFrame.from_dict(data=dic_mix, orient="index")
    df.columns = ["ILLUMINA", "NANOPORE"]
    
    return df


def legend (x, y, **kws):
    r = x.corr(y)
    ax = plt.gca()
    ax.annotate("r = {:.2f}".format(r),
                xy=(.1, .9), xycoords = ax.transAxes )


def ploting(dataframe):
    sns.set_context("talk")
    g = sns.PairGrid(dataframe, dropna=True)
    g.map_upper(sns.regplot, color=".3", scatter_kws={'s':6})
    g.map_lower(legend)
    g.map_lower(sns.regplot, color=".3", scatter_kws={'s':6})
    g.map_diag(sns.kdeplot)
#    plt.suptitle("Cut-off (nanopore >= 1)", fontsize=12)
    plt.show(g)

    plt.savefig(inputPath+"transcripts.png")
    print("Done transcripts")
