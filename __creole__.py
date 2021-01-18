#!/usr/bin/env python3
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from io import StringIO
import os
import pathlib
import subprocess
import sys
import tempfile

try:
    from io import StringIO # for Python 3
except:
    from StringIO import StringIO # for Python 2
#Run GeneMark
def run_gmst(gmst_path, inp_f, outp_f, cell_type="prok", strand="direct"): #euka
    if cell_type is "prok": # Check type of cell (prok or euka)
        try:
            subprocess.run(
                           [
                           'perl', gmst_path,
                           '--format', 'GFF',
                           '--strand', strand,
                           '--faa',
                           '--fnn',
                           '--prok', #this option is the same as:  --bins 1  --filter 0  --order 2  --order_non 2  --gcode 11 --width 6  --prestart 40 --fixmotif 0
                           '--output', outp_f, inp_f,
                           ],
                           stdout=subprocess.PIPE,
                           check=True,
                           universal_newlines=True,
                           )
        except subprocess.CalledProcessError as err:
            print('ERROR: ', err)
    else:
        try:
            subprocess.run(
                           [
                           'perl', gmst_path,
                           '--format', 'GFF',
                           '--strand', strand,
                           '--faa',
                           '--fnn',
                           '--output', outp_f, inp_f,
                           ],
                           stdout=subprocess.PIPE,
                           check=True,
                           universal_newlines=True,
                           )
        except subprocess.CalledProcessError as err:
            print('ERROR: ', err)
    return inp_f + ".faa"
#Run DIAMOND blastp
def run_diamond(diamond_path, diamond_db, inp_file_dm, out_file_dm):
    try:
        subprocess.run(
                       [
                       diamond_path,
                       'blastp',
                       '-d', diamond_db,
                       '-q', inp_file_dm,
                       '-o', out_file_dm,
                       '--strand', 'plus',
                       '--threads', '96',
                       '--outfmt', '6',
                       'qseqid',
                       'sseqid',
                       'pident',
                       'length',
                       'qlen',
                       'slen',
                       'mismatch',
                       'gapopen',
                       'qstart',
                       'qend',
                       'sstart',
                       'send',
                       'evalue',
                       'bitscore',
                       ],
                       stdout=subprocess.PIPE,
                       check=True,
                       universal_newlines=True,
                       )
    except subprocess.CalledProcessError as err:
        print('ERROR: ', err)

#join headers
def join_headers_orfs(orf_faa):
    get_headers = subprocess.run(["grep", ">", orf_faa], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    dic_ids = {}
    with tempfile.TemporaryFile() as fp:
        fp.write(get_headers.stdout)
        fp.seek(0)
        for i in fp.readlines():
            aux = i.decode("utf-8").split(">")[1].split()
            dic_ids[aux[0]] = str(aux[0])
    newseqs = []
    new_names = dic_ids
    for record in SeqIO.parse(orf_faa, 'fasta'):
        newid = new_names[record.id]
        newseqs.append(SeqRecord(record.seq, id=newid, description=''))
    SeqIO.write(newseqs, orf_faa + "_join_h", 'fasta')
    return orf_faa + "_join_h"

# map reads against the reference transcript using minimap2
def map_read_mm(read, ref, outp, mm):
    paf_file = pathlib.Path(outp + "_unfiltered_minimap2.paf")
    if (paf_file).exists():
        pass
    else:
        read_mapped = subprocess.run([mm, "-t7", "-cx", "map-ont", ref, read], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(paf_file, 'w') as sfile:
            sfile.write(read_mapped.stdout.decode("utf-8"))
        return read_mapped

# Save dictionary to text file
def dic_write(dicionary, output_file):
    try:
        with open(output_file + "_bestRead_mapped.paf", "w") as write_dic:
            for key, val in dicionary.items():
                val = str(val).replace("'", '').replace(",", '').replace("[", '')\
                .replace("]", '').replace("(", '').replace(")", '').replace(" ", '\t')
                write_dic.writelines([str(key) + "\t", str(val) + "\n"])
    except Exception as e:
        print(e)

def select_best_read_myscore(in_paf, output_ids):
    dic_bestMap = {}
    for line in open(in_paf):
        line = line.rstrip().split('\t')
        """ myScore = division of the number of matched bases between the query and target,
        by the length of the target included gaps * 100 """
        myScore = (float(line[10-1]) / float(line[11-1])) * 100
        if not line[0] in dic_bestMap:
            dic_bestMap[line[0]] = line[1:], myScore
        elif myScore > float(dic_bestMap.get(line[0])[-1]):
            dic_bestMap[line[0]] = line[1:], myScore
        with open(output_ids + "_bestRead_mapped.ids", "a+") as file_bm_ids:
            file_bm_ids.writelines(str(line[1-1]) + "\t" + str(line[6-1]))
    return dic_bestMap

def best_match_ids(line):
    dic_ids_bm = {}
    dic_ids_bm[line[1-1]] = str(line[6-1])
    return dic_ids_bm

def write_match_ids(dicionary, output_file):
    with open(output_file, "w") as file_bm_ids:
        for key, val in dicionary.items():
            write_dic.writelines([str(key) + "\t", str(val).split("") + "\n"])

def select_best_read_MapQ(in_paf, output_ids, threshold=0):
    dic_bestMapQ = {}
    for line in open(in_paf):
        line = line.rstrip().split('\t')
        mq = int(line[12-1])
        if threshold >= mq:
            sw = int(line[15-1].split(':')[2])
            nm = int(line[13-1].split(':')[2])
            if not line[0] in dic_bestMapQ:
                dic_bestMapQ[line[0]] = line[1:], nm, sw, mq
            elif mq > int(dic_bestMapQ.get(line[0])[-1]):
                dic_bestMapQ[line[0]] = line[1:], nm, sw, mq
            elif mq == int(dic_bestMapQ.get(line[0])[-1]) \
                and sw > int(dic_bestMapQ.get(line[0])[-2]):
                dic_bestMapQ[line[0]] = line[1:], nm, sw, mq
            elif mq == int(dic_bestMapQ.get(line[0])[-1]) \
                and sw == int(dic_bestMapQ.get(line[0])[-2]) \
                and nm < int(dic_bestMapQ.get(line[0])[-3]):
                dic_bestMapQ[line[0]] = line[1:], nm, sw, mq
    return dic_bestMapQ

def select_highest_MapQ(in_paf, output_ids):
    dic_highestMapQ = {}
    for line in open(in_paf):
        line = line.rstrip().split('\t')
        mq = int(line[12-1])
        sw = int(line[15-1].split(':')[2])
        nm = int(line[13-1].split(':')[2])
        if not line[0] in dic_highestMapQ and mq == 60:
            dic_highestMapQ[line[0]] = line[1:], nm, sw, mq
        elif mq == 60 and sw > int(dic_highestMapQ.get(line[0])[-2]):
            dic_highestMapQ[line[0]] = line[1:], nm, sw, mq
        elif mq == 60 and sw == int(dic_highestMapQ.get(line[0])[-2]) \
            and nm < int(dic_highestMapQ.get(line[0])[-3]):
            dic_highestMapQ[line[0]] = line[1:], nm, sw, mq
        with open(args.output_ids + "_bestRead_mapped.ids", "a+") as file_bm_ids:
            file_bm_ids.writelines(str(line[1-1]) + "\t" + str(line[6-1]))
    return dic_highestMapQ

def slct_best_read_sw_nm_score(in_paf, output_ids): #AS 15-1 | NM 13-1
    dic_AS_NM_score = {}
    for line in open(in_paf):
        line = line.rstrip().split('\t')
        sw = int(line[15-1].split(':')[2])
        nm = int(line[13-1].split(':')[2])
        if not line[0] in dic_AS_NM_score:
            dic_AS_NM_score[line[0]] = line[1:], sw, nm
        elif sw > int(dic_AS_NM_score.get(line[0])[-2]):
            dic_AS_NM_score[line[0]] = line[1:], sw, nm
        elif nm < int(dic_AS_NM_score.get(line[0])[-1]):
            dic_AS_NM_score[line[0]] = line[1:], sw, nm
        with open(args.output_ids + "_bestRead_mapped.ids", "a+") as file_bm_ids:
            file_bm_ids.writelines(str(line[1-1]) + "\t" + str(line[6-1]))
    return dic_AS_NM_score

def get_txt_from_anntProtein(annt_prot):
    dic_txts = {}
    for i in annt_prot:
        dic_txts[i.split("|")[1]] = i
    return dic_txts

def precentage_identity_readP_anntP(bm_ids, out_prec_iden, orfs_reads, annt_prot):
    with open(out_prec_iden + ".txt", "w") as pIdentity:
        pIdentity.write("read" + "\t" + "txt" + "\t" + "per_identity" + "\t" + "len_alig" + "\t" + "match" + "\t" + "mismatch" + "\n")
        indx_orfs_reads = SeqIO.index(orfs_reads, "fasta")
        indx_annt_protn = SeqIO.index(annt_prot, "fasta")
        txt_dic = get_txt_from_anntProtein(indx_annt_protn)
        in_both = 0
        no_read = 0
        no_txt = 0
        for map_ids in open(bm_ids):
            df = map_ids.split()
            txt = str(df[1]).split("|")[0]
            read = str(df[0])
            if read in indx_orfs_reads and txt in txt_dic:
                try:
                    in_both += 1
                    record1 = indx_orfs_reads[read]
                    record2 = indx_annt_protn[txt_dic[txt]]
                    records = (record1, record2)
                    handle = StringIO()
                    SeqIO.write(records, handle, "fasta")
                    muscle_cline = MuscleCommandline(clwstrict=True)#clwstrict , msf=True
                    data = handle.getvalue()
                    stdout, stderr = muscle_cline(stdin=data)
                    align = AlignIO.read(StringIO(stdout), "clustal")
                    target = str(align[0].seq)
                    query = str(align[1].seq)
                    match = 0
                    mismatch = 0
                    for t, q in zip(target, query):
                        if t == q:
                            match += 1
                        else:
                            mismatch += 1
                    pIdentity.writelines(str(read) + "\t" + str(txt) + "\t" + str((match * 100 / len(target))) + "\t" + str(len(target)) + "\t" + str(match) + "\t" + str(mismatch) + "\n")
                except Exception as e:
                    print(e)
            else:
                if read in indx_orfs_reads:
                    no_txt = no_txt + 1
                    with open(out_prec_iden + "_transcpNoFound.txt", "a+") as tnf:
                        tnf.writelines(str(txt) + "\n")
                elif txt in txt_dic:
                    no_read = no_read + 1
                    with open(out_prec_iden + "_readsNoFound.txt", "a+") as rnf:
                        rnf.writelines(str(read) + "\n")

    print("both: ", in_both, "read_no_found: ", no_read, "txt_no_found: ", no_txt)


def main():
    parser = argparse.ArgumentParser(description="Filter the best read alignment")
    parser.add_argument("-gmst", "--gmst_path", required=True, help="gmst.pl file")
    parser.add_argument("-reads", "--input_read", required=True, help="input file (fasta or fastq format)")
    parser.add_argument("-o", "--output_file", required=True, help="output file")
    parser.add_argument("-ct", "--cell_type", help="type of cell(prok/euka)")
    parser.add_argument("-st", "--strand", help="<string> sequence strand to predict genes in\
            (default: 'both'; supported: direct, reverse and both )")
    parser.add_argument("-diamond", "--diamond_path", required=True, help="diamond path")
    parser.add_argument("-diamond_bd", "--diamond_bd", required=True, help="diamond BD")
    parser.add_argument("-f", "--filter_type", type=int, choices=[1, 2, 3, 4], required=True,
                        help="choose a filter type: 1 - my score, 2 - mapping quality score, 3 - max mapping quality (60), 4 - DP alignment score ")
    parser.add_argument("-minimap2", "--minimap2_path", required=True, help="minimap2_path")
    parser.add_argument("-mr", "--map_reference", required=True, help="reference to map the read (fasta or fastq format)")
    parser.add_argument("-prot", "--protein_reference", required=True, help="annoted protein")
    args = parser.parse_args()

    #GMST
    run_gmst (args.gmst_path, args.input_read, args.output_file, args.cell_type, args.strand)
    print("hecho gmst")
    #DIAMOND
    with open(args.output_file + ".faaa", "w") as out:
        outp = args.output_file + ".faa"
        subprocess.call(['sed', '/^$/d', outp], stdout=out)
    run_diamond(args.diamond_path, args.diamond_bd, args.output_file + ".faaa", args.output_file + "_dt.m8")
    print("Heacho diamond")
    #CREOLE_4
    read_mapped = map_read_mm(args.input_read, args.map_reference, args.output_file, args.minimap2_path)
    if args.filter_type == 1:
        dic_write(select_best_read_myscore(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
    elif args.filter_type == 2:
        dic_write(select_best_read_MapQ(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
    elif args.filter_type == 3:
        dic_write(select_highest_MapQ(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
    elif args.filter_type == 4:
        dic_write(slct_best_read_sw_nm_score(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
    with open(args.output_file + "_bestRead_mapped.paf", "r") as bmpaf, \
        open(args.output_file + "_bestRead_mapped.id", "w") as ids_save:
        for line in bmpaf:
            line = line.split("\t")
            ids_save.writelines(str(line[1-1]) + "\t" + str(line[6-1]) + "\n")
    #perc identity
    bm_ids = args.output_file + "_bestRead_mapped.id"
    out_prec_iden = args.output_file + "_percentage_identity"
    precentage_identity_readP_anntP(bm_ids, out_prec_iden, join_headers_orfs(args.output_file + ".faaa"), args.protein_reference)

    print ("done", args.output_file)

if __name__ == "__main__":
    main()
