# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import runtools as rtools
import filtres as fil
import tempfile
from Bio import SeqIO
import argparse
import sys
import subprocess

class creole (rtools.ExternalTools):
    """
    Reference-free identification of open reading frames and encoded proteins from nanopore transcriptomic long reads. 
    """
    pass

    def __init__(self, input_read, output_file, threads, strand):
        """ 
        By default the strand is both(cDNA), other possible options are direct and reverse
        """
        super().__init__(input_read, output_file, threads, strand) #inheritance rtools.ExternalTools
        
    def __call__(self, input_read, output_file, threads=12, strand="both"):
        self.input_read = input_read
        self.output_file = output_file
        self.threads = threads
        self.strand = strand # strand - direct, reverse and both


    def join_headers_orfs(self, orf_faa):
        """
        Simplifies read ids from open reading frame
        """
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
        jheader = str(orf_faa)[:-4]+"join_h"
#        SeqIO.write(newseqs, jheader, 'fasta')
        return jheader
    
    def change_dic_write(func): #dicionary, output_file
        """
        Writes best read match to reference in PAF format
        """
        def __changes(dicionary, output_file):
            try:
                with open(output_file + "_bestRead_mapped.paf", "w") as write_dic:
                    for key, val in dicionary.items():
                        val = str(val).replace("'", '').replace(",", '').replace("[", '')\
                        .replace("]", '').replace("(", '').replace(")", '').replace(" ", '\t')
                        write_dic.writelines([str(key) + "\t", str(val) + "\n"])
            except Exception as e:
                print(e)
            func(dicionary, output_file)
        return __changes
    
    @change_dic_write
    def dic_write(self, dicionary, output_file):
        with open(output_file + "_bestRead_mapped.paf", "w") as write_dic:
                for key, val in dicionary.items():
                    val = str(val)
                    write_dic.writelines([str(key) + "\t", str(val) + "\n"])

    def precentage_identity_readP_anntP(self, bm_ids, out_prec_iden, orfs_reads, annt_prot):
        """
        Calculate the percentage of identity between read ORFs and annotated protain
        """
        with open(out_prec_iden + ".txt", "w") as pIdentity:
            fields = ("read","txt","len_target","len_query","match","mismatch","per_identity")
            pIdentity.write("{0}{6}{1}{6}{2}{6}{3}{6}{4}{6}{5}{7}".format(fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],"\t","\n"))
            indx_orfs_reads = dict()
            duplicate = {}
            with open(orfs_reads, "rU") as fasta, open ("temp_orf.faa", "w+") as fp:
                for read_record in SeqIO.parse(fasta, "fasta"):
                    if read_record.id in indx_orfs_reads.keys():
                        print("Duplicate id: %s" % read_record.id)
                        duplicate[read_record.id] =+1
                    else:
                        indx_orfs_reads[read_record.id] = read_record.seq
                        SeqIO.write(read_record, fp, "fasta")
                indx_orfs_reads = SeqIO.index("temp_orf.faa", "fasta") #Query

            indx_annt_protn = SeqIO.index(annt_prot, "fasta") #Target
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
                        record1 = indx_orfs_reads[read] #query
                        record2 = indx_annt_protn[txt_dic[txt]] #target
                        records = (record1, record2)
                        handle = StringIO()
                        SeqIO.write(records, handle, "fasta")
                        muscle_cline = MuscleCommandline(clwstrict=True)#clwstrict , msf=True
                        data = handle.getvalue()
                        stdout, stderr = muscle_cline(stdin=data)
                        align = AlignIO.read(StringIO(stdout), "clustal")
                        query = str(align[0].seq)
                        target = str(align[1].seq)
                        match = 0
                        mismatch = 0
                        for t, q in zip(target, query):
                            if t == q:  
                                match += 1
                            else:
                                mismatch += 1
                        pIdentity.writelines(str(read)+"\t"+str(txt)+"\t"+str(len(query))+"\t"+str(len(target))+"\t"+str(match)+"\t"+str(mismatch)+"\t"+str((match*100/len(target)))+"\n")
                    except Exception as e:
                        print(read, e)
                else:
                    if read in indx_orfs_reads:
                        no_txt = no_txt + 1
                        with open(out_prec_iden + "_transcpNoFound.txt", "a+") as tnf:
                            tnf.writelines(str(txt)+"\n")
                    elif txt in txt_dic:
                        no_read = no_read + 1
                        with open(out_prec_iden + "_readsNoFound.txt", "a+") as rnf:
                            rnf.writelines(str(read)+"\n")

        print("both: ", in_both, "read_no_found: ", no_read, "txt_no_found: ", no_txt)


def cmdline_args():
    try:
        parser = argparse.ArgumentParser(description="Parses command CREOLE")
        parser.add_argument("-gmst", "--gmst_path", required=False, help="points the direction of the gmst.pl file")
        parser.add_argument("-reads", "--input_read", required=True, help="input file (fasta or fastq format)")
        parser.add_argument("-o", "--output_file", required=True, help="output file")
        parser.add_argument("-ct", "--cell_type", type=str, choices=['prok', 'euka'], help="type of cell(prok/euka)")
        parser.add_argument("-st", "--strand", type=str, choices=['direct', 'reverse', 'both'], default="both", help="<string> sequence strand to predict genes in\
                (default: 'both'; supported: direct, reverse and both )#")
        #For Diamond
        parser.add_argument("-diamond", "--diamond_path", required=False, help="diamond path")
        parser.add_argument("-diamond_bd", "--diamond_bd", required=False, help="diamond BD")
        parser.add_argument("-t", "--threads", type=int, default=12, help="number of threads")
        #Criole
        parser.add_argument("-f", "--filter_type", type=int, default=2, choices=[1, 2, 3, 4], required=False,
                            help="choose a filter type: 1 - my score, 2 - mapping quality score, 3 - max mapping quality (60), 4 - DP alignment score ")
        #parser.add_argument("-i", "--input-paf", required=True, help="input file (paf format)")
        parser.add_argument("-minimap2", "--minimap2_path", required=False, help="minimap2_path")
        parser.add_argument("-txtRef", "--txt_reference", required=False, help="reference to map the read (fasta or fastq format)")
        parser.add_argument("-protRef", "--protein_reference", required=False, help="annoted protein")
        parser.add_argument("-v", "--verbose", dest='verbose', action='store_true', help="Verbose mode.")
        return (parser.parse_args())
    except argparse.ArgumentError as e:
        print("argparse error: ", e)
    
    

#    crle = creole(input_read=args.input_read, output_file=args.output_file, threads=args.treads, strand=orgs.strand)

if __name__ == "__main__":
    
    args = cmdline_args()
    
    if args.verbose:
        print("Verbose mode on")
    else:
        print("Verbose mode off")
    
    print(args.strand)
    crle = creole(input_read=args.input_read, output_file=args.output_file, threads=args.threads, strand=args.strand)
    crle.gmst(gmstPath=args.gmst_path)

#    crle.diamond(diamondPath=args.diamond_path, diamond_db=args.diamond_bd, faa_input="/home/jmaky/NetBeansProjects/CREOLE/creole/tests/joe.faa")
    
#    crle = creole(args.input_read, args.output_file, 12)
#    print(crle.gmst())
    
#    
#    
#    
#    run_gmst (args.gmst_path, args.input_read, args.output_file, args.cell_type, args.strand)
#    print("hecho gmst")
#    #DIAMOND
#
#    with open(args.output_file+".faaa","w") as out:
#        outp = args.output_file+".faa"
#        subprocess.call(['sed', '/^$/d', outp], stdout=out)
#
#    run_diamond(args.diamond_path, args.diamond_bd, args.output_file+".faaa", args.output_file+"_dt.m8")
#    print("Heacho diamond")
#    #CREOLE_4
#    #mapping read with minimap2
#    read_mapped = map_read_mm(args.input_read, args.map_reference, args.output_file, args.minimap2_path)
#    # r_m = read_mapped.stdout.decode("utf-8")
#    if args.filter_type == 1:
#        dic_write(select_best_read_myscore(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
#    elif args.filter_type == 2:
#        dic_write(select_best_read_MapQ(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
#        # paf_f = pathlib.Path(args.output_file + "_minimap2.paf")
#        # if (paf_f).exists():
#        #     pass
#        # else:
#        #     dic_write(select_best_read_MapQ(paf_f, args.output_file), args.output_file)
#    elif args.filter_type == 3:
#        dic_write(select_highest_MapQ(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
#    elif args.filter_type == 4:
#        dic_write(slct_best_read_sw_nm_score(args.output_file + "_unfiltered_minimap2.paf", args.output_file), args.output_file)
#
#    #select_best_read_myscore(args.output_file + "_full_mapped.paf")
#
#    with open(args.output_file + "_bestRead_mapped.paf", "r") as bmpaf, \
#        open(args.output_file + "_bestRead_mapped.id", "w") as ids_save:
#        for line in bmpaf:
#            line = line.split("\t")
#            ids_save.writelines(str(line[1-1]) + "\t" + str(line[6-1]) + "\n")
#
#
#    #perc identity
#    bm_ids = args.output_file + "_bestRead_mapped.id"
#    out_prec_iden = args.output_file + "_percentage_identity"
#    precentage_identity_readP_anntP(bm_ids, out_prec_iden, join_headers_orfs(args.output_file+".faaa"), args.protein_reference)
#
#    print ("done", args.output_file)
#

#if __name__ == "__main__":
#    print("ei")
