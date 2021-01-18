# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


import subprocess
import time
import sys
import os
import os.path
"""
This module shows documentation available for use with pydoc.
To generate HTML documentation for this module issue the
command:
       pydoc -w ExternalTools
"""

class Clrs:
    def __init__(self):
        self.HEADER = '\033[95m'
        self.OKBLUE = '\033[94m'
        self.OKCYAN = '\033[96m'
        self.OKGREEN = '\033[92m'
        self.WARNING = '\033[93m'
        self.FAIL = '\033[91m'
        self.ENDC = '\u001b[0m' #'\033[0m'
        self.BOLD = '\033[1m'
        self.UNDERLINE = '\033[4m'


class ExternalTools:
    """
    Run some external tools used by CREOLE
    """
    pass

    def __init__(self, input_read, output_file, threads, strand):
        self.input_read = input_read
        self.output_file = output_file
        self.threads = threads
        self.strand = strand # strand - direct, reverse and both
        

    def gmst(self, gmstPath, cellType="prok"):
        """
        Run GeneMarkS-T algorithm
        It is necessary to differentiate the type of cell:
            Eukaryotic (euka) or prokaryotic (prok) | by default: prok
            cellType=prok: this option is the same as:
            --bins 1  --filter 0  --order 2  --order_non 2  --gcode 11 --width 6  --prestart 40 --fixmotif 0
            in GeneMark-ES algorithm
        """
        
        start_time = time.time()
        cwd = os.getcwd()
        print("-------------------------------------------------\n\n")
        print("{}***CREOLE starts predicting ORFs using GeneMarkS-T Algorithm...{}".format(Clrs().OKGREEN, Clrs().ENDC))
        
        print("Start time\t: {}".format(time.ctime()).expandtabs(15))
        print("Working dir\t: {}".format(cwd).expandtabs(15))
        print("command line\t: --format {} --strand {} --faa --fnn --prok --output --verbose".format("GFF",self.strand,self.output_file))
        print("Gmst\t: {}".format(gmstPath).expandtabs(15))
        print("Cell\t: {}".format(cellType).expandtabs(15))
        print("Input\t: {}".format(self.input_read).expandtabs(15))
        print("Output\t: {}".format(self.output_file).expandtabs(15))
        print("Strand\t: {}".format(self.strand).expandtabs(15))
        print("Threads\t: {}".format(self.threads).expandtabs(15))

        if os.path.isfile(self.input_read):
            pass
        else:
            print("{}ERROR: {} not exist{}Please check and try again".format(Clrs().FAIL, "\n",self.input_read) +Clrs().ENDC)
            
        if cellType is "prok":
            default_par = [
                            'perl', gmstPath,
                            '--format', 'GFF',
                            '--strand', self.strand, #direct, reverse and both
                            '--faa',
                            '--prok', #this option is the same as:  --bins 1  --filter 0  --order 2  --order_non 2  --gcode 11 --width 6  --prestart 40 --fixmotif 0
                            '--output', self.output_file, self.input_read,
                            '--verbose',
                            ]
            default_par.append('--fnn')
            try:
                subjob=subprocess.run(
                                default_par,
                                stdout=subprocess.PIPE,
                                check=True,
                                universal_newlines=True,
                               )
                print(subjob.stdout)
            except subprocess.CalledProcessError as err:
                print(Clrs().FAIL +'ERROR: Please check the parameters and try again \n'+Clrs().ENDC, err)
            else:
                print(Clrs().OKGREEN + "\t********|ORFs successfully predicted|********".expandtabs(15) +Clrs().ENDC)
                print(Clrs().OKCYAN + "\truntime: %s seconds" % (time.time() - start_time),"|".expandtabs(2) + Clrs().ENDC+"\n")
                print("Orfs\t: {}{}".format(self.output_file, ".faa").expandtabs(15))
        else: #euka
            try:
                subprocess.run(
                               [
                               'perl', gmstPath,
                               '--format', 'GFF',
                               '--strand', strand,
                               '--faa', '--fnn',
                               '--output', self.output_file, self.input_read,
                               ],
                               stdout=subprocess.PIPE,
                               check=True,
                               universal_newlines=True,
                               )
            except subprocess.CalledProcessError as err:
                print('ERROR: ', err)
            else:
                print("ORFs successfully predicted")
        return self.input_read + ".faa"


    def diamond(self, diamondPath, diamond_db, faa_input): #diamond_path, diamond_db, inp_file_dm, out_file_dm
        """
        Run Diomond's blastp module
        """
        tsv_output = self.output_file + "_infiltred_diamond_out.tsv"
        print("{}\n{}\n{}\n{}".format(diamondPath, diamond_db, faa_input, tsv_output))
        
        
        if self.strand is "direct":
            self.strand = "plus"
        elif self.strand is "reverse":
            self.strand = "minus"
        try:
            subprocess.run(
                           [
                           diamondPath,
                           'blastp',
                           '-d', diamond_db,
                           '-q', faa_input,
                           '-o', tsv_output,
                           '--strand', self.strand, #both, plus, minus}
                           '--threads', self.threads,
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


    def map_read_mm2(self, minimap2Path, ref):
        """
        Run the Minimap map-ont module
        """
        paf_file = pathlib.Path(self.output_file + "_unfiltred_mm2.paf")
        if (paf_file).exists():
            print("The {} file already exists".format(self.output_file + "_unfiltred_mm2.paf"))
            pass
        else:
            threads = "-t" + self.threads
            read_mapped = subprocess.run([minimap2Path, threads, "-cx", "map-ont", ref, self.input_read], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with open(paf_file, 'w') as sfile:
                sfile.write(read_mapped.stdout.decode("utf-8"))
            return read_mapped
