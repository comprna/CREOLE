#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

if __name__ == "__main__":
    
    diam = "/home/jmaky/Documents/projects/creoleT3/datas/precision/diamond/cdna/hop1/hop1_canu_dt.m8"
    out = "new_canu_cdna_hop1.m8"
    
    with open(diam) as f, open(out, "w") as outf:
        name = ("qseqid",
                "sseqid",
                "pident",
                "length",
                "qlen",
                "slen",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore")
        outf.writelines([i + "\t" for i in name])
        
        cnt = 1
        transcripts = list()
        new_id = dict()
        duplicated = dict()
        import time, sys
        print("Loading...")
        for r in f:
            time.sleep(0.1)
            line = list(r.split())
            _keys = line[0]
            bits = line[len(line)-1]
            if not _keys in new_id.keys():
                new_id[_keys] = list(line[1:])
                duplicated[_keys] = [1, line[1]]
            elif float(bits) > float(new_id.get(_keys)[-1]):
                new_id[_keys] = list(line[1:])
                duplicated[_keys].append(line[1])
                duplicated[_keys][0] += 1
            else:
                duplicated[_keys].append(line[1])
                duplicated[_keys][0] += 1
#                print("{}: {} {} more than {} {}".format(_keys, new_id.get(_keys)[0], new_id.get(_keys)[-1], bits, line[1]))
            sys.stdout.write(u"\u001b[1000D" + str(cnt + 1) + " reads")
            sys.stdout.flush()
            cnt += 1
#            if cnt == 10500: break
        for k, v in new_id.items():
            txt = "{}{}{}".format(str(k), "\t", str(v).replace("[", "").replace("]", "").replace("'", "").replace(",", "\t"))
            outf.writelines("{}{}".format("\n", txt))
        with open("diamond_duplicate_transc.txt", "w") as dupl_txt:
            dupl_txt.writelines("reads \t duplication_count \t txt \n")
            for x, y in duplicated.items():
                dp = "{}{}{}".format(str(x), "\t", str(y).replace("[", "").replace("]", "").replace("'", "").replace(",", "\t"))
                dupl_txt.writelines("{}{}".format("\n", dp))
