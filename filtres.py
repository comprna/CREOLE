# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

class Filters:
    """
    Filter best read matched
    """
    pass

    def __init__(self, in_paf, output_ids):
        self.in_paf = in_paf
        self.output_ids = output_ids
        

    def best_read_myscore(self):
        """
       Division of the number of matched bases between the query and target,
            by the length of the target included gaps * 100
            
        """
        line = None
        dic_bestMap = {}
        for line in open(self.in_paf):
            line = line.rstrip().split('\t')
            myScore = (float(line[10-1]) / float(line[11-1])) * 100
            if not line[0] in dic_bestMap:
                dic_bestMap[line[0]] = line[1:], myScore
            elif myScore > float(dic_bestMap.get(line[0])[-1]):
                dic_bestMap[line[0]] = line[1:], myScore
            with open(self.output_ids + "_bestRead_mapped.ids", "a+") as file_bm_ids:
                file_bm_ids.writelines(str(line[1-1]) + "\t" + str(line[6-1]))
        return dic_bestMap

    def best_match_ids(self, line):
        """
        Set read Ids with Transcript Ids
        """
        dic_ids_bm = {}
        dic_ids_bm[line[1-1]] = str(line[6-1])
        return dic_ids_bm

    def write_match_ids(self, dicionary, output_file):
        """
        write best match Ids ind TAB format
        """
        with open(output_file, "w") as file_bm_ids:
            for key, val in dicionary.items():
                write_dic.writelines([str(key) + "\t", str(val).split("") + "\n"])

    def best_read_MapQ(self, threshold=0):
        """
        Quantify the probability that a read is misplaced. 
        In PAF format, MapQ generally ranges from 0 to 60, with 60 being the optimal value.
        """
        dic_bestMapQ = {}
        for line in open(self.in_paf):
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
            with open(self.output_ids + "_bestRead_mapped.ids", "a+") as file_bm_ids:
                file_bm_ids.writelines(str(line[1-1]) + "\t" + str(line[6-1]))
        return dic_bestMapQ

    def highest_MapQ(self):
        """
        Select the maximum MapQ in the PAF file
        """
        dic_highestMapQ = {}
        for line in open(self.in_paf):
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
            with open(self.output_ids + "bestRead_mapped.ids", "a+") as file_bm_ids:
                file_bm_ids.writelines(str(line[1-1]) + "\t" + str(line[6-1]))
        return dic_highestMapQ

    def best_read_sw_nm_score(self): #AS 15-1 | NM 13-1
        """
        Finds the optimal alignment by using specific scores for matches and mismatches through a scoring matrix
        """
        dic_AS_NM_score = {}
        for line in open(self.in_paf):
            line = line.rstrip().split('\t')
            sw = int(line[15-1].split(':')[2])
            nm = int(line[13-1].split(':')[2])
            if not line[0] in dic_AS_NM_score:
                dic_AS_NM_score[line[0]] = line[1:], sw, nm
            elif sw > int(dic_AS_NM_score.get(line[0])[-2]):
                dic_AS_NM_score[line[0]] = line[1:], sw, nm
            elif nm < int(dic_AS_NM_score.get(line[0])[-1]):
                dic_AS_NM_score[line[0]] = line[1:], sw, nm
            with open(self.output_ids + "_bestRead_mapped.ids", "a+") as file_bm_ids:
                file_bm_ids.writelines(str(line[1-1]) + "\t" + str(line[6-1]))
        return dic_AS_NM_score

    def get_txt_from_anntProtein(self, annt_prot):
        """
        Get transcripts from annotated protein
        """
        dic_txts = {}
        for i in annt_prot:
            dic_txts[i.split("|")[1]] = i
        return dic_txts

