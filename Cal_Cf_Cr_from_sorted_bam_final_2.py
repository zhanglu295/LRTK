#import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import sys
def flatten(listoflist):
    flat_list = []
    for alist in listoflist:
        for val in alist:
            flat_list.append(val)
    return flat_list


def process_sorted_bam(sam_file,output_file):
    f = open(sam_file,"r")
    fw= open(output_file,"w")
    mole_dict = defaultdict(list)
    mole_dict_2 = defaultdict(lambda: defaultdict(int))
    mole_dict_3 = defaultdict(lambda: defaultdict(list))
    mole_dict_save = defaultdict()
    curr = 0
    count_mole = 1
    for line in f:
        #print(curr)
        curr += 1
        data = line.rsplit()
        barcode_field = [s for s in data if "BX:Z:" in s]
        if barcode_field != []:
            barcode =  barcode_field[0].split(":")[2].split("-")[0]
            start_pos = int(data[3])
            chr_num = data[2]
            qname = data[0]
            if len(mole_dict[barcode]) == 0:
                mole_dict[barcode].append(start_pos) 
                mole_dict_2[barcode][qname] += 1 
                mole_dict_3[barcode][qname].append(start_pos) 
            elif len(mole_dict[barcode]) > 0:
                dist = start_pos - mole_dict[barcode][-1]
                if dist < 50000:   # 100kb
                    mole_dict[barcode].append(start_pos) 
                    mole_dict_2[barcode][qname] += 1
                    mole_dict_3[barcode][qname].append(start_pos)
                else:
                    for key,value in mole_dict_2[barcode].items():
                        if value == 1:
                            mole_dict_3[barcode].pop(key)
                    all_pos = mole_dict_3[barcode].values()
                    if len(all_pos) >= 2:  # 2 pairs
                       # print(all_pos)
                        mole_dict_save[(barcode,count_mole)] = all_pos
                        poss = flatten(all_pos)
                        start_ = min(poss)
                        end_ = max(poss)
                        mole_len = end_ - start_  + 151
                        fw.writelines(chr_num + "\t" + str(start_) + "\t" + str(end_) + "\t" + str(mole_len) + "\t" + str(len(poss)) + "\t" + str(barcode) + "\t" + str(count_mole) + "\n")
                        mole_dict[barcode] = [] 
                        mole_dict_2[barcode] = defaultdict(int)
                        mole_dict_3[barcode] = defaultdict(list) 
                        mole_dict[barcode].append(start_pos)
                        mole_dict_2[barcode][qname] += 1
                        mole_dict_3[barcode][qname].append(start_pos)
                        count_mole += 1
                    else:
                        mole_dict[barcode] = []
                        mole_dict_2[barcode] = defaultdict(int)
                        mole_dict_3[barcode] = defaultdict(list)
                        mole_dict[barcode].append(start_pos)
                        mole_dict_2[barcode][qname] += 1
                        mole_dict_3[barcode][qname].append(start_pos)


    # for the droplet only with one molecule
    for barcode,value in mole_dict.items():
        if len(value) > 0:
            for qname,num_reads in mole_dict_2[barcode].items():
                if num_reads == 1:
                    print(qname,num_reads)
                    mole_dict_3[barcode].pop(qname)
            all_pos = mole_dict_3[barcode].values()
            if len(all_pos) >= 2: # 2 pairs
                mole_dict_save[(barcode,count_mole)] = all_pos
                poss = flatten(all_pos)
                start_ = min(poss)
                end_ = max(poss)
                mole_len = end_ - start_  + 151
                fw.writelines(chr_num + "\t" + str(start_) + "\t" + str(end_) + "\t" + str(mole_len) + "\t" + str(len(poss)) + "\t" + str(barcode) + "\t" + str(count_mole) + "\n")

                count_mole += 1

    f.close()
    fw.close()
    return mole_dict



def read_pickle(mole_dict_save,output_file,chr_num):
    fw = open(output_file,"w")
    for key, raw_value in mole_dict_save.items():
        if type(raw_value[0]) == type([]):
            value = flatten(raw_value)
        else:
            value = raw_value
        start_pos = min(value)
        end_pos = max(value)
        mole_num = key[1]
        print(mole_num)
        mole_len = end_pos - start_pos  + 151
        fw.writelines(chr_num + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(mole_len) + "\t" + str(len(value)) + "\t" + str(key[0]) + "\t" + str(mole_num) + "\n")

    fw.close()



if __name__ == "__main__":
    chr_sam_file = sys.argv[1]
    chr_h5 = sys.argv[2]
    process_sorted_bam(chr_sam_file,chr_h5)
