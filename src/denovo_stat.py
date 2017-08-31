#!/usr/bin/env python
# coding: utf-8

# ### Import essential packages

# In[1]:
'''
from collections import defaultdict
import matplotlib
import re
import numpy
import operator
from scipy import stats as sta
import seaborn as sns
import math
from collections import defaultdict
from scipy.interpolate import spline
import operator
from itertools import groupby
import sys
import gzip
import getopt
'''
# ## Input Contig alignment

# In[2]:

class aligncontig(object):#class 'aligncontig' to store alignment information for contigs,
          def __init__(self):
            self.start=0 # start of contig alignment ,
            self.end=0   # end of contig alignment ,
            #self.length=0 # total length including gap,
            self.alignlength=0 #total length without misassemblies,
            self.forward=1
            self.misstart=[] #list for the starts of aligned contigs, which may include the breakpoints of misassemblies,
            self.misend=[]   #list for the ends of aligned contigs, which may include the breakpoints of misassemblies,
def alignment_contig_coordinate(tsv_file,chr_id,chr_len):
        f_tsv = open(tsv_file,'r')
        count = 0
        contig_dict = defaultdict(aligncontig)
        start_align=[]
        end_align=[]
        start_contig=[]
        end_contig=[]
        offset_array=[]
        offset_array.append(0)
        offset=0
        for length_chr in chr_len:
            offset=offset+length_chr
            offset_array.append(offset)
        alignlength=0
        processid=0
        for line in f_tsv:
            data=line.rsplit()
            if count==0:
                count+=1
                continue
            elif count%2 == 1:
                data=line.rsplit()
                if len(data)==4:
                    continue
                index_offset=chr_id.index(data[4])
                offset_add=offset_array[index_offset]
                if processid==0:
                    processid=int(data[5])
                    start_align.append(int(data[0])+offset_add)
                    end_align.append(int(data[1])+offset_add)
                    alignlength=alignlength+int(data[1])-int(data[0])+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
                elif processid==int(data[5]):
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
                else:
                    Contig=aligncontig()
                    Contig.start=min(start_align)
                    Contig.end=max(end_align)
                    Contig.misstart=start_align
                    Contig.misend=end_align
                    Contig.alignlength=alignlength
                    if start_contig[0]<end_contig[0]:
                        Contig.forward=1
                    else:
                        Contig.forward=2
                    contig_dict[processid]=Contig
                    start_align=[]
                    end_align=[]
                    start_contig=[]
                    end_contig=[]
                    alignlength=0
                    processid=int(data[5])
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3])) 
            count=count+1
        contig=aligncontig()
        Contig.start=min(start_align)
        Contig.end=max(end_align)
        Contig.misstart=start_align
        Contig.misend=end_align
        if start_contig[0]<end_contig[0]:
               Contig.forward=1
        else:
               Contig.forward=2
        Contig.alignlength=alignlength
        contig_dict[processid]=Contig
        return contig_dict


# Input contig alignment breaked by misassembly

def alignment_miss_contig(tsv_file,chr_id,chr_len):
        f_tsv = open(tsv_file,'r')
        count = 0
        contig_dict = defaultdict(aligncontig)
        start_align=[]
        end_align=[]
        start_contig=[]
        end_contig=[]
        alignlength=0
        processid=0
        contig_status=1
        indexid=0
        offset_array=[]
        offset_array.append(0)
        offset=0
        for length_chr in chr_len:
            offset=offset+length_chr
            offset_array.append(offset)
        for line in f_tsv:
            if count==0:
                count+=1
                continue
            elif count%2==0:
                data=line.rsplit()
                if data[0]=='relocation,':
                    contig_status=0            
            elif count%2 == 1:
                data=line.rsplit()
                if len(data)==4:
                    continue
                index_offset=chr_id.index(data[4])
                offset_add=offset_array[index_offset]
                if processid==0:
                    processid=int(data[5])
                    start_align.append(int(data[0])+offset_add)
                    end_align.append(int(data[1])+offset_add)
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
                    alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                elif processid==int(data[5]) and contig_status==1:
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
                else:
                    contig_status=1
                    Contig=aligncontig()
                    Contig.start=min(start_align)
                    Contig.end=max(end_align)
                    Contig.alignlength=alignlength
                    if start_contig[0]<end_contig[0]:
                        Contig.forward=1
                    else:
                        Contig.forward=2
                    Contig.misstart=start_align
                    Contig.misend=end_align
                    contig_dict[count]=Contig
                    start_align=[]
                    end_align=[]
                    start_contig=[]
                    end_contig=[]
                    alignlength=0
                    processid=int(data[5])
                    if int(data[0])<=int(data[1]):
                        start_align.append(int(data[0])+offset_add)
                        end_align.append(int(data[1])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    else:
                        start_align.append(int(data[1])+offset_add)
                        end_align.append(int(data[0])+offset_add)
                        alignlength=alignlength+abs(int(data[2])-int(data[3]))+1
                    start_contig.append(int(data[2]))
                    end_contig.append(int(data[3]))
            count=count+1
        contig=aligncontig()
        Contig.start=min(start_align)
        Contig.end=max(end_align)
        if start_contig[0]<end_contig[0]:
               Contig.forward=1
        else:
               Contig.forward=2
        Contig.misstart=start_align
        Contig.misend=end_align
        Contig.alignlength=alignlength
        contig_dict[processid]=Contig
        return contig_dict


# ## Calculate Contig Statistics
def CalculateNx_contig_align(contig_dict):
    NX_align_contig=[]
    contig_align_list=[]
    contig_align_list_new=[]
    for k,v in contig_dict.items():
       contig_align_list.append(v.alignlength)
    lastlength=0
    contig_align_list_new=sorted(contig_align_list,reverse=True)
    totallen=sum(contig_align_list_new)
    accu_len=0
    for contiglen in contig_align_list_new:
        lastlength=accu_len
        accu_len=accu_len+contiglen
        if accu_len>=0.1*totallen and lastlength<0.1*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.2*totallen and lastlength<0.2*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.3*totallen and lastlength<0.3*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.4*totallen and lastlength<0.4*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.5*totallen and lastlength<0.5*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.6*totallen and lastlength<0.6*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.7*totallen and lastlength<0.7*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.8*totallen and lastlength<0.8*totallen:
            NX_align_contig.append(contiglen/1000)
        if accu_len>=0.9*totallen and lastlength<0.9*totallen:
            NX_align_contig.append(contiglen/1000)
    return NX_align_contig

def CalculateNx_contig(infile,threshold):
    NX_contig=[]
    lastlength=0
    if '.gz' in infile:
         contig_fasta=gzip.open(infile,"r")
    else:
        contig_fasta=open(infile,"r")
    contig_list=[]
    contigseq=0
    index=0
    for contigx in contig_fasta:
        if '.gz' in infile:
            contig=contigx.decode()
        else:
            contig=contigx 
        if contig[0]!='>':
            contig_clean=contig.strip('\n')
            contigseq=contigseq+len(contig_clean)
        else:
            if contigseq>threshold*1000 and index!=0:
                contig_list.append(contigseq)    
            contigseq=0
        index+=1
    if contigseq>threshold*1000:
        contig_list.append(contigseq)
    contig_list_new=sorted(contig_list,reverse=True)
    totallen=sum(contig_list_new)
    accu_len=0
    for contiglen in contig_list_new:
        lastlength=accu_len
        accu_len=accu_len+contiglen
        if accu_len>=0.1*totallen and lastlength<0.1*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.2*totallen and lastlength<0.2*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.3*totallen and lastlength<0.3*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.4*totallen and lastlength<0.4*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.5*totallen and lastlength<0.5*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.6*totallen and lastlength<0.6*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.7*totallen and lastlength<0.7*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.8*totallen and lastlength<0.8*totallen:
            NX_contig.append(contiglen/1000)
        if accu_len>=0.9*totallen and lastlength<0.9*totallen:
            NX_contig.append(contiglen/1000)
    return NX_contig
def CalculateNAx_contig(misassemblist):
    lengthlist=[]
    NAX_contig=[]
    lastlength=0
    for key,value in misassemblist.items():
        lengthlist.append(value.alignlength)
    AA=sorted(lengthlist,reverse=True)
    totallength=sum(AA)
    length_agg=0
    for line in AA:
        lastlength=length_agg
        length_agg=length_agg+line
        if length_agg>=0.1*totallength and lastlength<0.1*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.2*totallength and lastlength<0.2*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.3*totallength and lastlength<0.3*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.4*totallength and lastlength<0.4*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.5*totallength and lastlength<0.5*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.6*totallength and lastlength<0.6*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.7*totallength and lastlength<0.7*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.8*totallength and lastlength<0.8*totallength:
            NAX_contig.append(line/1000)
        if length_agg>=0.9*totallength and lastlength<0.9*totallength:
            NAX_contig.append(line/1000)
    return NAX_contig
# ## Get chrid and chrlen

class contigcov:#class to store contig length and coverage,
        def __init__(self,length,coverage):
            self.length=length
            self.coverage=coverage
class uncovered:#class for uncovered genomic regions,
        def __init__(self, start, end, length):
            self.start=start
            self.end=end
            self.length=length


# In[20]:

def CalculateNCx_contig(alignment,ref_len):
        list_contig_cov=[]
        uncovered_region=[]
        NCx=[]
        nowlength=0
        lastlength=0
        laststart=0
        lastend=0
        index=0
        ref_len=ref_len*1000000
        total_uncover=ref_len
        total_uncoverold=0
        alignlen=[]
        init=uncovered(1,ref_len,ref_len)
        uncovered_region.append(init)
        for key,value in alignment.items():
            alignlen.append(value)
        alignlen.sort(key=operator.attrgetter('alignlength'),reverse=True)
        index=0
        for contig in alignlen:
            index+=1
            start=contig.misstart
            end=contig.misend
            if contig.start==laststart and contig.end==lastend:
                continue
            else:
                laststart=contig.start
                lastend=contig.end
            total_uncoverold=total_uncover
            for i in range(len(start)):
                total_uncover=0
                uncovered_regionold=uncovered_region
                uncovered_region=[]
                for region in uncovered_regionold: 
                    if region.start>=end[i]:
                        uncovered_region.append(region)
                        total_uncover+=region.length
                    elif region.end<=start[i]:
                        uncovered_region.append(region)                                                                                                           
                        total_uncover+=region.length
                    elif region.start<=end[i] and region.end>=end[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            length2=region.end-end[i]
                            unregion1=uncovered(region.start,start[i]-1,length1)
                            unregion2=uncovered(end[i]+1,region.end,length2)
                            uncovered_region.append(unregion1)
                            uncovered_region.append(unregion2)
                            total_uncover=total_uncover+length1+length2
                        elif start[i]<=region.start:
                            length1=region.end-end[i]
                            unregion=uncovered(end[i]+1,region.end,length1)   
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
                    elif region.end<=end[i] and region.end>=start[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            unregion=uncovered(region.start,start[i]-1,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
            lastlength=ref_len-total_uncoverold
            nowlength=ref_len-total_uncover 
            list_contig_cov.append(contigcov(contig.alignlength/1000,nowlength/ref_len))
            if lastlength<ref_len*0.1 and nowlength>=ref_len*0.1:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.2 and nowlength>=ref_len*0.2:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.3 and nowlength>=ref_len*0.3:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.4 and nowlength>=ref_len*0.4:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.5 and nowlength>=ref_len*0.5:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.6 and nowlength>=ref_len*0.6:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.7 and nowlength>=ref_len*0.7:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.8 and nowlength>=ref_len*0.8:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len*0.9 and nowlength>=ref_len*0.9:
                NCx.append(contig.alignlength/1000)
            if lastlength<ref_len and nowlength>=ref_len:
                NCx.append(contig.alignlength/1000)
        return NCx,list_contig_cov,round(nowlength/ref_len,4)


# ## Input Scaffold alignment

# In[21]:

class scaffolding(object):
          def __init__(self):
            self.scaffold_len=0
            self.contig_num=0
            self.pro_correct=0
            self.pro_error=0
            self.pro_skip=0
            self.contigalign_start=[]
            self.contigalign_end=[]
def analyze_scaffold_info(info_file,align_contig):
        f_info = open(info_file,"r")
        scaffold_dict1 = defaultdict(list) #contig id,
        scaffold_dict2 = defaultdict(list) #the order of contigs in scaffold based on assembly,
        scaffold_dict3 = defaultdict(list) #contig end position,
        scaffold_dict4 = defaultdict(list) #the order of contigs in scaffold based on alignment,
        scaffold_dict5 = defaultdict(list) #contig alignment length,
        scaffold_dict6 = defaultdict(list) #contig start list,
        scaffold_dict7 = defaultdict(list) #contig end list,
        total=0 
        for line in f_info:
            data = line.rsplit()
            scaffold_id = int(data[0])
            contig_id = int(data[1])
            contig_order = int(data[2])
            if contig_id in list(align_contig.keys()):
                scaffold_dict1[scaffold_id].append(contig_id)
                scaffold_dict2[scaffold_id].append(contig_order)
                scaffold_dict3[scaffold_id].append(align_contig[contig_id].end)
                scaffold_dict5[scaffold_id].append(align_contig[contig_id].alignlength)
                scaffold_dict6[scaffold_id].extend(align_contig[contig_id].misstart)
                scaffold_dict7[scaffold_id].extend(align_contig[contig_id].misend)
        scaffoldlist=[]
        for key,value in scaffold_dict3.items():
            if len(value) > 1:
                count_correct = 0
                count_err = 0
                count_skip = 0
                scaffold_dict4[key] = sorted(value)
                scaffold_info=scaffolding()
                scaffold_info.contig_num=len(scaffold_dict3[key])
                scaffold_info.contigalign_start=scaffold_dict6[key]
                scaffold_info.contigalign_end=scaffold_dict7[key]
                for i in range(len(scaffold_dict3[key])-1):
                    contig_id_a = scaffold_dict1[key][i]
                    contig_id_b = scaffold_dict1[key][i+1]
                    a = scaffold_dict3[key][i]
                    b = scaffold_dict3[key][i+1]
                    idx_a = scaffold_dict4[key].index(a)
                    idx_b = scaffold_dict4[key].index(b)
                    if idx_a + 1 ==  idx_b and align_contig[contig_id_a].forward ==1 and align_contig[contig_id_b].forward==1:
                        count_correct += 1
                    elif idx_a -1 ==idx_b and align_contig[contig_id_a].forward ==2 and align_contig[contig_id_b].forward==2:
                        count_correct += 1
                    elif idx_a +1 < idx_b and align_contig[contig_id_a].forward ==1 and align_contig[contig_id_b].forward==1:
                        count_skip += 1
                    elif idx_a -1 > idx_b and align_contig[contig_id_a].forward ==2 and align_contig[contig_id_b].forward==2:
                        count_skip += 1
                    else:
                        count_err += 1    
                    scaffold_info.scaffold_len=scaffold_info.scaffold_len+scaffold_dict5[key][i]/1000000
                scaffold_info.scaffold_len=scaffold_info.scaffold_len+scaffold_dict5[key][-1]/1000000
                scaffold_info.pro_correct=count_correct*100/(scaffold_info.contig_num-1)
                scaffold_info.pro_error=count_err*100/(scaffold_info.contig_num-1)
                scaffold_info.pro_skip=count_skip*100/(scaffold_info.contig_num-1)
                total+=scaffold_info.scaffold_len
            else:
                scaffold_info=scaffolding()
                scaffold_info.contigalign_start=scaffold_dict6[key]
                scaffold_info.contigalign_end=scaffold_dict7[key]
                scaffold_info.contig_num=len(scaffold_dict3[key])
                scaffold_info.scaffold_len=scaffold_info.scaffold_len+scaffold_dict5[key][-1]/1000000
                total+=scaffold_info.scaffold_len
                scaffold_info.pro_correct=-1
                scaffold_info.pro_error=-1
                scaffold_info.pro_skip=-1
            scaffoldlist.append(scaffold_info)
        return scaffoldlist


# ### block basic information (break scaffold by misscaffolding)

# In[22]:

class block(object):
          def __init__(self):
            self.block_num=0
            #self.block_len=0,
            self.block_alignlen=0
            self.block_startlist=[]
            self.block_endlist=[]
def analyze_block_info(info_file,alignment):
        f_info = open(info_file,"r")
        scaffold_dict1 = defaultdict(list)
        scaffold_dict2 = defaultdict(list)
        scaffold_dict3 = defaultdict(list)
        scaffold_dict4 = defaultdict(list)
        scaffold_dict5 = defaultdict(list)
        scaffold_dict6 = defaultdict(list)
        scaffold_dict7 = defaultdict(list)
        #total=0
        for line in f_info:
            data = line.rsplit()
            scaffold_id = int(data[0])
            contig_id = int(data[1])
            if contig_id in list(alignment.keys()):
                contig_order = int(data[2])
                scaffold_dict1[scaffold_id].append(contig_id)
                scaffold_dict2[scaffold_id].append(contig_order)
                scaffold_dict3[scaffold_id].append(alignment[contig_id].end)
                scaffold_dict5[scaffold_id].append(alignment[contig_id].alignlength)
                scaffold_dict6[scaffold_id].extend(alignment[contig_id].misstart)
                scaffold_dict7[scaffold_id].extend(alignment[contig_id].misend)
        block_list=[]
        for key,value in scaffold_dict3.items():
                count_correct = 1
                contig_id_a = scaffold_dict1[key][0]
                block_len=alignment[contig_id_a].alignlength/1000000
                block_startlist=alignment[contig_id_a].misstart[:]
                block_endlist=alignment[contig_id_a].misend[:]
                scaffold_dict4[key] = sorted(value)
                for i in range(len(scaffold_dict3[key])-1):
                    contig_id_a = scaffold_dict1[key][i]
                    contig_id_b = scaffold_dict1[key][i+1]
                    a = scaffold_dict3[key][i]
                    b = scaffold_dict3[key][i+1]
                    idx_a = scaffold_dict4[key].index(a)
                    idx_b = scaffold_dict4[key].index(b)
                    if idx_a + 1 ==  idx_b and alignment[contig_id_a].forward == 1 and alignment[contig_id_b].forward==1:
                        count_correct += 1
                        block_len+=alignment[contig_id_b].alignlength/1000000
                        block_startlist.extend(alignment[contig_id_b].misstart)
                        block_endlist.extend(alignment[contig_id_b].misend)
                    elif idx_a - 1 ==  idx_b and alignment[contig_id_a].forward == 2 and alignment[contig_id_b].forward==2:
                        count_correct += 1
                        block_len+=alignment[contig_id_b].alignlength/1000000
                        block_startlist.extend(alignment[contig_id_b].misstart)
                        block_endlist.extend(alignment[contig_id_b].misend)
                    else:
                        BL=block()
                        BL.block_num=count_correct
                        BL.block_alignlen=block_len
                        BL.block_startlist=block_startlist
                        BL.block_endlist=block_endlist   
                        block_list.append(BL)
                        if contig_id_b in alignment:
                            count_correct=1
                            block_len=alignment[contig_id_b].alignlength/1000000
                            block_startlist=alignment[contig_id_b].misstart[:]
                            block_endlist=alignment[contig_id_b].misend[:]
                BL=block()
                BL.block_num=count_correct
                BL.block_alignlen=block_len
                BL.block_startlist=block_startlist
                BL.block_endlist=block_endlist  
                block_list.append(BL)
        return block_list


# ## Evaluation for assemble scaffold

# ### Calculate NX for Block (Mb)

def CalculateNx_block(blocklist):
        Nx=[]
        total_length=0
        subtotal_length=0
        blocklist.sort(key=operator.attrgetter('block_alignlen'),reverse=True)
        for block in blocklist:
            total_length+=block.block_alignlen
        for block in blocklist:
            subtotal_length+=block.block_alignlen
            temp=subtotal_length-block.block_alignlen
            if temp<total_length*0.1 and subtotal_length>=total_length*0.1:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.2 and subtotal_length>=total_length*0.2:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.3 and subtotal_length>=total_length*0.3:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.4 and subtotal_length>=total_length*0.4:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.5 and subtotal_length>=total_length*0.5:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.6 and subtotal_length>=total_length*0.6:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.7 and subtotal_length>=total_length*0.7:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.8 and subtotal_length>=total_length*0.8:
                Nx.append(block.block_alignlen)
            if temp<total_length*0.9 and subtotal_length>=total_length*0.9:
                Nx.append(block.block_alignlen)
            if temp<total_length and subtotal_length>=total_length:
                Nx.append(block.block_alignlen)
        return Nx


# ### Calculate NX for Scaffold (Mb)

# In[24]:

def CalculateNx_scaffold_align(scaffoldlist):
        Nx=[]
        total_length=0
        subtotal_length=0
        scaffoldlist.sort(key=operator.attrgetter('scaffold_len'),reverse=True)
        #print(len(scaffoldlist))
        for scaffold in scaffoldlist:
            total_length+=scaffold.scaffold_len
        for scaffold in scaffoldlist:
            subtotal_length+=scaffold.scaffold_len
            temp=subtotal_length-scaffold.scaffold_len
            if temp<total_length*0.1 and subtotal_length>=total_length*0.1:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.2 and subtotal_length>=total_length*0.2:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.3 and subtotal_length>=total_length*0.3:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.4 and subtotal_length>=total_length*0.4:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.5 and subtotal_length>=total_length*0.5:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.6 and subtotal_length>=total_length*0.6:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.7 and subtotal_length>=total_length*0.7:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.8 and subtotal_length>=total_length*0.8:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length*0.9 and subtotal_length>=total_length*0.9:
                Nx.append(scaffold.scaffold_len)
            if temp<total_length and subtotal_length>=total_length:
                Nx.append(scaffold.scaffold_len)
        return Nx


# ### Calculate NGX for Scaffold

# In[25]:

def CalculateNGx_scaffold(scaffoldlist,ref_len):
    NGx=[]
    subtotal_length=0
    scaffoldlist.sort(key=operator.attrgetter('scaffold_len'),reverse=True)
    for scaffold in scaffoldlist:
        subtotal_length+=scaffold.scaffold_len
        temp=subtotal_length-scaffold.scaffold_len
        if temp<ref_len*0.1 and subtotal_length>=ref_len*0.1:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.2 and subtotal_length>=ref_len*0.2:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.3 and subtotal_length>=ref_len*0.3:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.4 and subtotal_length>=ref_len*0.4:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.5 and subtotal_length>=ref_len*0.5:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.6 and subtotal_length>=ref_len*0.6:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.7 and subtotal_length>=ref_len*0.7:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.8 and subtotal_length>=ref_len*0.8:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len*0.9 and subtotal_length>=ref_len*0.9:
            NGx.append(scaffold.scaffold_len)
        if temp<ref_len and subtotal_length>=ref_len:
            NGx.append(scaffold.scaffold_len)
    return NGx


# ### Calculate NGX for block

# In[26]:

def CalculateNGx_block(blocklist,ref_len):
        NGx=[]
        total_length=0
        subtotal_length=0
        blocklist.sort(key=operator.attrgetter('block_alignlen'),reverse=True)
        for block in blocklist:
            subtotal_length+=block.block_alignlen
            temp=subtotal_length-block.block_alignlen
            if temp<ref_len*0.1 and subtotal_length>=ref_len*0.1:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.2 and subtotal_length>=ref_len*0.2:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.3 and subtotal_length>=ref_len*0.3:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.4 and subtotal_length>=ref_len*0.4:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.5 and subtotal_length>=ref_len*0.5:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.6 and subtotal_length>=ref_len*0.6:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.7 and subtotal_length>=ref_len*0.7:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.8 and subtotal_length>=ref_len*0.8:
                NGx.append(block.block_alignlen)
            if temp<ref_len*0.9 and subtotal_length>=ref_len*0.9:
                NGx.append(block.block_alignlen)
            if temp<ref_len and subtotal_length>=ref_len:
                NGx.append(block.block_alignlen)
        return NGx


# ### Calculate NCX for scaffold

# In[27]:

def CalculateNCx_scaffold(align_scaffold,ref_len):
        list_contig_cov=[]
        uncovered_region=[]
        NCx=[]
        nowlength=0
        lastlength=0
        laststart=0
        lastend=0
        index=0
        ref_len=ref_len*1000000
        total_uncover=ref_len
        total_uncoverold=0
        init=uncovered(1,ref_len,ref_len)
        uncovered_region.append(init)
        align_scaffold.sort(key=operator.attrgetter('scaffold_len'),reverse=True)
        index=0
        A=0
        for scaffold in align_scaffold:
            index+=1
            start=scaffold.contigalign_start
            end=scaffold.contigalign_end
            total_uncoverold=total_uncover
            for i in range(len(start)):
                total_uncover=0
                uncovered_regionold=uncovered_region
                uncovered_region=[]
                for region in uncovered_regionold:
                    if region.start>=end[i]:
                        uncovered_region.append(region)
                        total_uncover+=region.length
                    elif region.end<=start[i]:
                        uncovered_region.append(region)    
                        total_uncover+=region.length
                    elif region.start<=end[i] and region.end>=end[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            length2=region.end-end[i]
                            unregion1=uncovered(region.start,start[i]-1,length1)
                            unregion2=uncovered(end[i]+1,region.end,length2)
                            uncovered_region.append(unregion1)
                            uncovered_region.append(unregion2)
                            total_uncover=total_uncover+length1+length2
                        elif start[i]<=region.start:
                            length1=region.end-end[i]
                            unregion=uncovered(end[i]+1,region.end,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
                    elif region.end<=end[i] and region.end>=start[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            unregion=uncovered(region.start,start[i]-1,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
            lastlength=ref_len-total_uncoverold
            nowlength=ref_len-total_uncover
            list_contig_cov.append(contigcov(scaffold.scaffold_len,nowlength/ref_len))
            if lastlength<ref_len*0.1 and nowlength>=ref_len*0.1:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.2 and nowlength>=ref_len*0.2:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.3 and nowlength>=ref_len*0.3:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.4 and nowlength>=ref_len*0.4:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.5 and nowlength>=ref_len*0.5:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.6 and nowlength>=ref_len*0.6:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.7 and nowlength>=ref_len*0.7:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.8 and nowlength>=ref_len*0.8:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len*0.9 and nowlength>=ref_len*0.9:
                NCx.append(scaffold.scaffold_len)
            if lastlength<ref_len and nowlength>=ref_len:
                NCx.append(scaffold.scaffold_len)
        return NCx,list_contig_cov


# ### Calculate NCX for block

# In[28]:

def CalculateNCx_block(align_block,ref_len):
        list_contig_cov=[]
        uncovered_region=[]
        NCx=[]
        nowlength=0
        lastlength=0
        laststart=0
        lastend=0
        index=0
        ref_len=ref_len*1000000
        total_uncover=ref_len
        total_uncoverold=0
        init=uncovered(1,ref_len,ref_len)
        uncovered_region.append(init)
        align_block.sort(key=operator.attrgetter('block_alignlen'),reverse=True)
        index=0
        A=0
        for block in align_block:
            index+=1
            start=block.block_startlist
            end=block.block_endlist
            total_uncoverold=total_uncover
            for i in range(len(start)):
                total_uncover=0
                uncovered_regionold=uncovered_region
                uncovered_region=[]
                for region in uncovered_regionold:
                    if region.start>=end[i]:
                        uncovered_region.append(region)
                        total_uncover+=region.length
                    elif region.end<=start[i]:
                        uncovered_region.append(region)  
                        total_uncover+=region.length
                    elif region.start<=end[i] and region.end>=end[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            length2=region.end-end[i]
                            unregion1=uncovered(region.start,start[i]-1,length1)
                            unregion2=uncovered(end[i]+1,region.end,length2)
                            uncovered_region.append(unregion1)
                            uncovered_region.append(unregion2)
                            total_uncover=total_uncover+length1+length2
                        elif start[i]<=region.start:
                            length1=region.end-end[i]
                            unregion=uncovered(end[i]+1,region.end,length1)   
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
                    elif region.end<=end[i] and region.end>=start[i]:
                        if start[i]>=region.start:
                            length1=start[i]-region.start
                            unregion=uncovered(region.start,start[i]-1,length1)
                            uncovered_region.append(unregion)
                            total_uncover=total_uncover+length1
            lastlength=ref_len-total_uncoverold
            nowlength=ref_len-total_uncover 
            list_contig_cov.append(contigcov(block.block_alignlen,nowlength/ref_len))
            if lastlength<ref_len*0.1 and nowlength>=ref_len*0.1:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.2 and nowlength>=ref_len*0.2:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.3 and nowlength>=ref_len*0.3:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.4 and nowlength>=ref_len*0.4:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.5 and nowlength>=ref_len*0.5:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.6 and nowlength>=ref_len*0.6:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.7 and nowlength>=ref_len*0.7:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.8 and nowlength>=ref_len*0.8:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len*0.9 and nowlength>=ref_len*0.9:
                NCx.append(block.block_alignlen)
            if lastlength<ref_len and nowlength>=ref_len:
                NCx.append(block.block_alignlen)
        return NCx,list_contig_cov
## reference id and length
def reflen_id(in_path_ref):
    gzornot=0
    if "gz" in in_path_ref:
        infile=gzip.open(in_path_ref,'r')
        gzornot=1
    else:
        infile=open(in_path_ref,'r')
    chrid=[]
    chrlen=[]
    singlen=0
    for line in infile:
        lineinfo=""
        if gzornot==1:
           lineinfo=line.decode()
        else:
           lineinfo=line
        if lineinfo[0]=='>':
            A_1=lineinfo[1:].replace('  ','__').strip('\n')
            A_2=A_1.replace(' ','_')
            B=A_2.replace(':','_')
            chrid.append(B)
            if singlen!=0:
                chrlen.append(singlen)
                singlen=0
        else:
            singlen=singlen+len(lineinfo[:])
    chrlen.append(singlen)
    return chrid,chrlen,sum(chrlen)/1000000
def helpinfo():
    helpinfo=\
    '''
        Calculate NX,NAX,NCX and NCAX for contigs and scaffolds
        Version: 1.0.0
        Dependents: Python (>=3.0)
        Last Updated Date: 2017-07-26
        Contact: zhanglu295@gmail.com

        Usage: python Cal_Denovo_Stat <options>

        Options:
                -t --tsv, tsv files located in ./contigs_reports/all_alignments_*.tsv from QUAST(multiple inputs separated by comma)
                -c --contig, compressed or uncompressed contig fasta files (multiple inputs separated by comma, only used to calculate Contig NX)
                -s --scaffold, conpressed or uncompressed scaffold fasta files (multiple inputs separated by comma, only used to calculate Scaffold NX)
                -r --reference, compressed or uncompressed reference genome fasta, if reference is missing, only alignment free statistics would be caluclated.
                -i --info, .info file(three columns: scaffold id, contig id, order of contig in the scaffold)
                -a --min_contig, the length of minimum contigs (kb)(only for calculating N50, the short contigs for calculating NA50,NC50 and NCA50 should be eliminated first by -m in QUAST)
                -b --min_scaffold, the length of minimum scaffold (kb)
                -p --prefix, the prefix of output file(default:Assembly_summary)
                -l --label, Human-readable assembly names (multiple labels separated by comma).

                -o --out, the path to output
                -h --help, help info
    '''
    print(helpinfo)
def main():
    tsvlist=[]
    contiglist=[]
    scaffoldlist=[]
    reference=None
    infolist=[]
    labellist=[]
    outpath=None
    min_contig=0
    min_scaffold=0
    prefix="Assembly_summary"
    opts, args = getopt.gnu_getopt(sys.argv[1:], 't:c:s:r:i:l:a:b:p:o:h', ['tsv', 'contig','scaffold' ,'reference', 'info' ,'label','min_contig','min_scaffold','prefix','out','help'])
    for o, a in opts:
        if o == '-t' or o == '--tsv':
                tsvlist=a.split(',') 
        if o == '-c' or o == '--contig':
                contiglist = a.split(',')
        if o == '-s' or o == '--scaffold':
                scaffoldlist = a.split(',')
        if o == '-r' or o == '--reference':
                reference = a
        if o == '-i' or o == '--info':
                infolist = a.split(',')
        if o == '-l' or o == '--label':
                labellist = a.split(',')
        if o =='-a' or o =='--min_contig':
                min_contig = float(a)
        if o =='-b' or o =='--min_scaffold':
                min_scaffold = float(a)
        if o =='-p' or o =='--prefix':
                prefix = a
        if o == '-o' or o == '--out':
                outpath = a
        if o == '-h' or o == '--help':
                helpinfo()
                sys.exit(-1)
    outfile=open(outpath+'/'+prefix+'.txt','w')
    if reference==None:
           outfile.write('Statistics')
           outfile.write('\t')
           if len(labellist)==0:
               for index in range(len(contiglist)):
                   outfile.write('Sample'+str(index+1)+'\t')
           else:
               for label in labellist:
                   outfile.write(label+'\t')
           outfile.write('\n')
           outfile.write('Contig N50(kb)'+'\t')
           for contig in contiglist:
               Contig_NX=CalculateNx_contig(contig,min_contig)
               outfile.write(str(round(Contig_NX[4],4))+'\t')
           outfile.write('\n')
           outfile.write('Scaffold N50(kb)'+'\t')
           for scaffold in scaffoldlist:
               Scaffold_NX=CalculateNx_contig(scaffold,min_scaffold)
               outfile.write(str(round(Scaffold_NX[4],4))+'\t')
           outfile.close()
    else:
           outfile.write('Statistics')
           outfile.write('\t')
           if len(labellist)==0:
               for index in range(len(contiglist)):
                    outfile.write('Sample'+str(index+1)+'\t')
           else:
               for label in labellist:
                   outfile.write(label+'\t')
           outfile.write('\n')
           print('Calculationg Contig N50...')
           outfile.write('Contig N50(kb)'+'\t')
           for contig in contiglist:
               Contig_NX=CalculateNx_contig(contig,min_contig)
               outfile.write(str(round(Contig_NX[4],4))+'\t')
           outfile.write('\n')
           print('Calculationg Scaffold N50...')
           outfile.write('Scaffold N50(kb)'+'\t')
           for scaffold in scaffoldlist:
               Scaffold_NX=CalculateNx_contig(scaffold,min_scaffold)
               outfile.write(str(round(Scaffold_NX[4],4))+'\t')
           outfile.write('\n')
           index=1
           print('Reading refernece genome...')
           (chrid,chrlen,singlen)=reflen_id(reference)
           list_align_contig=[]
           list_misassemble_contig=[]
           list_align_scaffold=[]
           list_align_block=[]
           print('Reading Alignment files...')
           for tsvid in tsvlist:
               align_contig=alignment_contig_coordinate(tsvid,chrid,chrlen)
               misassemble_contig=alignment_miss_contig(tsvid,chrid,chrlen)
               list_align_contig.append(align_contig)
               list_misassemble_contig.append(misassemble_contig)
           for infoid in range(len(infolist)):
               align_scaffold=analyze_scaffold_info(infolist[infoid],list_align_contig[infoid])
               align_block=analyze_block_info(infolist[infoid],list_align_contig[infoid])
               list_align_scaffold.append(align_scaffold)
               list_align_block.append(align_block)
           print('Calculating Contig N50 (Just keep segments can be aligned to the reference genome)')
           outfile.write('Contig N50 (alignment,kb)'+'\t')
           index=0
           for contigid in list_align_contig:
               NX_contig_align=CalculateNx_contig_align(contigid)
               outfile.write(str(round(NX_contig_align[4],4))+'\t')
           outfile.write('\n')
           outfile.write('Scaffold N50 (alignment, kb)'+'\t')
           print('Calculating Scaffold N50 (Just keep segments can be aligned to the reference genome)')
           for scaffoldid in list_align_scaffold:
               #print(len(scaffoldid))
               NX_scaffold_align=CalculateNx_scaffold_align(scaffoldid)
               outfile.write(str(round(NX_scaffold_align[4]*1000,4))+'\t')
           outfile.write('\n')
           print('Calculate Contig NA50...')
           outfile.write('Contig NA50(kb)'+'\t')
           for misassemble_contig in list_misassemble_contig:
               NA_contig=CalculateNAx_contig(misassemble_contig)
               outfile.write(str(round(NA_contig[4],4))+'\t')
           outfile.write('\n')
           print('Calculate Scaffold NA50...')
           outfile.write('Scaffold NA50(kb)'+'\t')
           for align_block in list_align_block:
               NX_block=CalculateNx_block(align_block)
               outfile.write(str(round(NX_block[4]*1000,4))+'\t')
           outfile.write('\n')
           coverage_list=[]
           print('Calculate Contig NC50...')
           outfile.write('Contig NC50(kb)'+'\t')
           #print(singlen)
           for contigid in list_align_contig:
               [NC_contig,list_NC_contig,coverage]=CalculateNCx_contig(contigid,singlen)
               coverage_list.append(coverage)
               if len(NC_contig)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NC_contig[4],4)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Scaffold NC50...')
           outfile.write('Scaffold NC50(kb)'+'\t')
           for scaffoldid in list_align_scaffold:
               [NC_scaffold,list_NC_scaffold]=CalculateNCx_scaffold(scaffoldid,singlen)
               if len(NC_scaffold)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NC_scaffold[4]*1000,4)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Contig NCA50...')
           outfile.write('Contig NCA50(kb)'+'\t')
           for misassemble_contig in list_misassemble_contig:
               [NCA_contig,list_NCA_contig,coverage]=CalculateNCx_contig(misassemble_contig,singlen)
               if len(NCA_contig)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NCA_contig[4],4)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Scaffold NCA50...')
           outfile.write('Scaffold NCA50(kb)'+'\t')
           for align_block in list_align_block: 
               [NC_block,list_NC_block]=CalculateNCx_block(align_block,singlen)
               if len(NC_block)<5:
                   outfile.write('0')
               else:
                   outfile.write(str(round(NC_block[4]*1000,4)))
               outfile.write('\t')
           outfile.write('\n')
           print('Calculate Scaffold genome coverage...')
           outfile.write('Genome Coverage(%)'+'\t')
           for coverage in coverage_list:
               outfile.write(str(round(coverage*100,2))+'%\t')
           outfile.close()
    return None
if __name__=="__main__":
    import sys
    if len(sys.argv) == 1:
         helpinfo()
    else:
        from collections import defaultdict
        #import matplotlib
        import re
        import numpy
        import operator
        from scipy import stats as sta
        #import seaborn as sns
        import math
        #from scipy.interpolate import spline
        from itertools import groupby
        import gzip
        import getopt
        main()

