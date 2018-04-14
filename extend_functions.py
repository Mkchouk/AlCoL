from functions import *
def Extend(pos_begLR__,pos_endLR__,Liste_non_overlap__,path__,dic_contig__):
    if len(path__) ==1:
        contig_name=path__[0]
        contig_seq=dic_contig__[contig_name]
        contig_beg=int(Liste_non_overlap__[0][7])
        contig_end=int(Liste_non_overlap__[0][8])
        if contig_beg>contig_end: #to treat the reverse complement
            diff_beg=pos_begLR__-1
            diff_end=int(Liste_non_overlap__[0][5])-pos_endLR__
            new_contig_beg=contig_beg+diff_beg+int(diff_beg*15/100) #it's contig end
            new_contig_end=contig_end-(diff_end+int(diff_end*15/100)) #it's contig begin
            if new_contig_end<1: #because it the contig beg (RevComp)
                new_contig_end=1
            if new_contig_beg>len(contig_seq):
                new_contig_beg=len(contig_seq)
            consensus_=contig_seq[new_contig_end-1:new_contig_beg+1]
            consensus=ReverseComplement(consensus_)
            return consensus
        else:
            diff_beg=pos_begLR__-1
            diff_end=int(Liste_non_overlap__[0][5])-pos_endLR__
            new_contig_beg=contig_beg-(diff_beg+int(diff_beg*15/100))
            new_contig_end=contig_end+diff_end+int(diff_end*15/100)
            if new_contig_beg<1:
                new_contig_beg=1
            if new_contig_end>len(contig_seq):
                new_contig_end=len(contig_seq)
            consensus=contig_seq[new_contig_beg-1:new_contig_end+1]
            return consensus
    else:
        print("FATAL ERROR in function: _Extend_ from _extend_functions.py_ beacause the lenght of path is >1")
        sys.exit(0)


def Extend_two(pos_begLR_,pos_endLR_,Liste_non_overlap_,path_,dic_contig_):
    """Process first contig in path_ The process just the first position in contig and first position in LR"""
    contig_name=path_[0]
    print("contig_name= ",contig_name)
    if contig_name not in dic_contig_.keys():
        contig_name=test_contig_name(contig_name,dic_contig_)
    
    contig_seq=dic_contig_[contig_name]
    contig_beg=int(Liste_non_overlap_[0][7])
    contig_end=int(Liste_non_overlap_[0][8])
    if contig_beg>contig_end: #to treat the Reverse Complement
        diff_beg=pos_begLR_-1
        new_contig_beg=contig_beg+diff_beg+int(diff_beg*15/100)
        if new_contig_beg>len(contig_seq):
            new_contig_beg=len(contig_seq)
        consensus_beg=contig_seq[contig_beg-1:new_contig_beg+1]
        consensus_beg_=ReverseComplement(consensus_beg)#we'll put this at the beg of the corrected LR
    else: #If contig_beg<contig_end
        diff_beg=pos_begLR_-1
        new_contig_beg=contig_beg-(diff_beg+int(diff_beg*15/100))
        if new_contig_beg<1:
            new_contig_beg=1
        consensus_beg_=contig_seq[new_contig_beg-1:contig_beg+1]
    
    """---------------------------------------------------------------------------------------"""
    """Process last contig in path_ to process the last position in contig and right position in LR"""
    lenght_path=len(Liste_non_overlap_)
    contig_name=""
    contig_seq=""
    contig_name=path_[lenght_path-1]
    print("contig_name= ",contig_name)
    if contig_name not in dic_contig_.keys():
        contig_name=test_contig_name(contig_name,dic_contig_)
    
    contig_seq=dic_contig_[contig_name]
    contig_beg=int(Liste_non_overlap_[lenght_path-1][7])
    contig_end=int(Liste_non_overlap_[lenght_path-1][8])
    if contig_beg>contig_end: #to treat the Reverse Complement
        diff_end=int(Liste_non_overlap_[0][5])-pos_endLR_
        new_contig_end=contig_end-(diff_end+int(diff_end*15/100))
        if new_contig_end<1:
            new_contig_end=1
        consensus_end=contig_seq[new_contig_end-1:contig_end+1]
        consensus_end_=ReverseComplement(consensus_end)#we'll put this at the end of the corrected LR
    else:
        diff_end=int(Liste_non_overlap_[0][5])-pos_endLR_
        new_contig_end=contig_end+diff_end+int(diff_end*15/10)
        if new_contig_end>len(contig_seq):
            new_contig_end=len(contig_seq)
        consensus_end_=contig_seq[contig_end-1:new_contig_end+1]
    return consensus_beg_,consensus_end_

def LR_Consensus_extend(L_non_overlap,pos_begLR_,pos_endLR_,path_,dic_contig):
    length=len(L_non_overlap)
    consensus=""
    for i in range(1,length):
        nbr_nucl=L_non_overlap[i][2]-L_non_overlap[i-1][2]
        #print(nbr_nucl)
        cpt=0
        len_substring=0
        for j in L_non_overlap[i-1][9]:
            len_substring=len_substring+1
            if j != "-":
                cpt=cpt+1
            if cpt==nbr_nucl:
                break
        #print("len_substring= ",len_substring)
        for c_nucl in range(len_substring):
            consensus=consensus+L_non_overlap[i-1][10][c_nucl]
        #print("len consensus= ",len(consensus))
    consensus=consensus+L_non_overlap[length-1][10]
    print("Final len consensus= ",len(consensus))
    #HERE we add the function Extend_two(...)
    consensus_begin,consensus_ending=Extend_two(pos_begLR_,pos_endLR_,L_non_overlap,path_,dic_contig)
    final_consensus=consensus_begin+consensus+consensus_ending
    return final_consensus

def test_contig_name(contig__name,dic__contig_):#To test if contig_name in dic_contig or  not
    if contig__name not in dic__contig_.keys():
        contig__name=contig__name[:-2]
        if contig__name not in dic__contig_.keys():
            contig__name=contig__name[:-3]
            if contig__name not in dic__contig_.keys():
                contig__name=contig__name[:-4]
            else:
                print("Fatal Error: Contig_name not in dic_contig!!")
                sys.exit(0)
    return contig__name
