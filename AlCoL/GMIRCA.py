import os
import sys
import time
import operator
from functions import *
from graph_functions import *

def Find_pos_deb_end_LR(L_non_overlap):
    Pos_begLR=L_non_overlap[0][2]
    length=len(L_non_overlap)
    Pos_endLR=L_non_overlap[length-1][3]
    return Pos_begLR,Pos_endLR

def LR_Consensus(L_non_overlap):
    length=len(L_non_overlap)
    consensus=""
    for i in range(1,length):
        nbr_nucl=L_non_overlap[i][2]-L_non_overlap[i-1][2]
        #print(nbr_nucl)
        cpt=0
        len_substring=0
        for j in L_non_overlap[i-1][7]:
            len_substring=len_substring+1
            if j != "-":
                cpt=cpt+1
            if cpt==nbr_nucl:
                break
        #print("len_substring= ",len_substring)
        for c_nucl in range(len_substring):
            consensus=consensus+L_non_overlap[i-1][8][c_nucl]
        #print("len consensus= ",len(consensus))
    consensus=consensus+L_non_overlap[length-1][8]
    print("Final len consensus= ",len(consensus))
    return consensus

def dic_LR_LRfile(LR_file_name):
    dic_contig={}
    for line in LR_file_name:
        if line[0]=='>':
            dic_contig_id=line[1:-1]
        else:
            if line[len(line)-1]=='\n' or line[len(line)-1]=='\r':
                contig_seq=line[:-1]
            if contig_seq[len(contig_seq)-1]=='\n' or contig_seq[len(contig_seq)-1]=='\r':
                contig_seq=contig_seq[:-1]
            dic_contig[dic_contig_id]=contig_seq
    return dic_contig

def Consensus_LR_without_gaps(consensusLR):
    consensus_LR_wt_gaps=""
    for nucl in consensus:
        if nucl !="-":
            consensus_LR_wt_gaps=consensus_LR_wt_gaps+nucl
    return consensus_LR_wt_gaps


def Correct_LR(Pos_begLR,Pos_endLR,dic_LR_from_lrfile,lr_name,consensus_lr,dic_corrected_lr,last_pos_end):
    Long_read=dic_LR_from_lrfile[lr_name]
    print("coverage= ",(pos_endLR-pos_begLR+1)/len(dic_LR_from_LRfile[Liste_non_overlap[0][0]]))
    print("Len Long read= ",len(Long_read))
    #Corrected_LR=Long_read[:Pos_begLR]+consensus_lr+Long_read[Pos_endLR+1:]
    if lr_name in dic_corrected_lr.keys() and last_pos_end==-1:
        print("FATAL ERROR in Function Correct_LR!!!")
        sys.exit(0)

    if lr_name in dic_corrected_lr.keys():
        temp_Corrected_LR=dic_corrected_lr[lr_name]
        temp_Corrected_LR=temp_Corrected_LR+Long_read[last_pos_end:Pos_begLR]+consensus_lr
        return temp_Corrected_LR
    else:
        #Corrected_LR=Long_read[:Pos_begLR]+consensus_lr
        Corrected_LR=consensus_lr
        return Corrected_LR
#    print("Len CORRECTED Long read= ",len(Corrected_LR))
    

#************************ main Program**************************#
file1 = open("../../simulated_reads/cerevisae/resultatblastn_sim_cerevisae_without_Evalue_improving", "r")
LR_file = open("../../PBSIM/sim_cerevisae/sd_sim_cerevisae.fa", "r")
dic_LR_from_LRfile=dic_LR_LRfile(LR_file)
#print(dic_LR_from_LRfile)
precedent=""
Liste=[]
dic_contigs={}
dic_corrected_LR={}
for line1 in file1:
    L=[]
    L=file_to_list(line1)
    #print("L=",L)
    if len(Liste)==0:
        Liste.append(L)

    if len(Liste)>0 and Liste[0][0]==L[0]:
        Liste.append(L)
    else:
        
        Liste=sorted(Liste, key=operator.itemgetter(2))
        print("-------------------------------Processing of Long Read ",Liste[0][0],"------------------------------")
        #print("\nListe=", Liste)
        print("\n************************Sorted List********************\n")
        #print_Liste(Liste)
        Liste_r=delete_internal_overlaps(Liste)
        print("\n**************Sorted with deleting internal overlaps List**********\n")
        #print_Liste(Liste_r)
        contig_name_for_graph_begin=Graph_number(Liste_r)
        #print("Contig_begin= ",contig_name_for_graph_begin)
        dic_contigs=dic_for_contigs_creation(Liste_r)
        #print("dic_contigs= ",dic_contigs)
        file_creation(Liste_r)
        dic_graph=Graph_creation("graph_file.dat")
        #print("dic_graph= ",dic_graph)
        Last_pos_end=-1
        #determine the path using graph number
        for first_node in contig_name_for_graph_begin:
            print("**********************NEW PATH**********************")
            path=Path(dic_graph,first_node)
            print("path= ",path)
            Liste_path=[]
            Liste_path= Liste_from_contigs_path(path,dic_contigs)
            #print("Liste_path= ",Liste_path)
           # print("Liste_path")
          #  print_Liste(Liste_path)
            Liste_non_overlap=[]
            Liste_non_overlap=del_overlaps(Liste_path)
           # print("Liste_non_overlap")
           # print_Liste(Liste_non_overlap)
            
            pos_begLR,pos_endLR=Find_pos_deb_end_LR(Liste_non_overlap)
            print("pos_begLR= ",pos_begLR,"  pos_endLR= ",pos_endLR)
            consensus=LR_Consensus(Liste_non_overlap)
            consensus_LR=Consensus_LR_without_gaps(consensus)
            print("len consensus_LR without gaps= ",len(consensus_LR))
            LR_name=Liste_non_overlap[0][0]
            corrected_LR=Correct_LR(pos_begLR,pos_endLR,dic_LR_from_LRfile,LR_name,consensus_LR,dic_corrected_LR,Last_pos_end)
            
            if LR_name not in dic_corrected_LR.keys():
                dic_corrected_LR[LR_name]=corrected_LR
            else:
                del dic_corrected_LR[LR_name]
                dic_corrected_LR[LR_name]=corrected_LR
            Last_pos_end=pos_endLR
            print("Len CORRECTED Long read= ",len(dic_corrected_LR[LR_name]))
        
        Liste=[]
        Liste.append(L)

file2 = open("Corrected_LR.fasta", "w")
for lr_name in dic_corrected_LR:
    file2.write(">"+str(lr_name)+"\n")
    file2.write(str(dic_corrected_LR[lr_name])+"\n")
print("\nListe=", Liste)
#del_overlaps(Liste)
