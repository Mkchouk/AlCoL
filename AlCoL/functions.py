import os
import sys
import time
import operator

def file_to_list(line):
    line=line[:-1]
    L_elem=[]
    L_elem=line.split(",")
    L_elem[2]=int(L_elem[2])
    L_elem[3]=int(L_elem[3])
    L_elem[4]=float(L_elem[4])
    L_elem[5]=int(L_elem[5])
    L_elem[6]=int(L_elem[6])
    return L_elem

def reverse_liste_overlaps(Liste,index_del):
    index_del=index_del.reverse()

def del_overlaps(Liste):
    length=len(Liste)
    index_delete=[]
    for i in range(1,length):
        if Liste[i-1][2]<=Liste[i][2] and Liste[i-1][3]>=Liste[i][3]:
            index_delete.append(i)
        if Liste[i-1][2]>=Liste[i][2] and Liste[i-1][3]<=Liste[i][3]:
            index_delete.append(i-1)
    reverse_liste_overlaps(Liste,index_delete)
   # print("index_delete=",index_delete)
    Liste=delete_intern_overlap(Liste,index_delete)
    return Liste

def delete_intern_overlap(Liste_p,index_del_liste):
    for i in index_del_liste:
        del Liste_p[i]
    return Liste_p

def print_Liste(Liste_p):
    for list in Liste_p:
        print(list[:7])