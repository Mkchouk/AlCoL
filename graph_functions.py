from functions import *

def file_creation(Liste):
    graph={}
    length=len(Liste)
    file2 = open("graph_file.dat", "w")
    i=0
    j=1
    out_from_j=True
    while i<length:
        j=i+1
        while j<length:
            #print("i= "+str(i)+" and "+"j= "+str(j))
            if Liste[i][3]>=Liste[j][2]:
                weight=calcul_weight(Liste[i][4],Liste[i][6],Liste[i][5],Liste[j][4],Liste[j][6])
                #print(str(Liste[i][1])+" --> "+str(Liste[j][1])+" ==> "+str(weight))
                file2.write(str(Liste[i][1])+" "+str(Liste[j][1])+" "+str(weight)+'\n')
                j=j+1
            else:
                i=i+1
                j=length+1
        if j==length:
            break

def calcul_weight(match_percent_node1,alignment_len_node1,query_len,match_percent_node2,alignment_len_node2):
    P_match_node1=match_percent_node1/100
    P_len_node1=alignment_len_node1/query_len
    P_match_node2=match_percent_node2/100
    P_len_node2=alignment_len_node2/query_len
    weight_node1=(P_match_node1+P_len_node1)/2
    weight_node2=(P_match_node2+P_len_node2)/2
    return (weight_node1+weight_node2)/2

def Graph_creation(file_name):
    file=open(file_name, "r")
    dic_for_graph={}
    for line in file:
        line=line[:-1]
        L=[]
        L=line.split(" ")
        L[2]=float(L[2])
        #print(L)
        if L[0] not in dic_for_graph.keys():
            dic_for_graph[L[0]]=[[L[1],L[2]]]
            #print(str(L[0])+" added")
        else:
            dic_for_graph[L[0]].append([L[1],L[2]])
            dic_for_graph[L[0]]=sorted(dic_for_graph[L[0]], key=operator.itemgetter(1),reverse=True)
    return dic_for_graph

def Path(dic_graph_p,f_node):
    path=[]
    while True:
        if f_node in dic_graph_p.keys():
            #print("f_node= "+f_node)
            path.append(f_node)
            L=dic_graph_p[f_node]
            f_node=L[0][0]
        elif f_node not in dic_graph_p.keys():
            path.append(f_node)
            print("Path Found!!")
            break
        else:
            print("FATAL ERROR There is a Problem!!")
            break
    return path
    
def dic_for_contigs_creation(Liste_):
    dic_contigs={}
    i=1
    for contigs in Liste_:
        if contigs[1] not in dic_contigs.keys():
            dic_contigs[contigs[1]]=contigs
        else:
            i=i+1
            new_key=contigs[1]+"_"+str(i)
            contigs[1]=new_key
            dic_contigs[new_key]=contigs
            #print("new_key= ",new_key)
            #print("FATAL ERROR in function dic_for_contigs_creation in graph_functions.py file!!!")
    return dic_contigs

def Liste_from_contigs_path(path,dic_contigs):
    Liste_path=[]
    for contig in path:
        if contig in dic_contigs.keys():
            Liste_path.append(dic_contigs[contig])
    return Liste_path

def delete_internal_overlaps(Liste_):
    length=len(Liste_)
    index_delete=[]
    i=1
    index=0
    while i<length:
        if Liste_[index][2]<=Liste_[i][2] and Liste_[index][3]>=Liste_[i][3]:
            index_delete.append(i)
            i=i+1
        else:
            index=i
            i=i+1
    reverse_liste_overlaps(Liste_,index_delete)
    #print("index__internal_delete=",index_delete)
    Liste_=delete_intern_overlap(Liste_,index_delete)
    return Liste_
    
"""    for i in range(1,length):
        print("i= ",i)
        if Liste_[i-1][2]<=Liste_[i][2] and Liste_[i-1][3]>=Liste_[i][3]:
            index_delete.append(i)"""

def Graph_number(Liste_):
    length=len(Liste_)
    if length==0:
        print("La liste est vide")
        return []
    else:
        Contig_name_for_graph_begin=[]
        Contig_name_for_graph_begin.append(Liste_[0][1])
        for i in range(1,length):
            if Liste_[i-1][3] < Liste_[i][2]:
                Contig_name_for_graph_begin.append(Liste_[i][1])
        return Contig_name_for_graph_begin

