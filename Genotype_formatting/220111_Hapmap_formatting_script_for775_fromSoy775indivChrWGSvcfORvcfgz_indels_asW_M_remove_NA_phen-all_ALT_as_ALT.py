# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 15:21:17 2020

@author: Host
"""

# For NIck Acutool vcf.gz or Soy775.vcf.gz file only!
# If ref>1 R assign in hap map 
# all/one for alt alleles

import gzip
import math
import re
from datetime import datetime 
import pandas as pd
import numpy as np

##Function vcf line to hhmp
def ln_to_hmp(line,i_indels="one"):
    variables = [line[3]]+line[4].split(",") 

    indel="n"
    for a in variables:
        if len(a)>1: indel="y"
    ln = line.copy()  
    if i_indels=="one" and len(variables)>2:
        d = dict ()
        for a in ln[9:]:
            if "0" not in a[:3] and "*" not in a[:3] and "." not in a[:3] and a[0]==a[2]:
                d[a[:3]]=d.get(a[:3],0)+1
        try:
            max_key = max(d, key=d.get)
            max_key_pos = int(max_key[0])
        except:
            max_key_pos=1
            print(d)
        d_var = dict()
       
        if indel=="y":
            new_variables = [("N") for a in variables]
            if len(variables[0])>1: new_variables[0]="R"
            else: new_variables[0]=variables[0]
            if len(variables[0])>len(variables[max_key_pos]):new_variables[max_key_pos]="M"
            elif len(variables[0])<len(variables[max_key_pos]):new_variables[max_key_pos]="W"
            else:
                if len(variables[max_key_pos])==1:new_variables[max_key_pos]=variables[max_key_pos]
                else: new_variables[max_key_pos]="W"
            for a,n in zip(new_variables,range(0,len(new_variables))): d_var[n]=a
        else: 
            for a,n in zip(variables,range(0,len(variables))): d_var[n]=a
        line=ln[:9]
        
    else:
        d_var = dict()
        if indel=="y":
            new_variables = []
            if len(variables[0])>1: new_variables.append("R")
            else: new_variables.append(variables[0])
            for a in variables[1:]:
                if len(variables[0])>len(a):new_variables.append("M")
                elif len(variables[0])<len(a):new_variables.append("W")
                else:
                    if len(a)==1:new_variables.append(a)
                    else: new_variables.append("W")
            for a,n in zip(new_variables,range(0,len(new_variables))): d_var[n]=a
        else: 
            for a,n in zip(variables,range(0,len(variables))): d_var[n]=a
        
    line=ln[:9]
    for a in ln[9:]:
            if a[0]==a[2] and "." not in a[:3] and "*" not in a[:3]: 
                try: line.append(d_var[int(a[0])])
                except: 
                    print(a[0],a[:3],d_var)
                    line.append("N")
                
            else: 
                line.append("N")
                
   
    return line

#Assigning the chromosome and positions values
i_POS = input ("Please, assign a genomic region of your desired hapmap." 
             "Use following format: chromosome,start position,end position (e.g. 2,105,900), use only digits. ", )
pos = i_POS.split(",")

i_vcfgz = input("Please, assign a genome file in a vcf or vcf.gz format." )
i_indels = input("Do you want to assign as Alt all INDELs or only the most frequent one?(all/one)." )

#Chromosome and positions values check and completing
if pos[0].isnumeric () is False:
    print ("Wrong chromosome value, use only digit!")
    raise Exception("exit")

pos[1]=0 if pos[1].isnumeric() is False else int(pos[1])
pos[2]=math.inf if pos[2].isnumeric() is False else int(pos[2])


res_name = input ("Name your hapmap. Do not use a file extension.",)
fres = input ("Do you want to incorporate date and time in formate YYMMDDHHMM (e.g. 16:45 8.10.2020 - 2010081645) into the file result name? (y/n)", )
now  = datetime.now()

x=0
lst_PI = list()
i_phn = input ("Do you want to remove isolines without phenotype? (y/n)",)
if i_phn =="y":
    phn_nam = input ("Assign the phenotype file.",)
    i_pn_f = input ("Do you want to create new phenotype file without missing values? (y/n)",)
    if i_pn_f=="y": 
        pn_f_nam = input ("Name your phenotype file. Do not use a file extension.",)
        pn_f = open(pn_f_nam+".txt","w")
    f = open (phn_nam, "rt")
    for line in f:
        if x==0: 
            x=1
            continue
        line=line.split()
        
        if line[1]!="NA" and line[1]!="-" and line[1]!="." and line[1]!="" and line[1]!=" ":
            if i_pn_f=="y":
                pn_f.write(line[0])
                pn_f.write("\t")
                pn_f.write(line[1])
                pn_f.write("\n")
            lst_PI.append(line[0])
            if "PI" or "FC" in line[0]:
                if re.search("PI\d+_[0-9]",line[0]): 
                    a=line[0].replace("_", "-")
                    lst_PI.append(a.replace("PI", "PI_"))
                elif re.search("PI_\d+-[0-9]",line[0]):  
                    a=line[0].replace("PI_", "PI")
                    lst_PI.append(a.replace("-", "_"))
                elif re.search("PI\d+",line[0]):  lst_PI.append(line[0].replace("PI", "PI_"))
                elif re.search("PI_\d+",line[0]):  lst_PI.append(line[0].replace("PI_", "PI"))

    if i_pn_f=="y": pn_f.close()            
                
#Selecting data according the chromosome and position values
data = list ()
Head = list()
hd_n = list()
if "vcf.gz" in i_vcfgz: file = gzip.open (i_vcfgz, "rt")
else: file = open (i_vcfgz, "rt")
for line in file:
        if line.startswith("#CHROM"): 
            
            if i_phn =="y":
                ln=line.split()
                hd=ln[:9]
                
               
                for b in lst_PI:
                    m=9
                    for a in ln[9:]:
                        if b in a and m not in hd_n:
                            hd.append(a)
                            hd_n.append(m)
                        m=m+1
                Head.append(hd)
                
            else: Head.append(line.split())
                    
            break
for line in file:
        line=line.split()
        if int(pos[0])<10:
            if line[0]==pos[0] or line[0]=="0"+pos[0] or line[0]=="Chr0"+pos[0] or line[0]=="Gm0"+pos[0] or line[0]=="CH0"+pos[0] or line[0]=="Ch0"+pos[0]: 
                #print (line)
                try:
                    if int(line[1])>= pos[1] and int(line[1])<=pos[2]:
                        
                        if i_phn =="y":
                            ln=line[:9]
                            for a in hd_n: ln.append(line[a])
                            data.append(ln)
                                                                      
                        else:
                            data.append(line)
                            
                        #print (line)
                except: continue
        else:
            if line[0]==pos[0] or line[0]=="Chr"+pos[0] or line[0]=="Gm"+pos[0] or line[0]=="CH"+pos[0] or line[0]=="Ch"+pos[0]: 
                try:
                    if int(line[1])>=pos[1] and int(line[1])<=pos[2]: 
                        if i_phn =="y":
                            ln=line[:9]
                            for a in hd_n: ln.append(line[a])
                            data.append(ln)
                                               
                        else:
                            data.append(line)
                       
                except: continue
file.close()
#(Head)
#print (data)
df = pd.DataFrame(data,columns=Head[0])
dt = list()
df[Head[0][1]]=pd.to_numeric(df[Head[0][1]])
pos = np.array(df[Head[0][1]])
upos = np.unique(pos)
for a in upos:
    ndf = df[df[Head[0][1]]==a]
    if len(ndf.index)>1:
        var =list(ndf[Head[0][4]])
        v=None
        for a in var:
            if v is None: v=str(a)
            else: v = v+","+str(a)
        ln =list(ndf.iloc[0])
        ln[4]=v
    else:
        ln =list(ndf.iloc[0])
    dt.append(ln)

data = list()
for line in dt:
    ln=ln_to_hmp(line,i_indels)
    data.append(ln)
    
#Data identifying: max, nim, position values interval, no. of positions
pos_mx = None
pos_mn = None
for line in data:
    if pos_mx is None or pos_mx<int(line[1]): pos_mx=int(line[1])
    if pos_mn is None or pos_mn>int(line[1]): pos_mn=int(line[1])
print ("\n","First position in the hapmap is: ",pos_mx,"\n",
       "Last position in the dataset is: ",pos_mn, "\n",)
print ("The genomic region size is: ", pos_mx-pos_mn, "\n")
print ("The no. of variants/positions in the hapmap is: ", len(data), "\n")

#Data oraganazing
for line in Head:
    lst1 =["RS","Ref","CHR","POS","strand","assembly","center","protLSID","assayLSID","panelLSID","Qccode"]
    n=9
    lst2 = list ()
    while n < (len(line)):
        lst2.append(line[n].rstrip())
        n=n+1 
    lst_head  = [lst1+lst2]

for line in data:
    ln = line
    a = re.findall("\d+", line[0])#Chromosome number extraction
    if not re.search("\d+",line[2]): line[2] = line[1]#Missing ID replace with position
    lst1 =[line[2],line[3],str(int(a[0])),line[1],"+","NA","NA","NA","NA","NA","NA"]
    n=9
    lst2 = list ()
    while n < (len(line)):
        lst2.append(line[n])
        n=n+1 
    data[data.index(ln)] = lst1+lst2

data = lst_head + data

#Data check 1
l = len (data[0])
for line in data:
    if l != len (line):
        print ("ERROR! Missing data.")


#Result file creation
if fres == "y": res_file = open (res_name+"_"+now.strftime("%Y")[2:]+now.strftime("%m%d%H%M")+".txt", "w")
else: res_file = open (res_name+".txt","w")

for line in data:
        n=0
        while n<(len(line)-1):
                res_file.write(str(line[n]))
                res_file.write("\t")
                n=n+1
        res_file.write(line[len(line)-1])
        res_file.write("\n")
res_file.close()

#Readme file creation
readme = ("First position in the hapmap: ",pos_mx,"\n",
          "Last position in the hapmap: ",pos_mn, "\n","The genomic region interval is (bp): ", 
          pos_mx-pos_mn,"\n","The no. of variants/positions in the hapmap is: ",len(data)-1,"\n")

if fres == "y": readme_file = open (res_name+"_readme_"+now.strftime("%Y")[2:]+now.strftime("%m%d%H%M")+".txt", "w")
else: readme_file = open (res_name+"_readme.txt","w")
for a in readme: 
    readme_file.write(str(a))
    readme_file.write("\n")
readme_file.close()

    

