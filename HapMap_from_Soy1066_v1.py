# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 18:31:05 2023

@author: Jana
"""
###############################################################
#  Import
import gzip
import pandas as pd

###############################################################
#  Fill in
file_vcf=gzip.open("Soy1066_Chr11.vcf.gz","rt")
rng=[11,20000,20100] #HapMap positio definiio Chromosome, Start, Stop
name="Hmp_Soy1066_Chr11_20-20K.hmp.txt" #HapMap name


###############################################################

##Function vcf line to hhmp
def ln_to_hmp(line,i_indels="one"):
    variables = [line[3]]+line[4].split(",") 

    indel="n"
    for a in variables:
        if len(a)>1: indel="y"
            
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
        
    res=list()
  
    for a in line[9:]:
            if a[0]==a[2] and "." not in a and "*" not in a: 
                try:
                    x=d_var[int(a[0])]
                    res.append(x)
                except: 
                    print(a[0],a[:3],d_var)
                    res.append("N")
                
            else: 
                res.append("N")
    
                
   
    return res

hd0=["RS","Ref","CHR","POS","strand","assembly","center","protLSID","assayLSID","panelLSID","Qccode"]


start="n"
lst=list()
for row in file_vcf:
    if row.startswith("#CHROM"): 
           ln=row.split()
           hd=hd0+ln[9:]
           hdl=len(hd)
           start="y"
    elif start=="y":
        row=row.split()
        if int(row[1])>rng[1] and int(row[1])<rng[2]:
            ln=ln_to_hmp(row,i_indels="all")
            ln0=[row[2],row[3],rng[0],row[1],"+","NA","NA","NA","NA","NA","NA"]
            lst.append(ln0+ln)
            lnl=len(ln0+ln)
 

        
df=pd.DataFrame(lst,columns=hd)
mn=df["POS"].min()
mx=df["POS"].max()
cnt=len(df["POS"])

print(f"The first pos is {mn}.")
print(f"The last pos is {mx}.")
print(f"The pos count is {cnt}.")

df.to_csv(name,sep="\t",index=False)


