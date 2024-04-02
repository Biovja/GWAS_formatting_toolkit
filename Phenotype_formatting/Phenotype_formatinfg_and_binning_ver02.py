 # -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 09:59:44 2023

@author: biovja00
"""
###############################################################################
#User predefined info
## Please, complete the required information.
###############################################################################

# Original phenotype file name
file = "pheno.txt"
# Original phenotype file phenotype column number (n-1, for phenotype in 2nd column use number 1)
coln = 1
# Original phenotype file  column separator (examples: "\t",",",";")
cols = "\t"
# Original phenotype file sample name format (G - Grin name format, S - Soy1066 name format , S7 -Soy775 name format )
namform = "G"
# Original phenotype type (Qual - qualitative, Quant - quantitative)
phentyp = "Quant"
# Original phenotype, Multiple phenotype separator (if none ignore)
msep = ","
#Header in the Original phenotype file ("Y","N")
header="Y"

#Result phenotype file name
res = "pheno_res.txt"
#Result phenotype file formate ("AccuCalc","MADis","AccuTool")
form = "AccuCalc"

# You can start the script run 
###############################################################################
###############################################################################
# Script body
###############################################################################
# importt
import matplotlib.pyplot as plt
import pandas as pd


###############################################################################
# Rename sample
if namform=="G" and (form=="AccuCalc" or form=="MADis"):
    key = open("Soy1066_key.txt","rt",encoding="utf-8")
    d_key=dict()
    for line in key:
        line=line.split()
        d_key[line[2]]=line[0]
    key.close()
    #print(d)
    
if namform=="G" and form=="AccuTool":
    key = open("Soy1066_key.txt","rt",encoding="utf-8")
    d_key=dict()
    for line in key:
        line=line.split()
        d_key[line[2]]=line[1]
    key.close()
    #print(d)    

def rename_sam(name_list):
    name_list_res=list()
    for nam in name_list:
        try: nam=d_key[nam]
        except: 
            nam=nam  
            #print(f"Sample name {nam} not identified in name translation key file.")
        name_list_res.append(nam)   
    return name_list_res

###############################################################################
# Import data
file=open(file,"rt")
name_list,phen_list=list(),list()
n=0
for line in file:
    if header=="Y" and n==0: n=1
    else:
        line=line.split(cols)
        name_list.append(line[0].strip())
        
        if msep in line[coln]:
            phn=line[coln].split(msep)
            df=pd.DataFrame(phn,columns=["Phenotype"])
            if phentyp=="Quant":  
                df["Phenotype"]=pd.to_numeric(df["Phenotype"],errors='coerce')
                phn=df["Phenotype"].mean()
            elif phentyp=="Qual":
                phn=df["Phenotype"].value_counts().idxmax()
            else: 
                raise Exception("Incorrect phenotype format type.")
        
        else: phn=line[coln].strip()
        phen_list.append(phn)

if namform=="G": name_list=rename_sam(name_list)
df=pd.DataFrame(name_list,columns=["Sample"])
df["Phenotype"]=phen_list
    
###############################################################################
# Data statistic
print (f"The phenotype file contains {len(df.Phenotype)} .")
if phentyp=="Quant":
    df["Phenotype"]=pd.to_numeric(df["Phenotype"],errors='coerce')
    print (f"The minimal value of the trait is {df.Phenotype.min()} and the maximal value is {df.Phenotype.max()}.")
    print (f"The average of the trait is {df.Phenotype.mean()} and the median value is {df.Phenotype.median()}.")
    print (f"The 25% quantile (1st) of the trait is {df.Phenotype.quantile(0.25)} and the 75% quantile (3rd) of is {df.Phenotype.quantile(0.75)}.")
    delt = df.Phenotype.max()-df.Phenotype.min()
    dfx=df.dropna()
    if delt <3:plt.hist(dfx.Phenotype, bins = 20)
    elif delt <5:plt.hist(dfx.Phenotype, bins = 30)
    elif delt <7:plt.hist(dfx.Phenotype, bins = 50)
    elif delt <9:plt.hist(dfx.Phenotype, bins = 70)
    else:plt.hist(dfx.Phenotype, bins = 100)
    plt.show()
elif phentyp=="Qual":
    print("The phenotype file includes the following phenotypes/phenotype counts.")
    count = df["Phenotype"].value_counts()
    print(count)
else: 
    raise Exception("Incorrect phenotype format type.")

###############################################################################
# Transform phenotype
bin_phen_list=list()
if phentyp=="Quant":
    lim = input("Please, Set limit values for WT and MU categorization in format: WT_start,WT_stop,MUT_start,MUT_stop. Example:1,2,3,3.5 ")
    lim=[float(x) for x in lim.split(",")] 
    df["WT"]=(df.Phenotype>lim[0])&(df.Phenotype<lim[1])
    df["WT"]=df["WT"].astype(int)
    df["MUT"]=(df.Phenotype>lim[2])&(df.Phenotype<lim[3])
    df["MUT"]=df["MUT"].astype(int)
    df["MUT"]=df["MUT"].replace(1,2)
    df["Trait"]=df["WT"]+df["MUT"]
    df["Trait"]=df["Trait"].replace(0,"NA")
    df["Trait"]=df["Trait"].replace(1,0)
    df["Trait"]=df["Trait"].replace(2,1)
    df=df.drop(columns=["WT","MUT","Phenotype"])
elif phentyp=="Qual":
    catW = input("Please, Set values for WT categorization in format: WT1,WT2. Example: Br,Bl ")
    catW=[x for x in catW.split(",")]
    catM = input("Please, Set values for MUT categorization in format: MUT1,MUT2. Example: Tn,Y ")
    catM=[x for x in catM.split(",")]
    df["Trait"]=df["Phenotype"]
    for p in catW:
        df["Trait"]=df["Trait"].replace(p,0)
    for p in catM:
        df["Trait"]=df["Trait"].replace(p,1)
    df.loc[~df.Trait.isin([1,0]), 'Trait'] = 'NA'
    df=df.drop(columns=["Phenotype"])

###############################################################################
# Phormate phenotype
if form=="AccuTool":
    fileA=open("Soy775_Header.txt","rt")
    lst2 = [(a.strip()) for a in fileA]
    df2=pd.DataFrame(lst2,columns=["Sample"])
    df2["Trait"]="NA"
    df = pd.concat([df, df2]).drop_duplicates(subset="Sample", keep="first")
elif form=="MADis":
    df.loc[~df["Trait"].isin(["NA"])]
elif form=="AccuCalc":
    fileA=open("Soy1066_22_Header.txt","rt")
    for line in fileA:
        lst2=line.split()[9:]
    df2=pd.DataFrame(lst2,columns=["Sample"])
    df2["Trait"]="NA"
    df = pd.concat([df, df2]).drop_duplicates(subset="Sample", keep="first")
else:
    raise Exception("Incorrect result phenotype file formate.")
###############################################################################
# Save results
df.to_csv(res, sep="\t", index=False)




