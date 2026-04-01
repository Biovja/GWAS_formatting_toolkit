# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 13:45:39 2020

@author: Host
"""

import re
import matplotlib.pyplot as plt
from datetime import datetime


#Input 1 - processing of the grin_results file into a list creation of the Head list
i0 = input ('Enter GRIN input file in .txt format: ')
try: file_data = open (i0,"rt", encoding="utf-8")
except:
    i0 = input ('Error. Enter GRIN input file in .txt format: ')
    file_data = open (i0,"rt",encoding="utf-8")
lst_data = list ()

for line in file_data:
    if re.search("[0-9]-[0-9]", line): line=line.replace("-","_")
    a = line.split("\t")
    l = len (a)
    n=0
    while n<l:
        if a[n] == "": a[n] = "NA"
        elif a[n] == "\n": del a[n]
        n=n+1
    l = len (a)
    if re.search("https|Link|link",a[l-1]): del a[l-1]
    lst_data.append (a)
    
lst_data_head =lst_data[0]
del lst_data[0]

#lst_data transformation
n = 0
for line in lst_data:
    if re.search ("PI[0-9]{5}$|PI[0-9]{5}-|PI[0-9]{5}[A-Z]|PI[0-9]{5}\s",line[0]):
        line[0] = line[0].replace("PI","PI0")
        lst_data[n] = line
    elif re.search ("FC[0-9]{5}$|FC[0-9]{5}-|FC[0-9]{5}[A-Z]|FC[0-9]{5}\s",line[0]):
        line[0] = line[0].replace("FC","FC0")
        lst_data[n] = line 
    n = n+1


#Choosing of the trait of interest, adding of the trait into a dictionary
lst_trait = lst_data_head[1:]
print ("\n","Choose a trait of interest for formatting: ", "\n")
m = 1
for a in lst_trait:
    print ("If you want to choose trait: ",a, "print number: ",m)
    m=m+1

i = input ("Type a number to choose a trait of interest for processing: ", )
try: n = int (i)
except:
    m = 1
    for a in lst_trait:
        print ("If you want to choose trait: ",a, "print number: ",m)
        m=m+1
    i = input ("Type a number to choose a trait of interest for processing: ", )
    n = int (i)

lst_data_trait = list() #Choosing of the trait - adding of the trait into a list
for ls in lst_data: lst_data_trait.append([ls[0],ls[n]])

ref = "NA"
for line in lst_data_trait: 
    if "PI518671" in line: ref = line[1]
    
print ("The ref phenotype from Grin file is: ",ref)

trait = lst_data_head[n]
v="str"
for line in lst_data_trait:
    if re.search("[0-9]", line[1]): v="num"
    
#Choosing of result file name
res = input("Please, enter a file name to save your numerical phenotype (do not use file name extension!). The results file will be saved as txt file and readme.txt file with details of the run: ")
fres = input ("Do you want to incorporate date and time in formate YYMMDDHHMM (e.g. 16:45 8.10.2020 - 2010081645) into the file result name? (y/n)", )
now  = datetime.now()

# Info for a readme file 
rme1 = ("The ",i0," Grin file was open for sorting.","\n"
        "\t", "Total no. of lines in in the GRIN file: ", "\n", "\t","\t", str(len(lst_data)), "\n"
        "\t", "Following traits were detected in the GRIN file: ", "\n","\t","\t", str(lst_trait), "\n",
        "\t", "\t", "The ", lst_data_head[n], " trait was choosen.", "\n",
        "\t","The ref phenotype from Grin file is: ",ref, "\n")

#Writting the ReadMe file part
if fres == "y": file_rme = open (res+"_readme_"+now.strftime("%Y")[2:]+now.strftime("%m%d%H%M")+".txt", "w")
else: file_readme = open (res+"_readme.txt", "w")

for a in rme1: file_rme.write(a)

#Processing of the trait values (sorting values for processing), adding into a list/counting max or min)
#Processing of the number trait values 
print ("Total no. of accessions in the GRIN file: ",len (lst_data_trait))

if v == "num":
    n=0
    for line in lst_data_trait:
        if ";" in line[1]:
            m = line[1].split(";")
            sm = 0
            for a in m: #Average for more then one value in line
                sm = sm+float(a)
            lst_data_trait[n]=[line[0],round(sm/(len (m)),2)]
        n=n+1    
    
    mx_data = None
    mn_data = None
    n_NA = 0
    for line in lst_data_trait:#Counting of the minimum and maximum values in grin file
        if line[1] == "NA": n_NA = n_NA+1
        else:
                if mx_data is None or float(line[1]) > float(mx_data): mx_data = line[1]
                if mn_data is None or float(line[1]) < float(mn_data): mn_data = line[1]
    print ("The GRIN file contains ", n_NA," lines with NA phenotype and ", (len(lst_data_trait)-n_NA)," lines with detected phenotype.")
    print ("The minimum value of the trait is ",mn_data, "and the maximum value is ",mx_data)
        
    rme2 = ("\t","The GRIN file contains ",str(n_NA)," lines with NA phenotype and ", str(len(lst_data_trait)-n_NA)," lines with detected phenotype.","\n",
            "\t", "The minimum value of the trait is ",str(mn_data), " and the maximum value is ",str(mx_data),".", "\n", "\n" )
        
    import csv #Writing the processed grin data into a file
    with open ("grin_data-avr.csv", "w", newline="") as file_csv:
            file_csv = csv.writer(file_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            file_csv.writerows(lst_data_trait)
     
    #Grin result trait value histogram
    delt = float(mx_data)-float(mn_data)
    lst_data_hist = list()    
    for line in lst_data_trait: 
        if line[1]!="NA":lst_data_hist.append(float(line[1]))
    if delt <3:plt.hist(lst_data_hist, bins = 20, color='royalblue')
    elif delt <5:plt.hist(lst_data_hist, bins = 30, color='royalblue')
    elif delt <7:plt.hist(lst_data_hist, bins = 50, color='royalblue')
    elif delt <9:plt.hist(lst_data_hist, bins = 70, color='royalblue')
    else:plt.hist(lst_data_hist, bins = 100, color='royalblue')
    plt.xlabel(trait)
    plt.ylabel('Count')
    if fres == "y": plt.savefig (trait +"_hist_grin_"+now.strftime("%Y")[2:]+now.strftime("%m%d%H%M")+".png")
    else: plt.savefig (trait +"_hist_grin.png") 
    plt.show() 
    ih1 = input ("Look at the panel 'Plots'. You should see a histogram of values for the chosen trait from Grin input file. When you are ready to continue, press enter.",)

else:
    d_data_trait = dict()
    n=0
    for line in lst_data_trait:
        if ";" in line[1]:
            a = line[1]
            lst_data_trait[n] = [line[0],a[:a.find(";")]]
        n=n+1
        
    for line in lst_data_trait:    
        if line[1] not in d_data_trait: d_data_trait[line[1]]=1
        else: d_data_trait[line[1]]=d_data_trait[line[1]]+1
    print ("The GRIN file contain ", d_data_trait["NA"]," lines with NA phenotype and ", (len(lst_data_trait)-(d_data_trait["NA"]))," lines with founded phenotype.")
    print ("Following phenotypes/counts were detected in the GRIN file: ", "\n", d_data_trait)
        
    rme2 = ("\t","The GRIN file contain ", str(d_data_trait["NA"])," lines with NA phenotype and ", str(len(lst_data_trait)-(d_data_trait["NA"]))," lines with founded phenotype.","\n",
        "\t","Following phenotypes/counts were detected in the GRIN file: ", "\n", "\t", "\t", str(d_data_trait), "\n", "\n")
    


for a in rme2: file_rme.write(a) #Readme file part 2
    
#input 2 - processing of the phenotype color template for sorting
file_sort = open ("Output_vzor.txt")
lst_sort = list ()
for line in file_sort: lst_sort.append (line.split())

#creating the Head list
lst_sort_Head = lst_sort [0]
del (lst_sort [0])
#print (lst_sort)

#Compare lists 1
n=0
for line in lst_sort:
            for lna in lst_data_trait:
                if lna[0] in line[0]: 
                    try: line[1]= lna[1].strip()
                    except: line[1]= lna[1]
                    lst_sort[n] = line
                
            n=n+1

#lst_data transformation 2
n = 0
for line in lst_data_trait:
    if re.search ("PI0",line[0]): line[0] = line[0].replace("PI0","PI_")
    elif re.search ("PI[0-9]",line[0]): line[0] = line[0].replace("PI","PI_")
    elif re.search ("FC0",line[0]): line[0] = line[0].replace("FC0","FC_")
    elif re.search ("FC[0-9]",line[0]): line[0] = line[0].replace("FC","FC_")
    lst_data[n] = line 
    n = n+1

#Compare lists 2
n=0
for line in lst_sort:
            for lna in lst_data_trait:
                if line[0] == lna[0]: 
                    try: line[1]= lna[1].strip()
                    except: line[1]= lna[1]
                    lst_sort[n] = line
                                           
            n=n+1


#Processing of the trait values in merged dataset 
print ("Total no. of accessions in the merged dataset: ",len (lst_sort))

if v =="num":
    mx_sort = None
    mn_sort = None
    tn_NA = 0
    for line in lst_sort:#Founding of the minimum and the maximum value in merged dataset
            if line[1] == "NA" : tn_NA = tn_NA+1
            else:
                if mx_sort is None or float(mx_sort) < float(line[1]): mx_sort = line[1]
                if mn_sort is None or float(mn_sort) > float(line[1]): mn_sort = line[1]
    
    print ("The minimum value of the trait is ",mn_sort, " and the maximum value is ",mx_sort,"\n")   
    print ("Number of accessions with NA value: ",tn_NA,"and number of accessions with detected trait value:", (len(lst_sort)-tn_NA))
            
    rme4 = ("Total no. of lines in G2G_phenotype_file_template775.txt: ", "\n", "\t", str(len(lst_sort)),"\n", "\n",
            "Total no. of accessions in the merged dataset: ","\n", "\t",str(len (lst_sort)),"\n",
            "The minimum value of the trait is ",str(mn_sort), " and the maximum value is ",str (mx_sort),".", "\n",
            "Number of accessions with NA value: ", str(tn_NA),"\n")
    for a in rme4: file_rme.write(a)
        
    #Merged dataset trait value histogram
    lst_sort_hist = list()    
    for line in lst_sort: 
        if line[1]!="NA": lst_sort_hist.append(float(line[1]))
    delt = float(mx_sort)-float(mn_sort)
    if delt<3: plt.hist(lst_sort_hist, bins = 10, color='royalblue')
    elif delt<5: plt.hist(lst_sort_hist, bins = 25, color='royalblue')
    elif delt<7: plt.hist(lst_sort_hist, bins = 35, color='royalblue')
    elif delt<9: plt.hist(lst_sort_hist, bins = 50, color='royalblue')
    else: plt.hist(lst_sort_hist, bins = 70, color='royalblue')
    plt.xlabel(trait)
    plt.ylabel('Count')
    if fres == "y": plt.savefig (trait +"_hist_"+now.strftime("%Y")[2:]+now.strftime("%m%d%H%M")+".png")
    else: plt.savefig (trait +"_hist.png")
    plt.show() 
    ih1 = input ("Look at the panel Plots. You should see a histogram for the trait value from merged dataset. When you are ready to continue, press enter.",)
        
    #Creating of lists for fractions
    f1 = list ()#WT
    f2 = list ()#MUT 1
    f3 = list ()#MUT 2
    f4 = list ()#MUT 3
    f5 = list ()#MUT 4
    f6 = list ()#MUT 5
    f7 = list ()#MUT 6
    f8 = list ()#MUT 7
    lst_fPI = [f1,f2,f3,f4,f5,f6,f7,f8]
    
    i1 = input ("Do you want to assign numeric phenotypes for your traits value ranges? (y/n)", )
    
    if i1 == "y":   
        #Select fractions for WT and MUT phenotypes 
        lst_fname = ["WT","MUT 1","MUT 2","MUT 3", "MUT 4","MUT 5","MUT 6","MUT 7"]
        lst_h = list ()
        lst_hNA = ["NA","NA","NA","NA","NA","NA","NA","NA"]
        lst_range = list ()
        print ("Assigning of the traits values range for fractions creation and values of the fractions. Use only digits! ", 
                  "For assigning of the traits values for fraction creation use numbers (dot as decimal separator). "
                  "Write smaller number first. Separate the number by a dash. If you are finished, write exit instead of trait value range. " 
                  "For assignin of fraction values use number. Recommended for 1 WT and 2,3,4... for MUT. Use NA for fractions you are not interested in.")
        n=1
        m=0
        for name in lst_fname:
                ti1 = ("Assign range of trait values for ",name," (use format 10.3-15.2): ")
                i1 = input (ti1,)
                if i1 == "Exit" or i1 == "exit" or i1 == "EXIT":
                    break
                lst_range.append(i1.split("-"))
                rn = i1.split("-")
                ti2 = ("Assign value for ",name," (recommended ",n,"): ")
                i2 = input (ti2,)
                lst_h.append(i2)
                n=n+1
                for line in lst_sort:
                    if line[1]!="NA" and float(line[1])>=float(rn[0]) and float(line[1])<=float(rn[1]):lst_fPI[m].append(line[0])
                m=m+1
        lst_h = lst_h+lst_hNA
        for line in lst_sort:#Assigning of the values for fractions
                line[1]="NA"
            
        z = -1
        for lsts in lst_fPI:
                z=z+1
                for PI in lsts:
                    for line in lst_sort:
                        if line[0]==PI:line[1]=lst_h[z]
             
        n=0
        for a in lst_range:
                print ("Fraction ",lst_fname[n],": Fraction trait value range: ",a[0],"-",a[1]," Fraction value: ",lst_h[n])
                 
                rme = ("Fraction ",lst_fname[n],": Fraction trait value range: ",a[0]," - ",a[1]," Fraction value: ",lst_h[n], "\n")
                for b in rme: file_rme.write(str(b))
                n=n+1

#Count samples from merged files
else:
    d_psort = dict ()
    for line in lst_sort:
        n=0
        if line[1] in d_psort: d_psort[line[1]]=d_psort[line[1]]+1
        else: d_psort[line[1]]=1
    print  ("The merged file contain ", d_psort["NA"]," lines with NA phenotype and ", (len(lst_sort)-d_psort["NA"])," lines with detected phenotype.")
    print ("Following phenotypes/counts are available for merged lines data set:", "\n", d_psort)
    
    rme3 = ("Total no. of lines in G2G_phenotype_file_template775.txt: ", "\n", "\t", str(len(lst_sort)),"\n", "\n",
        "Total no. of accessions in the merged dataset: ","\n", "\t",str(len (lst_sort)), "\n",
            "\t","The merged file contain ",str(d_psort["NA"])," lines with NA phenotype and ", str(len(lst_sort)-d_psort["NA"])," lines with detected phenotype.","\n",
            "Following phenotypes/counts are available for merged lines complete dataset: ", "\n","\t", str (d_psort),"\n")
    for a in rme3:
        file_rme.write (a)
    
    #Select phenotypes, assign values
    print ("\n","For assigning values use only digits!")
    print ("\n","Assign value for phenotypes, 1 is recommended for WT and 2,3,4... for MUT.", "\n", "Use exit to end the assigning. Use leters NA for phenotypes you are not interesting in.")
    d_val = dict ()
    
    lst_psort = list () #Creating of phenotype list
    for  line in lst_sort:
        if line[1] in lst_psort: continue
        else: lst_psort.append (line[1])
    
    for ls in lst_psort:
        if ls == "NA": continue
        else:
            h = input (("Assign value for ",ls," (1 is recommended for WT and 2,3,4... are recommended for MUT, for phenotypes you are not interested in use NA): "), )
            if h == "exit" or h == "EXIT":
                for line in lst_sort:
                    if re.search("[0-9]+",line[1]): continue
                    else: line[1]="NA"
                break
            else:
                if h in d_val and h!="NA":
                    print ("The value",h, "was alredy used.", d_val )
                    h = input (("Assign value for ",ls," (1 is recommended for WT and 2,3,4... are recommended for MUT, for phenotypes you are not interested in use NA): "), )
                    d_val[h]=ls
                    if h == "exit" or h == "EXIT":
                            for line in lst_sort:
                                if re.search("[0-9]+",line[1]): continue
                                else: line[1]="NA"
                    else:
                        d_val[h]=ls
                        for line in lst_sort:
                                if line[1]==ls: line[1]=h
                else:
                    d_val[h]=ls
                    for line in lst_sort:
                                if line[1]==ls: line[1]=h
    
    rme4 = ("Folowing phenotype/values pairs were selected", "\n", "\t", str (d_val), "\n",)
    for a in rme4:
        file_rme.write(a)
    
      
#Counting of the NA and WT/MUT samples in sorted dataset for SNPViz
d_result_count = {"NA":0}
for line in lst_sort:
    if line [1]=="NA":d_result_count["NA"]=d_result_count["NA"]+1
    else:
        if line [1] in d_result_count: d_result_count[line[1]]=d_result_count[line[1]]+1
        else: d_result_count[line[1]]=1

print ("\n","Total no. of accessions in final file according their values:", d_result_count) 
rme5 = ("\n","Total no. of accessions in final file according their values:","\n", "\t", str(d_result_count))
for a in rme5: file_rme.write (a)
file_rme.close()

#Data oraganazing

#Writing the results into an csv file
name_csv = res.strip()+".csv"

if fres == "y": o1 = open (res+"_"+now.strftime("%Y")[2:]+now.strftime("%m%d%H%M")+".txt", "w")
else: o1 = open (res+".txt", "w")
o1.write("<Trait>")
o1.write("\t")
o1.write("Phenotype-")
o1.write(trait)
o1.write("\n")
for a in lst_sort:
    o1.write(a[0])
    o1.write("\t")
    o1.write(str(a[1]))
    o1.write("\n")
o1.close()

        



