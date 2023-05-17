# Set working directory
#import os
#os.chdir('C:/Users/maria/Phenotype_converter')

import re
import matplotlib.pyplot as plt

#Input 1 - Process grin_results file into a list with list
i0 = input ('Enter GRIN file in .txt format: ')
try: file_data = open (i0)
except:
    i0 = input ('Error. Enter GRIN file in .txt format: ')
    file_data = open (i0)
file_sort = open ("G2G_phenotype_file_template778.txt")
lst_data = list ()
lst_sort = list ()

for line in file_data:
    if re.search("[0-9]-[0-9]", line):
        line=line.replace("-","_")
    a = line.split("\t")
    l = len (a)
    n=0
    while n<l:
        if a[n] == "":
            a[n] = "NA"
        elif a[n] == "\n":
            del a[n]
        n=n+1
    l = len (a)
    if re.search("https|Link|link",a[l-1]):
        del a[l-1]
    lst_data.append (a)
  
lst_data_head =lst_data[0]
del lst_data[0]

# Choose trait of interest, add the trait into a dictionary
lst_trait = lst_data_head[1:]
print ("\n","Select a trait: ", "\n")
m = 1
for a in lst_trait:
    print ("If you want to choose the following trait: '",a, "' Type a number: ",m)
    m=m+1

i = input ("Choose a trait of interest for the processing:", )
try:n = int (i)
except:
    m = 1
    for a in lst_trait:
        print ("If you want to choose trait ",a, "print number: ",m)
        m=m+1
    i = input ("Type a number: ", )
    n = int (i)

d_data = dict() 
for ls in lst_data:d_data[ls[0]]=ls[n]
trait = lst_data_head[n]
# Readme file
rme1 = ("The ",i0," Grin file was opend for sorting.","\n"
        "\t", "Total no. of lines in in the GRIN file: ", "\n", "\t","\t", str(len(lst_data)), "\n"
        "\t", "Following traits were detected in the GRIN file: ", "\n","\t","\t", str(lst_trait), "\n",
        "\t", "\t", "The ",trait , " trait was choosen.", "\n")

# Process phenotype in Grin dataset and assign values to phenotypes
v="str"
for key,val in d_data.items():
    if re.search("[0-9]", val):
        v="num"
    else:
        continue

d_data_trait = dict()
if v=="str": #Process the grin result data for string trait
    for key,val in d_data.items():
        if ";" not in val:continue
        else:d_data[key]=val[:(val.find(";"))]
    for key,val in d_data.items():    
        if val not in d_data_trait:d_data_trait[val]=1
        else:d_data_trait[val]=d_data_trait[val]+1
    n_NA = 0
    n_pn = 0
    for key,val in d_data_trait.items():
        if key == "NA":n_NA=n_NA+val
        else:n_pn=n_pn+val
    print ("Total no. of accessions in the GRIN file: ",len (d_data))
    print ("Following phenotypes/counts were detected in the GRIN file: ", "\n", d_data_trait,"\n"
           "Total no. of accessions with NA phenotype is ", n_NA," and total no. with detected phenotype is ", n_pn,".")
    rme2 = ("\t","Following phenotypes/counts were detected in the GRIN file: ", "\n", "\t", "\t", str(d_data_trait),"\n",
            "\t","Total no. of lines with missing (NA) phenotype is: ", str(n_NA)," and total no.of lines with detected phenotype is: ",str(n_pn),"\n", "\n")
else:#Assign numerical value to phenotype from grin result
    for key,val in d_data.items():
        if val == "NA":continue
        elif ";" not in val:d_data[key]=float(val)
        else:#Average phenotype if more than one phenotype is vavilable for a trait
            m = val.split(";")
            sm = 0
            for n in m:sm = sm+float(n)
            d_data[key]=round(sm/(len (m)),2)
    print ("Total no. of accessions in the GRIN file: ",len (d_data))
    mx_data = None
    mn_data = None
    n_NA = 0
    for key,val in d_data.items():#Determine minimum and values of the trait
        if val == "NA":n_NA = n_NA+1
        else:
            if mx_data is None or val > mx_data:mx_data = val
            if mn_data is None or val < mn_data :mn_data = val
    print ("Total no. with NA phenotype is ", n_NA," and total no. with detected phenotype is ", (len (d_data)-n_NA),"\n",
        "The minimum value of the trait is ",mn_data, "and the maximum value is ",mx_data)
    rme2 = ("\t","Total no. of lines with missing (NA) phenotype is: ",str(n_NA)," and total no. of lines with detected phenotype is: ", str(len (d_data)-n_NA),"\n",
        "\t","The minimum value of the trait is: ",str(mn_data), " and the maximum value is: ",str(mx_data),".", "\n", "\n" )
    import csv
    lst_csv = list ()
    for key,val in d_data.items():lst_csv.append([key,val])
   
    with open ("grin_data-avr.csv", "w", newline="") as file_csv:
        file_csv = csv.writer(file_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        file_csv.writerows(lst_csv)

    #Grin result trait value histogram
    delt = mx_data-mn_data
    lst_data_hist = list()    
    for key,val in d_data.items():
        if val!="NA":lst_data_hist.append(val)
    if delt <3:plt.hist(lst_data_hist, bins = 20, color='royalblue')
    elif delt <5:plt.hist(lst_data_hist, bins = 30, color='royalblue')
    elif delt <7:plt.hist(lst_data_hist, bins = 50, color='royalblue')
    elif delt <9:plt.hist(lst_data_hist, bins = 70, color='royalblue')
    else:plt.hist(lst_data_hist, bins = 100, color='royalblue')
    plt.savefig (i0+"_"+trait+"_hist_plot.png")
    plt.show() 
    ih1 = input ("Look at the panel Plots. You should see a histogram for the trait value from Grin results file. When you are ready to continue, press enter.",)


#Get 778 lines for merged data set for Accuracy tool from GRIN file
for line in file_sort:
    if re.search ("Phenotype", line):
        lst_sort_head = line.split()
    else:
        lin = re.findall ("(.+)\sNA",line)
        lst_sort = lst_sort + lin
 
#Format 778 lines for Accuracy tool format
d_sort = dict ()
for ls in lst_sort:
       
    if re.search("PI\d+-[0-9]|FC\d+-[0-9]",ls):
        b = ls.replace("-","_")
        a = re.findall("PI\d+_[0-9]|FC\d+_[0-9]", b)
        for i in a:
            if i in d_data:d_sort[ls]=d_data[i]
            else:d_sort[ls]="NA" 
        
    elif re.search("PI\d+[A-Z](?![a-z])|FC\d+[A-Z](?![a-z])|PI\d+_[0-9]|FC\d+_[0-9]",ls): 
        a = re.findall ("(PI\d+[A-Z]|FC\d+[A-Z]|PI\d+_[0-9]|FC\d+_[0-9])", ls)
        for i in a:
            if i in d_data:d_sort[ls]=d_data[i]
            else:d_sort[ls]="NA" 
               
    elif re.search("PI\d+[A-Z](?=[a-z])|PI\d+_[A-Z](?=[a-z])|FC\d+[A-Z](?=[a-z])|FC\d+_[A-Z](?=[a-z])", ls): 
        a = re.findall ("(PI\d+|FC\d+)", ls)
        for i in a:
            if i in d_data:d_sort[ls]=d_data[i]
            else:d_sort[ls]="NA"  
    
    elif re.search("PI_\d+", ls):
        b = ls.replace("PI_","PI")
        a = re.findall ("(PI\d+_[0-9]|PI\d+[A-Z]|PI\d+)", b)
        for i in a:
            if i in d_data:d_sort[ls]=d_data[i]          
            else:d_sort[ls]="NA"
               
    elif re.search("FC_\d+", ls):
        b = ls.replace("FC_","FC")
        a = re.findall ("(FC\d+_[0-9]|FC\d+[A-Z]|FC\d+)", b)
        for i in a:
            if i in d_data:d_sort[ls]=d_data[i]          
            else:d_sort[ls]="NA"
        
    elif re.search("PI\d+|FC\d+|SNL\d+", ls): 
        a = re.findall ("(PI\d+|FC\d+|SNL\d+)", ls)
        for i in a:
            if i in d_data:d_sort[ls]=d_data[i]
            else:d_sort[ls]="NA"
        
    else:
        if ls in d_data:d_sort[ls]=d_data[i]
        else:d_sort[ls]="NA"

#Result file name
res = input("Please, enter a file name to save your numerical phenotype for Accuracy tool (do not use file name extension!). The results file will be saved as csv file and readme.txt file with details of the run: ")

#Writting of the readme file
file_readme = open (res.strip()+"_readme.txt", "w")
for a in rme1: file_readme.write (a)
for a in rme2: file_readme.write (a)

#Process the trait phenotypes in merged 778 dataset - sort values for processing
d_sort_trait = dict()
if v=="str": #Process the string phenotypes
    for key,val in d_sort.items():#Counting the no. of individual phenotype in merged data    
        if val not in d_sort_trait:d_sort_trait[val]=1
        else:d_sort_trait[val]=d_sort_trait[val]+1
    print ("Total no. of accessions in the merged dataset: ",len (d_sort))
    print ("Following phenotypes/counts were detected in the merged dataset: ", "\n", d_sort_trait)
    
    # The ref line
    try: Ref = print ("The phenotype of the Reference line (PI518671) is: ",d_data["PI518671"])
    except: Ref =  print ("The phenotype of the Reference line (PI518671) has not been found:")
    
    #Select phenotypes, assign values  
    p1=input ("Select WT phenotype by typing a phenotype key:  ",) 
    h1=input ("Assign value for WT (1 is recommended):  ",)
    p2=input ("Select MUT phenotype by typing a phenotype key:  ",) 
    h2=input ("Assign value for MUT (2 is recommended):  ",)    

    #Assign NA to all but WT and MUT phenotypes, verify total no. of lines
    listd = list ()
    for key,val in d_sort.items():
        if val==p1:d_sort[key]=h1
        elif val==p2:d_sort[key]=h2
        else:d_sort[key]="NA"
            
    print ("Total no. of lines for Accuracy tool: ", len (d_sort))
    n_NA = 0
    for key,val in d_sort.items():
        if val == "NA":n_NA=n_NA+1 #Total no. of NA for Accuracy tool
    print ("Total no. of lines with NA phenotype for Accuracy tool: ", n_NA)
    
    rme3 = ("Total no. of lines in G2G_phenotype_file_template778.txt: ", "\n", "\t", str(len(lst_sort)),"\n", "\n",
        "Total no. of accessions in the merged dataset: ","\n", "\t",str(len (d_sort)), "\n",
        "Following phenotypes/counts were detected in the merged dataset: ", "\n","\t", str (d_sort_trait),"\n"
        "\t",str(Ref),"\n",
        "\t","Selected WT phenotype: ",  str (p1), " and phenotype value: ",str (h1), "\n",
        "\t","Selected MUT phenotype: ",  str (p2), " and phenotype value: ",str (h2), "\n","\n"
        "Total no. of lines for Accuracy tool: ", str(len (d_sort)),"\n"
        "\t","Total no. of lines with NA phenotype for Accuracy tool: ", str(n_NA),"\n")
    for a in rme3: file_readme.write (a)
else:
    print ("Total no. of accessions in the merged dataset: ",len (d_sort))
    mx_sort = None
    mn_sort = None
    tn_NA = 0
    for key,val in d_sort.items():#Count maximum and minimum phenotype values in 778 merged dataset
        if val == "NA":tn_NA = tn_NA+1
        else:
            if mx_sort is None or val > mx_sort:mx_sort = val
            if mn_sort is None or val < mn_sort :mn_sort = val
    print ("The minimum value of the trait is ",mn_sort, "and the maximum value is ",mx_sort)   
    print ("Number of accessions with NA value: ", tn_NA)
    #The Ref 
    try: Ref = print ("The trait value of the Reference line (PI518671) is: ",d_data["PI518671"])
    except: Ref =  print ("The trait value of the Reference line (PI518671) has not been found:") 
    
    rme3 = ("Total no. of lines in G2G_phenotype_file_template778.txt: ", "\n", "\t", str(len(lst_sort)),"\n", "\n",
        "Total no. of accessions in the merged dataset: ","\n","\t", str(len (d_sort)), "\n",
        "\t","The minimum value of the trait is ",str(mn_sort), " and the maximum value is ",str(mx_sort),"." "\n",
        "\t",str(Ref),"\n")
    for a in rme3: file_readme.write (a)
    
    #778 merged dataset trait value histogram
    lst_sort_hist = list()    
    for key,val in d_sort.items():
        if val!="NA":lst_sort_hist.append(val)
    delt = mx_sort-mn_sort
    if delt<3: plt.hist(lst_sort_hist, bins = 10, color='royalblue')
    elif delt<5: plt.hist(lst_sort_hist, bins = 20, color='royalblue')
    elif delt<7: plt.hist(lst_sort_hist, bins = 35, color='royalblue')
    elif delt<9: plt.hist(lst_sort_hist, bins = 50, color='royalblue')
    else: plt.hist(lst_sort_hist, bins = 70, color='royalblue')
    plt.savefig ("G2G_phenotype_file_template778"+"_"+trait+"_hist_plot.png")
    plt.show() 
    ih1 = input ("Look at the panel Plots. You should see a histogram for the trait value from merged dataset. When you are ready to continue, press enter.",)
    
    #Sort trait numerical phenotypes into fractions (MUT/WT/NA)
    i1 = input("Assign phenotype value range for WT (use format: 10.3-15.8): ", )
    wt = i1.split("-")
    if float (wt[0]) < float (wt[1]) and len (wt)==2: ok = 0         
    else: 
        i1 = input("The format of the assigned values is wrong. Assign phenotype value range for WT (use format: 10.3-15.8 - lower number dash higher number, use dot as decimal separator): ", )
        wt = i1.split("-")
    h1 = input("Assign value for WT (1 is recommended): ", )
    i2 = input ("Assign phenotype value range for MUT (use format: 10.3-15.8): ", )
    mut = i2.split("-")
    if float(mut[0])<float(mut[1]) and len(mut)==2:ok = 0
    else: 
        i2 = input("The format of the assigned value is wrong. Assign phenotype value range for MUT (use format: 10.3-15.8 - lower number dash higher number, use dot as decimal separator): ", )
        mut = i2.split("-")
    h2 = input("Assign value for MUT (2 is recommended): ", )
    lst_usedkey = list ()
    for key,val in d_sort.items():
        if val == "NA":continue
        elif val>=float(wt[0]) and val<=float(wt[1]):
            d_sort[key]=h1
            lst_usedkey.append(key)
        elif val>= float(mut[0]) and val<=float(mut[1]):
            d_sort[key]=h2   
            lst_usedkey.append(key)
    for key,val in d_sort.items():
        if key not in lst_usedkey:d_sort[key]="NA"
    #print ("Number of accessions in the individual fractions: ",d_th)
    
    rme5=("\t","Number of accessions with NA value: ", str(tn_NA), "\n", 
            "\t","\t","Selected range of phenotype values for WT fraction: ",str(wt[0]),"-",str(wt[1])," and fractions value: ",str (h1), "\n", 
            "\t","\t","Selected range of phenotype values for MUT fraction: ",str(mut[0]),"-",str(mut[1])," and fractions value: ",str (h2), "\n","\n")
    for a in rme5: file_readme.write(a)

#Counting of the NA and WT/MUT samples in sorted dataset for Accuracy tool
d_result_count = {"NA":0,"WT":0,"MUT":0}
for key,val in d_sort.items():
    if val=="NA": d_result_count["NA"] = d_result_count["NA"]+1
    elif val==h1: d_result_count["WT"] = d_result_count["WT"]+1
    elif val==h2: d_result_count["MUT"] = d_result_count["MUT"]+1
print ("\n","Total no. of accessions in final file for the Accuracy tool according to their values:", d_result_count) 

rme4 = ("Total no. of accessions in final file for the Accuracy tool according to their values:","\n", "\t", str(d_result_count))
         
#Write results into csv soubor 
import csv
lst_csv = list ()
for key,val in d_sort.items():lst_csv.append([key,val])
   
with open ((res.strip()+".csv"), "w", newline="") as file_csv:
    file_csv = csv.writer(file_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    file_csv.writerow(lst_sort_head)
    file_csv.writerows(lst_csv)
    
#Write to readme file
for a in rme4: file_readme.write (a)
file_readme.close()