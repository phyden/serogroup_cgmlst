#/usr/bin/python

import sys, os.path

#input1: file prepared by shellscript.
inputfile=sys.argv[1]

#define variables
groups=[]

outputs=["IIa","IIb","IIc","IVb","spp.","IVb variant"]
lookup={}
targetlist=[]

#read interpreting scheme file, where sequencetypes for each locus and
#serogroup are listed. first row reads locus names, beginning marked by "#"
with open("interpreting_scheme_pcr.csv", "r") as interpreter:
  for line in interpreter:
    if line[0]=="#":
    #If line starts with #, line contains the name of loci
      #split and remove new-line char of the last entry
      targets=line.split(",")
      temp=targets[-1]
      targets[-1]=temp[:-1]
      continue
    #else: store sequence types
    chopped=line.split(",")
    temp=chopped[-1]
    chopped[-1]=temp[:-1]
    groups.append(chopped)

#print groups

#now check sample file
with open(inputfile, "r") as infile:
  for line in infile:

    #define and reset serogroup rating of this sample:
    sample=[0,0,0,0,0,0]

    data=line.split(",")
    temp=data[-1]
    data[-1]=temp[:-1]

    #first line read in locus names, header.
    if line[0]=="#":
      c=0
      #store locus names
      for loci in range(len(data)):
        #some loci might be listed twice, once in cgMLST, once for PCR_target.
        #however, use only first
        if lookup.has_key(data[loci]):
          continue
        lookup[data[loci]]=loci
      #move dict information to list
      for loci in range(1,len(targets)):
        targetlist.append(lookup[targets[loci]])
      continue


    #groups list of list which stores sequence type informations 
    #loop through serogroups (literal in "output")
    sg=0
    for lists in groups:
      for i in range(len(lists)-1):
        #loop through targets, if match: 1
        if lists[i+1]==data[targetlist[i]]:
          sample[sg]=sample[sg]+1
        elif lists[i+1]=="+" and data[targetlist[i]]=="-":
	  sample[sg]=sample[sg]+1
        elif lists[i+1]=="+" and int(data[targetlist[i]])>0:
          sample[sg]=sample[sg]+1
      #next serogroup
      sg+=1 

    #create list of touples from 2 lists: sample contains rating for each serogroup for this sample,
    #output contains serogroups as literal --> sort to return the highest ranked.
    new=zip(sample,outputs)
    new.sort(reverse=True) 

    #return highest ranked serogroup. ideally all 5 genes match. otherwise return undefined serogroup.
    sgroup=new[0]
    if sgroup[0] == 5:
      sys.stdout.write("%s,%s\n" % (data[0],sgroup[1]))
    else:
      sys.stderr.write("%s has an undefined serogroup\n" % (data[0]))
