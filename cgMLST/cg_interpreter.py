#/usr/bin/python

import sys, os.path

#input1: file prepared by shellscript.
inputfile=sys.argv[1]

#define variables
groups=[]
outputs=["IIa","IIb","IIc","IVb","spp."]
lookup={}
targetlist=[]

#read interpreting scheme file, where sequencetypes for each locus and
#serogroup are listed. first row reads locus names, beginning marked by "#"
with open("cg_interpretation_scheme.tsv", "r") as interpreter:
  for line in interpreter:

    #If line starts with #, line contains the name of loci
    if line[0]=="#":
      #split and remove new-line char of the last entry
      targets=line.split("\t")
      temp=targets[-1]
      targets[-1]=temp[:-1]
      continue

    #else: store sequence types
    chopped=line.split("\t")
    temp=chopped[-1]
    chopped[-1]=temp[:-1]
    groups.append(chopped)

#now check sample file
with open(inputfile, "r") as infile:
  for line in infile:

    #define and reset serogroup rating of this sample:
    sample=[0,0,0,0,0]

    data=line.split(",")
    temp=data[-1]
    data[-1]=temp[:-1]

    #first line read in locus names, header.
    if line[0]=="#":
      c=0
      #store locus names
      for loci in range(1,len(data)):
 
       #some loci might be listed twice, once in cgMLST, once for PCR_target.
        #however, use only first
        if lookup.has_key(data[loci]):
          continue
        
        #store indices for locus names in split rows in dict!
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
        #loop through targets, if match: add value

        #some targets might have more than one sequencetype possible --> split in more, loop through more
        more=lists[i].split(",")

	#possibility to change value added if more than one sequencetype could fit.
        if len(more)>1:
	  value=1
  	else:
  	  value=1


        for j in range(len(more)):
          #cases:
          #*: wildcard, don't get any reward or penalty
          #+ matches + (only possible if ST is "? (new)" on export)
          #complete match: "value" added as reward
          #not found == failed: "value" as reward
          #+ matches number: no reward, no penalty
          #else: mismatch! penalty of -1
          if more[j]=="*":
            break
          elif more[j]=="+" and data[targetlist[i-1]]=="+":
            #sample[sg]=sample[sg]+2
            break
          elif more[j]==data[targetlist[i-1]]:
            sample[sg]=sample[sg]+value
            break
          elif more[j]=="-" and data[targetlist[i-1]]==0:
            sample[sg]=sample[sg]+value
            break
          elif more[j]=="+" and int(data[targetlist[i-1]]>0):
            #sample[sg]=sample[sg]+1
            break
    	else:
          sample[sg]=sample[sg]-1
      
      sg+=1 #next serogroup

    
    #create list of touples from 2 lists: sample contains rating for each serogroup for this sample,
    #output contains serogroups as literal --> sort to return the highest ranked.
    new=zip(sample,outputs)
    new.sort(reverse=True)

    #return highest ranked serogroup. if all serogroups below zero, return "failed to classify" consider revision.
    sgroup=new[0]
    if sgroup[0] < 0:
       sys.stderr.write("failed to classify %s" % data[0])
    sys.stdout.write("%s,%s\n" % (data[0],sgroup[1]))

