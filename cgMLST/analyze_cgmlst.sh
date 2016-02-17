#!/usr/bin/bash

#usage: analyze_cgmlst.sh <file1> <file2>
#where file1 contains all information and file2 contains at the sampleIDs to be analyzed with the first row as header
#file1: comma separated raw export file from SeqSphere+ (v3.1), Line-Ending: UNIX!
#file2: comma separated file, first column needs to be sampleID

inputfile=$1
selectionfile=$2

#header contains all locus names/target names
#created for the interpreter: starts with#, then all file info
printf "#" > temp.csv
head -1 $inputfile >> temp.csv

#first row in samplefile is ignored (header), first field (semicolon separated) to another temp file.
tail -n+2 $selectionfile | cut -d";" -f1 > temp_selection.csv

#select for rows with matching sampleID, write to temp file after header.
for files in $(cat temp_selection.csv)
do
  fname=$(echo $files | cut -d";" -f1)
  grep $fname $inputfile >> temp.csv
done

#remove " (which is used by SeqSphere-Export), change "not found" to 0 and failed to -, new to +
sed 's/\;/,/g;s/"//g;s/? (not found)/0/g;s/? (failed)/-/g;s/? (new)/+/g;' temp.csv > outfile.csv

#run the interpreting python script.
python cg_interpreter.py outfile.csv
rm temp.csv temp_selection.csv outfile.csv

