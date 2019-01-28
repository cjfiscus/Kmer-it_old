#!/bin/bash 

A="value" 

echo $A

# checking if variable is empty
if [ -z "$A" ] ; then 
	echo "empty"  
else
	echo "has a value"
fi

temp=$(basename $(mktemp))
echo "$temp"

sample="sample1"
outputfile="$sample""$temp"
echo "$outputfile"
