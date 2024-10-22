tissueIDs='20219 21263 27642 29528 30559 30570 30577 34643 35230 35589 35590 35591 35597 35598 35599 35600 35602 35603 35604 35606 35608 35610 35611 35612 35613 35614 35617 35618'
buiIDs='21970 26229 26245 26243 26247 26232 25264 26254 26222 26224 26277 28962 28963 29506 29507 31351 33964 31352 33965 31353 33388 33668'
expNums='HI142 HI80 HI196 HI158 HI220 HI197 HI218 HI191 HI192 HI226 HI185 HI189 HI216'

function parsenames () {
	
	for pattern in $@

	do
		echo '--------------'
		
		echo This is for $pattern
		
		echo These are all the files that have $pattern in the file name.
		find . -type f -name *$pattern* | sort | uniq
		
		echo -e '\n'
		echo List of file names that should have $pattern in the contents of the file
		grep -ril . --include=\*.{fsa,xls,seq} -e $pattern | sort | uniq
	
		echo -e '\n'
		echo These are FASTA IDs from .fsa files that have $pattern in their contents or in file name.
		find . -type f -name *$pattern*.fsa -exec cat {} \; | grep \> | sort | uniq >> temp2.txt
		grep -ril . --include=\*.fsa -e $pattern | sort | uniq > temp.txt
		while IFS= read -r line
		do
			cat "$line" | grep \> >> temp2.txt
		done < temp.txt
		cat temp2.txt | sort | uniq

		echo '--------------'
		
		rm temp.txt
		rm temp2.txt

	done
}

# echo "This is for tissue IDs" > findingmissingfiles_TissueID.txt
# parsenames $tissueIDs >> findfiles_byTissueID.txt

# echo "This is for BUI IDs" > findingmissingfiles_BUI.txt
# parsenames $buiIDs >> findfiles_byBUI.txt

# echo "This is for ExpNums" > findingmissingfiles_ExpNum.txt
# parsenames $expNums >> findfiles_byExpNum.txt

parsenames 8106