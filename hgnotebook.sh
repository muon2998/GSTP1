# Defined path has been: C:\Users\HGupta\InoueGSTP1methylation\Seqdata
# 1) find: Find all the FASTA files (.fsa) in all of the subdirectories of current directory & these into a single file.
# 2) sed: replace all instances of "-" with "_"; "Label=" with ""; "F6R2" with ""; "__" with "_" [these occur in the names of FASTA records].
# 3) uniq: Every time a replacement occurs, sed repeats the line, so uniq just removes the duplicated lines.
# 4) Store these 26821 lines (6415 records) in combinedSeqFiles.fasta

find . -type f -name "*.fsa" -exec cat {} \; | sed 's/Label=//p' | sed 's/-/_/pg' | sed 's/F6R2//p' | sed 's/__/_/pg' | uniq > combinedSeqFiles.fasta