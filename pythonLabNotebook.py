#!/C:/Users/HGupta/AppData/Local/Programs/Python/Python38-32/python

# Entire script is run from the following directory: C:/Users/HGupta/InoueGSTP1methylation/Seqdata

import sys, subprocess, time, os
import numpy as np, pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

beginning_time = time.time()

# TODO Finish the usage document and command line arguments/options.
"""
Ideas for options:
    Help documentation - required input files. Structure of folders/files in working directory. No argument.
    Find the run_dmm.bat in the DNAmethylp map folder and run it from command line. No argument.
    Option to select the path for your .fsa files convertd from Excel (part of xlsTofsa function). Required argument.
"""


# bashCommand function runs Bash command from Python. If desired, return the output as a Python list
def bashCommand(findCMD, showOutput=False):
    """
    Run command for Bash from within Python. If desires, return the output as a Python list.
    :param findCMD: shell command. (semicolon does not need to be escaped)
    :param showOutput: Optional parameter. Default is to just run the shell ccommand. If True, outputs a list.
    :return: list that contains output from shell command. Each element corresponds to a line from the bash output.
    """
    out = subprocess.Popen(findCMD, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Get standard out and error
    (stdout, stderr) = out.communicate()

    # Save found files to list
    if showOutput == True:
        fileList = stdout.decode().splitlines()
        print("\nThere are %d file(s) identified from the given Bash command:\n%s" % (len(fileList), findCMD))
        return fileList


start_time = time.time()

# xls_files store the paths to all the .xls files that are located in the Sequenced_data folder.
xls_files = bashCommand('find Sequenced_data -type f -name "*.xls"', showOutput=True)
with open("xls_files.txt", "w+") as file:
    file.writelines("%s\n" % xls for xls in xls_files)


# Open .xls files given a list of their paths, read the sheet called "GSTP1", and save these sheets as .fsa files.
def xlsTofsa(xls_files, path="fsa_files/"):
    # Remove all the files that are stored in the path directory
    bashCommand('rm -r %s*' % path)

    for file in xls_files:
        file = file.rstrip()  # Remove white spaces and carriage returns at end of string
        df = pd.read_excel(file, sheet_name='GSTP1', header=None,
                           usecols="A")  # Read in GSTP1 sheet from each Excel workbook

        # Convert the 1-column pandas dataframe from Excel into a 1D numpy array. Easier to work with.
        numpy_array = df.to_numpy()

        # Save this 1D numpy array as ".fsa" file (FASTA format) in the given path folder.
        # Naming for the file: ExcelFileName.fsa
        try:
            fileName = "".join([path, file.split("/")[-1].replace(".xls", ".fsa")])
            if not os.path.exists(fileName):
                np.savetxt(fileName, numpy_array, fmt="%s")
            else:
                fileName = fileName.replace(".fsa", "_2.fsa")
                np.savetxt(fileName, numpy_array, fmt="%s")
        except IOError:
            os.mkdir(path)
            xlsTofsa(xls_files, path)
        if len(df.columns) > 1:
            # Keep track of files with more than 1 column.
            print(file, "has more than 1 column of data. Please check it.")
    print("\nYour Excel files converted to FASTA format are now stored in: %s" % path)


# Convert the GSTP1 sheet in all the Excel files in Sequenced_data folder to FASTA-formatted files (.fsa)
xlsTofsa(xls_files, path="fsa_files/")

print("--- %.3f seconds. Parse and store the .xls files as .fsa files. --- " % (time.time() - start_time))

start_time = time.time()
"""
The shell script below does the following:
# 1) find: Find all the .fsa files in all of subdirectories of Seqdata and combine into single file.
# 2) sed: In the FASTA IDs, replace all instances of "-" with "_"; remove "Label=" and F6R2; replace "__" with "_"
# 3) uniq: Every time a replacement occurs, sed repeats the line, so uniq just removes the duplicated lines.
# 4) Store these FASTA-formatted lines in combinedSeqFiles.fasta
# Note that an escape character is not required before the semicolon here, but normally needed when opening Git Bash.
"""
bashCommand(
    "find . -type f -name '*.fsa' -exec cat {} ; | sed 's/Label=//p' | sed 's/-/_/pg' | sed 's/F6R2//p' | sed 's/__/_/pg' | uniq > combinedSeqFiles.fasta")

print("--- %.3f seconds. Parse .fsa files and concatenate to combinedSeqFiles.fasta ---" % (time.time() - start_time))


# Remove duplicate sequences based on the actual sequences and ouput a deduplicated FASTA file.
def sequence_cleaner_byID(fasta_file, min_length=0, por_n=100):
    """
    sequence_cleaner takes in a FASTA file and removes duplicate sequences based on the nt sequence, not FASTA ID.
    It then writes a new FASTA with the duplicate sequences removed, and assigning the first FASTA ID found.

    :param fasta_file: Required argument of a FASTA formatted file.
    :param min_length: Optional. minimum length of the sequence be specified. Default is set to 0.
    :param por_n: Optional. Maximum percentage of 'N' nucleotides in each FASTA sequence. Maximum is set to 100%.
    :return: Formal return is the name of the deduplicated file. It also writes a deduplicated FASTA file to folder.
    """

    print("\nCleaning the %s file by removing duplicates by FASATA description." % fasta_file)
    sequenceList = {}  # Create an empty dictionary that will contain sequences and IDs
    duplications = {}  # Create a dictionary that will keep track of how many records had duplicate sequences.
    # Using the Biopython fasta parse we can read our fasta_file
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence (with the .seq method) in each record and capitalize it
        sequence = str(seq_record.seq).upper()

        # Check if the current sequence is according to the user parameters
        if (
                len(sequence) >= min_length
                and (sequence.count("N") / len(sequence)) * 100 <= por_n
                and seq_record.id != seq_record.description  # Remove if there is limited information in FASTA name
        ):
            # If the sequence passed in the test for the length and 'N' percentage parameters and it isn't already in
            # the dictionary, then the sequence and its ID are added to the dictionary. Sequence as key; ID as value.
            if seq_record.description not in sequenceList:
                sequenceList[seq_record.description] = sequence  # .description rather than .id keeps full FASTA name
            else:
                # If it is already in the dictionary, sequence and last fasta ID are placed in the duplications dictionary.
                duplications[seq_record.description] = sequence

    # Manually ensure that 4 sequences are assigned as follows based on checking duplicate occurrences of FASTA ID:
    for seq_record in SeqIO.parse("manual_sequence_corrections.fasta", "fasta"):
        sequenceList[seq_record.description] = str(seq_record.seq).upper()

    # Create a file in the same directory where you ran this script named deduplicated_[your_fasta_file] to write
    # the cleaned sequence.
    deduplicated_file = "deduplicated_" + fasta_file
    with open(deduplicated_file, "w+") as output_file:
        # Read the deduplicated dictionary and write it into the directory in a FASTA format.
        for seq_description in sequenceList:
            output_file.write(">" + seq_description + "\n" + sequenceList[seq_description] + "\n")
    print("\nThere were %d sequence records in the input file: %s." % (
        len(list(SeqIO.parse(fasta_file, "fasta"))), fasta_file))
    print("\nDuplicate sequences successfully removed! Please check %s" % deduplicated_file)
    print("There were %d sequences that had one or more duplicates and were removed." % len(duplications))
    print("There should be %d records in %s." % (len(sequenceList), deduplicated_file))
    return deduplicated_file


# Run sequence_cleaner on combinedSeqFiles.fasta.
FASTA_input_file = "combinedSeqFiles.fasta"
deduplicated_file = sequence_cleaner_byID(FASTA_input_file)

# Import the InoueDataFrame from GSTP1methylation_Analysis and store as a dataframe.
# Remove all the methylation data from and comments from the dataframe
inoueDataFrame = pd.read_excel("../InoueDataFrame.xls", sheet_name="Adjusted").iloc[:, 0:16]
print(inoueDataFrame.info())


def getUniqueValue(numpySeries):
    """
    :param numpySeries: Input is a numpy Series (1 column from a dataframe).
    :return: Expect the input numpy Series to have only 1 unique value and returns that value.
    If not a numpy Series or more than 1 unique value, then exits the script with error.
    """
    if isinstance(numpySeries, pd.Series):
        uniqueValues = numpySeries.unique().shape[0]
        if uniqueValues <= 1:
            return numpySeries.unique()[0]
        else:
            sys.exit("\n Your input \n%s\nhas more than 1 unique value." % numpySeries)
    else:
        sys.exit("\nInput is not a pandas Series.")


manualTissueDxCorrections = pd.read_excel("manual_tissuedx_corrections.xlsx")

# These manual sequences have not been through the sed command, so use python to do this.
for i, j in zip(["-", ">", "__"], ["_", "", "_"]):
    manualTissueDxCorrections['FASTA'] = manualTissueDxCorrections['FASTA'].str.replace(i, j)

listFastaIDs = []
for seq_record in SeqIO.parse(deduplicated_file, "fasta"):
    FASTA = seq_record.description
    FASTA_id = seq_record.id
    words = seq_record.id.split("_")
    identifier = words[0]
    chromat_id = seq_record.description.split()[1].replace("Chromat_id=", "")

    TissueDxNum = np.nan
    # Check if the TissueDxNum exists in Inoue data frame and if so, assign that as the tissue ID.
    for tissueid in set(inoueDataFrame['TissueDxNum']):
        if seq_record.id.find(str(tissueid)) == -1:
            pass
        else:
            TissueDxNum = int(tissueid)
    # Check if the current FASTA sequence is in the manual corrections list and assign the appropriate TissueDxNum
    if FASTA in set(manualTissueDxCorrections['FASTA']):
        TissueDxNum = getUniqueValue(
            manualTissueDxCorrections.loc[manualTissueDxCorrections['FASTA'] == FASTA, 'TissueDxNum'])

    BUI = np.nan
    tissueType = np.nan
    # To assign BUI, match the TissueDxNum in the FASTA ID to the corresponding BUI from Inoue Data Frame
    if TissueDxNum is not np.nan:
        try:
            BUI = getUniqueValue(inoueDataFrame.loc[inoueDataFrame['TissueDxNum'] == int(TissueDxNum), 'BUI'])
            tissueType = getUniqueValue(
                inoueDataFrame.loc[inoueDataFrame['TissueDxNum'] == int(TissueDxNum), 'TissueType'])
        except IndexError:
            print("\nThere is an issue with getting the BUI and tissuetype for %s" % FASTA_id)
            pass

    # If there is no TissueDxNum in the FASTA ID, then check if the BUI is written instead and assign that value.
    for batchID in set(inoueDataFrame['BUI']):
        if str(seq_record.id).find(str(batchID)) == -1:
            pass
        else:
            BUI = int(batchID)

    clone_number = np.nan
    # Clone Number is assigned as the number given immediately before "M13" in the FASTA ID, separated by "_".
    clone_number = ''.join([words[words.index(element) - 1] for element in words if "M13" in element])
    try:
        clone_number = int(clone_number)
    except ValueError:
        pass

    experimentNum = np.nan
    for descriptor in seq_record.description.split():
        if "Name=" in descriptor:
            startingPosition = descriptor.upper().find("HI")
            endingPosition = descriptor.upper().find("_", startingPosition)
            experimentNum = descriptor[startingPosition:endingPosition]
            if experimentNum == "HI":
                experimentNum = np.nan
            if "Inoue" in descriptor:
                experimentNum = "Inoue"
            if "Demarzo" in descriptor:
                experimentNum = "Demarzo"
            if "hinoue2" in descriptor:
                startingPosition = 5
                endingPosition = descriptor.find("_", startingPosition)
                experimentNum = descriptor[startingPosition:endingPosition]

    unknownID = seq_record.description.split()[-1].replace("Id=", "")
    sequence = str(seq_record.seq).upper()

    # For clones that have a BUI but not TissueDxNum, we will search their FASTA id for Ca, At, Nor, PIN
    if TissueDxNum is np.nan:
        abbreviations = {'Ca': 'cancer', 'At': 'atrophy', 'Nor': 'normal', 'PIN': 'PIN', 'LNCaP': 'LNCaP', 'WBC': 'WBC'}
        for key in list(abbreviations.keys()):
            if key in seq_record.id:
                tissueType = abbreviations[key]
                try:
                    TissueDxNum = int(getUniqueValue(inoueDataFrame.loc[
                                                         (inoueDataFrame['BUI'] == BUI) &
                                                         (inoueDataFrame['TissueType'] == tissueType),
                                                         'TissueDxNum']))
                except IndexError:
                    pass

    specimenID = np.nan
    if TissueDxNum is not np.nan:
        try:
            specimenID = getUniqueValue(inoueDataFrame.loc[inoueDataFrame['TissueDxNum'] == TissueDxNum, 'specimenID'])
        except IndexError:
            print("\nThere was an error getting the specimenID for %s" % FASTA_id)
            pass

    # keepClone = np.nan
    # if experimentNum is not np.nan:
    #     try:
    #         if experimentNum in set(inoueDataFrame.loc[inoueDataFrame['TissueDxNum'] == TissueDxNum, 'ExpNum']):
    #             keepClone = "Yes"
    #     except IndexError:
    #         keepClone = "No"
    #         pass
    #
    # uniqueClone = "Yes"

    listFastaIDs.append([
        FASTA,
        identifier,
        TissueDxNum,
        BUI,
        tissueType,
        clone_number,
        specimenID,
        chromat_id,
        experimentNum,
        sequence,
        unknownID,
        # keepClone,
        # uniqueClone
    ])

sampleMapDataFrame = pd.DataFrame(listFastaIDs,
                                  columns=['FASTA', 'Identifier', 'TissueDxNum', 'BUI', 'TissueType', 'CloneNum',
                                           'specimenID', 'Chromat_id', 'ExpNum', 'Sequence', 'Unknown ID'
                                           ])

print("\n%d FASTA records have been examined to extract the following information:."
      "\nTissueID, BUI, TissueType,CloneNum, specimenID, chromat_id, expNum, sequence, plateID." % (len(sampleMapDataFrame.index)))

# Exclude if TissueDxNum has null values in their columns.
sampleMapDataFrame = sampleMapDataFrame[sampleMapDataFrame['TissueDxNum'].notnull()]
# sampleMapDataFrame = sampleMapDataFrame[~sampleMapDataFrame['TissueType'].isin(['WBC','LNCaP'])]
print("\n%d records are in the SampleMapDataFrame after removing rows w/o TissueDxNum." % (len(sampleMapDataFrame.index)))

# Add DNAMethylMap Data to the sample Map Dataframe
methylMap = pd.read_csv("./MethylMaps/2020_11_02_all_MethylMap.csv")
sampleMapDataFrame = sampleMapDataFrame.merge(methylMap, left_on='FASTA', right_on='Full clone name', how='left')
sampleMapDataFrame = sampleMapDataFrame.drop(['Ref seq', 'Sample', 'Full clone name', 'Short clone name'], axis=1)

# Remove CpG columns, labeled 1-39
for x in range(1,40):
    del sampleMapDataFrame["%s" % str(x)]

# Remove clones if BS conversion is less than 95
sampleMapDataFrame = sampleMapDataFrame[sampleMapDataFrame['%Bs-conversion'] >= 95]
print("\n%d records are in SampleMapDataFrame after removing <95%% BS conversion." % (len(sampleMapDataFrame.index)))


# # Then, keep clone if TissueDxNum and CloneNUm exist only once in the dataframe
# sampleMapDataFrame.loc[sampleMapDataFrame.duplicated(['TissueDxNum','CloneNum'], keep=False),"UniqueClone"] = "No"
# sampleMapDataFrame.loc[sampleMapDataFrame["UniqueClone"] == "Yes", "KeepClone"] = "Yes"
#
# # Keep Clone if it has ExpNum hinoue2 or hinoue222
# for expnumhinoue in ["hinoue2","hinoue222"]:
#     sampleMapDataFrame.loc[sampleMapDataFrame["ExpNum"] == expnumhinoue, "KeepClone"] = "Yes"
#
# # Discard clone if the Experiment number is blank (i.e. null)
# sampleMapDataFrame.loc[sampleMapDataFrame['ExpNum'].isnull(),"KeepClone"] = "No"

# Sort the Dataframe by some columns
sampleMapDataFrame = sampleMapDataFrame.sort_values(['TissueDxNum', 'BUI', 'CloneNum', 'ExpNum'])
print(sampleMapDataFrame.info())

# Output the sampleMapDataFrame as a tab-separated CSV file.
sampleMapDataFrame.to_csv('sampleMapDataFrame.csv', index=False, sep=',')

# Output is a 1-column pandas Series that tells you the number of occurrences for each FASTA description
duplicatesByID = sampleMapDataFrame.pivot_table(index=['FASTA'], aggfunc='size')
duplicatesBySequence = sampleMapDataFrame.pivot_table(index=['Sequence'], aggfunc='size')

# Counts number of clones that have duplicate FASTA descriptions or duplicate FASTA sequences
num_duplicatesID = len(duplicatesByID[duplicatesByID > 1])
num_duplicatesSeq = len(duplicatesBySequence[duplicatesBySequence > 1])

print("\n%d records that are in the SampleMapDataFrame."
      "\nThere are %d duplicate IDs and %d duplicate sequences."
      % (len(sampleMapDataFrame.index), num_duplicatesID, num_duplicatesSeq,))

# Write the unique TissueDxNumbers to a file
used_tissuedxnumbers = sampleMapDataFrame.loc[sampleMapDataFrame['TissueDxNum'].notnull(), 'TissueDxNum'].unique()
print("\nThere are %d TissueDxNum in our data." % len(used_tissuedxnumbers))
present_tissueIds = sampleMapDataFrame.pivot_table(index=['TissueDxNum'], aggfunc='size')
present_tissueIds.to_csv('Used_TissueDxNum.csv', sep=',', index=True)

print("\nList of output files:"
      "\n\txls_files.txt -- List of paths to all .xls files in 'Sequenced_data' directory"
      "\n\tfsa_files/ -- Folder containing individual FASTA files, which were converted from the 'GSTP1' worksheet of the .xls files."
      "\n\tdeduplicated_combinedSeqFiles.fasta -- List of cleaned FASTA records with duplicates (by FASTA name removed)"
      "\n\tsampleMapDataFrame.csv -- Sequence recorda and extracted clinical information for sequences that had a TissueDxNum and a %% BS conversion >= 95%%"
      "\t\tFASTA	Identifier	TissueDxNum	BUI	TissueType	CloneNum	specimenID	Chromat_id	ExpNum	Sequence	Unknown ID	%%Methylation	%%Bs-conversion	%%Alignment	E-value	#CpGs in ref	#CpGs in clone"
      "\n\tUsed_TissueDxNum.csv -- List of unique Tissue IDs (TissueDxNum) found.")
"""
Notes about data:
- BUI 33388 has sequence data for PIN from ExpNum HI228. However, this PIN lesion doesn't exist in Inoue's dataframe.
Decided to remove this from our final dataframe by only calling TissueDxNum.notnull() instead of BUI.notnull() too
- 29 pairs of duplicate FASTA descriptions exist. Ran these through the DNAMethylMap. Differences in 4 situations.
Thus, these will randomly be removed when using the sequence cleaner that removes sequence records by duplicated FASTA descriptions.
- 62 pairs of duplicate sequences exist. 1 of 2 duplicates in each pair lacks a experiment number, so it gets removed. 

"""
print("\n--- %.3f seconds in total ---" % (time.time() - beginning_time))
