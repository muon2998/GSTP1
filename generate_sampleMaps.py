#!/C:/Users/HGupta/AppData/Local/Programs/Python/Python38-32/python
import numpy as np, pandas as pd

"""
From the hand-curated clones, extract sequence and clinical information into a dataframe.

Column A: Full clone name
Column B: Identifier (relevant info extracted from FASTA description)
Column C: TissueDxNum (Tissue ID)
Column D: BUI (Block ID. A single block comes from one patient.)
Column E: TissueType (Normal, Atrophy, PIN, Cancer)
Column F: CloneNum (Clone Number)
Column G: specimenID (patient identifying ID)
Column H: Chromat_id (Chromatogram ID for sequencing run)
Column I: ExpNum (Experiment number)
Column J: Sequence (DNA sequence of clone)
Column K: Unknown ID (actually, corresponds to plate ID during sequencing run)
Column L: KeepClone (Yes or No -- hand-curation of clones based on matching to original dataframe)
"""

file = "Penultimate_InoueFRAME.xlsx"
allClones = pd.read_excel(file, sheet_name="ManualAssignmentALL", usecols="A:L") # Import Excel file
allClones = allClones.loc[allClones['KeepClone'] == "Yes"] # Only keep clones in dataframe that have a Yes in column L
allClones = allClones.drop(['KeepClone'], axis=1) # Remove column L as we don't need it anymore
print(allClones.info())

"""
Import additional pathology data from Inoue's dataframe. Some columns will be used to merge with sequencing data.
Columns for matching: BUI, TissueDxNum
Different columns: CAinSlide, NextTo, Merging, Capture, Frozen, DNA_nanograms, PrimaryGrade, SecondaryGrade, GleasonSum, Comment, Vasan_Comments
"""
inoueDataFrame = pd.read_excel("../InoueDataFrame.xls", sheet_name="Adjusted", usecols="A:R") # Import Excel file
# Remove repetitive columns that we don't want to merge on.
inoueDataFrame = inoueDataFrame.drop(['specimenID', 'TissueType', 'CloneNum', 'ExpNum'], axis=1)
inoueDataFrame = inoueDataFrame.drop_duplicates() # Keep only unique rows for merging
### For this tissue (SurgPath S05-70941, SpecimenID 23333, BUI 21924, TissueDxNum 34664):
### ExpNum HI222 and HI218 are both valid experiment numbers for this tissue, but dropping duplicates removes one.
print("There are %d tissue lesions in the reference/Inoue's dataframe. " % len(inoueDataFrame))

# Merge the dataframe containing additional clinical information with the first sequencing dataframe by matching BUI and TissueDxNum
mergedData = allClones.merge(inoueDataFrame, how='left')

# Import DNA Methyl Map data for all clones and merge these data with the dataframe based on 'Full clone name'.
methylMap = pd.read_csv("./MethylMaps/2020_11_02_all_MethylMap.csv")
mergedData = mergedData.merge(methylMap, left_on='Full clone name', right_on='Full clone name', how='left')
mergedData = mergedData.drop(['Ref seq', 'Sample', 'Short clone name'], axis=1) # Remove unnecessary columns

print(mergedData.info())
mergedData.to_excel("Ultimate_merged_GSTP1Dataframe.xlsx", index=False)

"""
Generate sample map files for use in the DNA Methyl Map program (Java application)

Goal 1: For each specimen, generate a separate sample map .txt file.
Goal 2: For each tissue, generate a sample map where each line is a tissue region.
Goal 3: Generate 1 whole sample map that contains all the tissues, where each line is a patient.
Output: sample_map_specimenID.txt, which is a tab-delimited file.
    - Column 1: Exact FASTA ID -- "Full clone name"
    - Column 2: Sample name -- specimenID_TissueType
    - Column 3: Shortened clone name -- Clone_TissueDxNum_CloneNum
    - Column 4: Reference sequence to be aligned to. If you want to align to each clone to every reference, put "."
"""

num_clones = 0
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def make_sample_map(dataframe, fileName = "unnamedfiled.txt"):
    """
    :param dataframe: Input dataframe for which to use to make a sample map.
    :param sampleName: List of elements that should be in sample name
    :param fileName: Name for the sample map.
    :return: No return value, but writes a tab-separated .txt file in the output
    """
    outputDataframe = dataframe
    outputDataframe['Samplename'] = outputDataframe['TissueDxNum'].map(str)
    outputDataframe['CloneName'] = "Clone_" + outputDataframe.TissueDxNum.map(str) + "_" + outputDataframe.CloneNum.map(str)
    outputDataframe['Ref'] = "."
    outputDataframe = outputDataframe[['Full clone name', 'Samplename', 'CloneName', 'Ref']]
    outputDataframe.to_csv(fileName, sep='\t', index=False, header=False)

# for patient in set(mergedData['specimenID']):
#     directory = "./Sample_Maps/"
#     filename = directory + "sample_map_" + str(patient) + ".txt" # filename is ./Sample_Maps/sample_map_specimenID
#     # Look at the main dataframe and pull the rows that match this specimenID in the current iteration.
#     outputDataframe = mergedData[mergedData['specimenID'] == patient]
#     make_sample_map(outputDataframe, filename)
#     num_clones += file_len(filename)

print("\nThere are %d clones with all the written files" % num_clones)
num_clones = 0
for tissue in set(mergedData['TissueType']):
    directory = "./Sample_Maps/"
    filename = directory + "sample_map_" + str(tissue) + ".txt"
    outputDataframe = mergedData[mergedData['TissueType'] == tissue]
    make_sample_map(outputDataframe, filename)
    num_clones += file_len(filename)

print("\nThere are %d clones with all the written files" % num_clones)

make_sample_map(mergedData, fileName="./Sample_Maps/sample_map_curated_bySpecimenID.txt")