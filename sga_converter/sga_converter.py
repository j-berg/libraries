"""
AUTHOR: Jordan A. Berg, University of Utah
IMPORTANT INFORMATION:
Must have columns in csv named 'Table', 'Row', and 'Col'
Will convert format scheme from 96 well plates to 384 grouped plates
See attached image for example
"""

#Import dependencies
import sys
import seqnorm as sn
import pandas as pd
import math

#User input file name here to convert
file = str(sys.argv[1])

#read in data frame to convert
df = pd.read_csv(file, sep=',', index_col=False)

#Remove rows that don't have info in Plate && Row && Col columns in that given row
df = df.dropna(axis = 0,subset=['Plate','Row','Col'],thresh=3)

#Convert datatypes for concerned columns
df["Plate"]=df["Plate"].astype(int)
df["Row"]=df["Row"].astype(str)
df["Col"]=df["Col"].astype(int)

df['Plate_384'] = 0
df['Row_384'] = 0
df['Col_384'] = 0

table_list1 = [1,5,9,13,17,21,25,29,33,37,41,45,49,71,75]
table_list2 = [2,6,10,14,18,22,26,30,34,38,42,46,50,72,76]
table_list3 = [3,7,11,15,19,23,27,31,35,39,43,47,51,73,77]

for index, row in df.iterrows():

    """#Convert Plate number
    if isinstance(row['Plate'], str) == True:
        continue

    if math.isnan(row['Plate']):
        continue

    elif math.isnan(row['Col']):
        continue"""

    if int(row['Plate']) < 5:
        df.loc[index,'Plate_384'] = 1
    elif int(row['Plate']) > 4 and int(row['Plate']) < 9:
        df.loc[index,'Plate_384'] = 2
    elif int(row['Plate']) > 8 and int(row['Plate']) < 13:
        df.loc[index,'Plate_384'] = 3
    elif int(row['Plate']) > 12 and int(row['Plate']) <= 16:
        df.loc[index,'Plate_384'] = 4
    elif int(row['Plate']) > 16 and int(row['Plate']) <= 20:
        df.loc[index,'Plate_384'] = 5
    elif int(row['Plate']) > 20 and int(row['Plate']) <= 24:
        df.loc[index,'Plate_384'] = 6
    elif int(row['Plate']) > 24 and int(row['Plate']) <= 28:
        df.loc[index,'Plate_384'] = 7
    elif int(row['Plate']) > 28 and int(row['Plate']) <= 32:
        df.loc[index,'Plate_384'] = 8
    elif int(row['Plate']) > 32 and int(row['Plate']) <= 36:
        df.loc[index,'Plate_384'] = 9
    elif int(row['Plate']) > 36 and int(row['Plate']) <= 40:
        df.loc[index,'Plate_384'] = 10
    elif int(row['Plate']) > 40 and int(row['Plate']) <= 44:
        df.loc[index,'Plate_384'] = 11
    elif int(row['Plate']) > 44 and int(row['Plate']) <= 48:
        df.loc[index,'Plate_384'] = 12
    elif int(row['Plate']) > 48 and int(row['Plate']) <= 52:
        df.loc[index,'Plate_384'] = 13
    elif int(row['Plate']) > 69.9 and int(row['Plate']) <= 74:
        df.loc[index,'Plate_384'] = 14
    elif int(row['Plate']) > 74 and int(row['Plate']) <= 78:
        df.loc[index,'Plate_384'] = 15
    else:
        pass


    #Convert Col
    if int(row['Plate']) in table_list1 or int(row['Plate']) in table_list3:
        if int(row['Col']) == 1:
            df.loc[index,'Col_384'] = 1
        elif int(row['Col']) == 2:
            df.loc[index,'Col_384'] = 3
        elif int(row['Col']) == 3:
            df.loc[index,'Col_384'] = 5
        elif int(row['Col']) == 4:
            df.loc[index,'Col_384'] = 7
        elif int(row['Col']) == 5:
            df.loc[index,'Col_384'] = 9
        elif int(row['Col']) == 6:
            df.loc[index,'Col_384'] = 11
        elif int(row['Col']) == 7:
            df.loc[index,'Col_384'] = 13
        elif int(row['Col']) == 8:
            df.loc[index,'Col_384'] = 15
        elif int(row['Col']) == 9:
            df.loc[index,'Col_384'] = 17
        elif int(row['Col']) == 10:
            df.loc[index,'Col_384'] = 19
        elif int(row['Col']) == 11:
            df.loc[index,'Col_384'] = 21
        elif int(row['Col']) == 12:
            df.loc[index,'Col_384'] = 23
        else:
            pass
    else:
        if int(row['Col']) == 1:
            df.loc[index,'Col_384'] = 2
        elif int(row['Col']) == 2:
            df.loc[index,'Col_384'] = 4
        elif int(row['Col']) == 3:
            df.loc[index,'Col_384'] = 6
        elif int(row['Col']) == 4:
            df.loc[index,'Col_384'] = 8
        elif int(row['Col']) == 5:
            df.loc[index,'Col_384'] = 10
        elif int(row['Col']) == 6:
            df.loc[index,'Col_384'] = 12
        elif int(row['Col']) == 7:
            df.loc[index,'Col_384'] = 14
        elif int(row['Col']) == 8:
            df.loc[index,'Col_384'] = 16
        elif int(row['Col']) == 9:
            df.loc[index,'Col_384'] = 18
        elif int(row['Col']) == 10:
            df.loc[index,'Col_384'] = 20
        elif int(row['Col']) == 11:
            df.loc[index,'Col_384'] = 22
        elif int(row['Col']) == 12:
            df.loc[index,'Col_384'] = 24
        else:
            pass


    #Convert Row
    if int(row['Plate']) in table_list1 or int(row['Plate']) in table_list2:
        if str(row['Row']) == 'A':
            df.loc[index,'Row_384'] = 'A'
        elif str(row['Row']) == 'B':
            df.loc[index,'Row_384'] = 'C'
        elif str(row['Row']) == 'C':
            df.loc[index,'Row_384'] = 'E'
        elif str(row['Row']) == 'D':
            df.loc[index,'Row_384'] = 'G'
        elif str(row['Row']) == 'E':
            df.loc[index,'Row_384'] = 'I'
        elif str(row['Row']) == 'F':
            df.loc[index,'Row_384'] = 'K'
        elif str(row['Row']) == 'G':
            df.loc[index,'Row_384'] = 'M'
        elif str(row['Row']) == 'H':
            df.loc[index,'Row_384'] = 'O'
        else:
            pass
    else:
        if str(row['Row']) == 'A':
            df.loc[index,'Row_384'] = 'B'
        elif str(row['Row']) == 'B':
            df.loc[index,'Row_384'] = 'D'
        elif str(row['Row']) == 'C':
            df.loc[index,'Row_384'] = 'F'
        elif str(row['Row']) == 'D':
            df.loc[index,'Row_384'] = 'H'
        elif str(row['Row']) == 'E':
            df.loc[index,'Row_384'] = 'J'
        elif str(row['Row']) == 'F':
            df.loc[index,'Row_384'] = 'L'
        elif str(row['Row']) == 'G':
            df.loc[index,'Row_384'] = 'N'
        elif str(row['Row']) == 'H':
            df.loc[index,'Row_384'] = 'P'
        else:
            pass

df.to_csv('96to384well_converted_' + file,sep=',', header=True,index=False)
