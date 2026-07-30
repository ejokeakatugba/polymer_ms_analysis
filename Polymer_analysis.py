# required package: pandas, numpy, openpyxl

import pandas as pd
from numpy import *

# import excel file with mass vs intensity list
df0=pd.read_excel('C:\MP\MP94.xlsx',skiprows=range(1,8),usecols=[0])
df1=pd.read_excel('C:\MP\MP94.xlsx',skiprows=range(1,8),usecols=[1])

# input the searching polymer repeating unit
formula = 'CH2'
# set absolute intensity threshold
bkg = 50
# set mass tolerance in amu
mass_tolerance = 0.005
################################################################################################################################################################
# exact mass of elements
C = 12
H = 1.007825
O = 15.994915
N = 14.003074
S = 31.972072
# check if m/z with the searching mass difference exist in the mass list
# if exist, append it to the new list in the corresponding dictionary
def check_mass(mass, list_index, list_mass, list_inten):
    next_mass = mass + common_difference
    for a in mz:
        if abs(next_mass-a)<=mass_tolerance:
            index=mz.index(a)
            list_index.append(index)
            list_mass.append(mz[index])
            list_inten.append(inten[index])
            mz[index]=0
            inten[index]=0
            check_mass(a, list_index, list_mass, list_inten)

def translate_index(index_list, raw_list, new_list):
        for i in index_list:
            new_list.append(raw_list[i])

# calculate exact mass of the repeating unit
def mw(formula):
    expanded = ''
    for character in formula:
        if character.isdigit():
            expanded += expanded[-1] * (int(character) - 1)
        else:
            expanded += character
    lexp = list(expanded)
    for i in range(len(lexp)):
        if lexp[i] == 'C':
            lexp[i] = C
        if lexp[i] == 'H':
            lexp[i] = H
        if lexp[i] == 'O':
            lexp[i] = O
        if lexp[i] == 'N':
                lexp[i] = N
        if lexp[i] == 'S':
            lexp[i] = S
    return sum(lexp)

# extract list from excel
mz0=list(df0[df0.columns[0]])
inten0=list(df1[df1.columns[0]])
mz=[]
inten=[]
mz1=[]
inten1=[]

# filter out all intensity<bkg
for i in range(0,len(mz0)):
    if inten0[i]>bkg:
        mz.append(mz0[i])
        inten.append(inten0[i])

for i in range(0,len(mz0)):
    if inten0[i]>bkg:
        mz1.append(mz0[i])
        inten1.append(inten0[i])

common_difference = mw(formula)

# creat dictionary for repeating unit lists
d_index = {}
d_mass = {}
d_inten = {}

print('\nsearching for repeating unit:',formula)
print('\nexact mass of the repeating unit:',common_difference)

# search and add repeating unit lists to dictionary
for mass in mz:
        index=mz.index(mass)
        d_index[mass] = [index]
        d_mass[mass] = [mass]
        intensity=inten[index]
        d_inten[mass] = [intensity]
        check_mass(mass, d_index[mass], d_mass[mass], d_inten[mass])
#print('mass dictionary:', d_mass)
# sort dictionary based on the number of continuous fragments, from long to short
sorted_index = sorted(d_index.items(), key=lambda x: len(x[1]), reverse=True)


# to check the longest fragment sequency, add the following two lines:
# max_len = len(sorted_index[0][1])
# print('max length:', max_len)
sorted_mass = []
for i in range(0,len(sorted_index)):
    new_mass = []
    translate_index (sorted_index[i][1], mz1, new_mass)
    new_mass_set = (sorted_index[i][0],new_mass)
    sorted_mass.append(new_mass_set)
# print fragments mass list, in long to short order:
print('\nfragments mass list:', sorted_mass)

sorted_inten = []
for i in range(0,len(sorted_index)):
    new_inten = []
    translate_index (sorted_index[i][1], inten1, new_inten)
    new_inten_set = (sorted_index[i][0],new_inten)
    sorted_inten.append(new_inten_set)
# print fragments intensity list, in long to short order:
# print('correspond intensity list:',sorted_inten)
# count total number of fragments, and the total intensity of all fragments
total_len = total_inten = 0
for i in range(0,len(sorted_inten)):
    if len(sorted_index[i][1])>3: # only count if more than 3 continuous fragments
        total_len = total_len + len(sorted_index[i][1])
        total_inten = total_inten + sum(sorted_inten[i][1])

print('\ntotal number of fragment peaks:', total_len)
print('\ntotal intensity of fragment peaks:', int(total_inten))
