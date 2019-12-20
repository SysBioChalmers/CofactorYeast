#!/usr/bin/python
################################################################################
# extractECfiles
# Reads all EC files and organizes them into one file.
#
# Yu Chen. Last edited: 2019-12-17
################################################################################
#INPUTS:
#1) Path in which the EC files are stored (from script createECfiles.py):
input_path = '/Users/cheyu/Documents/GitHub/CofactorYeast/Data/brenda_parser_cofactors/EC_files'
#2) Path in which you wish to store the final table:
output_path = '/Users/cheyu/Documents/GitHub/CofactorYeast/Data/brenda_parser_cofactors'

################################################################################

#Main Script

#Read all EC file names:
import os
os.chdir(input_path)
dir_files = os.listdir(input_path)
dir_files.sort()

output = ''
for ec in dir_files:
    if ec[0] is 'E':
        print 'Processed file ' + ec
        fid = open(ec,'r')
        for line in fid:
            output = output + ec[0:len(ec)-4] + '\t' + line

#Write output:
os.chdir(output_path)
fid = open('cofactorsBRENDA.txt','w')
fid.write(output)
fid.close()