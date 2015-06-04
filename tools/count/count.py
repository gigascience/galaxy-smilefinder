"""
count.py
Prepare data from HGDP for Smilefinder
Wilfried Guiblet
"""

# Prepare the data from HGDP to Smilefinder
# Example of use : python count.py HGDP_FinalReport_Forward.txt hgdp-popfile.txt HGDP_Map.txt > output.csv

import sys
import optparse
import math
import os

# Get path to the galaxy count tool folder
parser = optparse.OptionParser()
parser.add_option("", "--tool_dir", dest="tool_dir")
opts, args = parser.parse_args()

# position map
SnpPositionMap = {}
positionFH = open(opts.tool_dir + '/HGDP_Map.txt', 'r')
for line in positionFH:
    (snp, chromo, pos) = line.split('\t')
    pos = pos.strip()
    SnpPositionMap[snp] = (chromo, pos)
positionFH.close()

# extract desired indexes
genoFH = open(opts.tool_dir + '/HGDP_FinalReport_Forward.txt', 'r')
line = genoFH.readline()
individuals = line.split('\t')[1:]
individuals = [x.strip() for x in individuals]

# populations map
populationsFH = open(sys.argv[1], 'r')
populationsMap = {}
for line in populationsFH:
    (population, individual) = line.split('\t')
    population = population.strip()
    individual = individual.strip()
    if population in populationsMap:
        populationsMap[population].append(individual)
    else:
        populationsMap[population] = [individual]
populationsFH.close()

# Output
outfile = open(sys.argv[2], 'w')

####

# # position map
# SnpPositionMap = {}
# positionFH = open(sys.argv[3], 'r')
# for line in positionFH:
#     (snp, chromo, pos) = line.split('\t')
#     pos = pos.strip()
#     SnpPositionMap[snp] = (chromo, pos)
# positionFH.close()

# # populations map
# populationsFH = open(sys.argv[2], 'r')
# populationsMap = {}
# for line in populationsFH:
#     (population, individual) = line.split('\t')
#     population = population.strip()
#     individual = individual.strip()
#     if population in populationsMap:
#         populationsMap[population].append(individual)
#     else:
#         populationsMap[population] = [individual]
# populationsFH.close()

# # extract desired indexes
# genoFile = sys.argv[1]
# genoFH = open(genoFile, 'r')
# line = genoFH.readline()
# individuals = line.split('\t')[1:]
# individuals = [x.strip() for x in individuals]

####

indexMap = {}
for population in populationsMap:
    indexMap[population] = []
    for individual in populationsMap[population]:
        # extract index
        index = individuals.index(individual)
        indexMap[population].append(index)

# indexMap contains the indexes of desired individuals
# in each population. indexMap-> population:index

# count
outfile.write('Name\tChr\tPos')
for population in populationsMap:
    outfile.write('\t' + population + ' MA\t'+ ' Frequency\t' + ' n\t' + 'He\t' + 'Ho')
outfile.write('\t'+'Fst'+'\n')

for line in genoFH:
    l = line.split('\t')
    currentSnp = l[0]
    populationDataMap = {}
    for population in populationsMap:
        alleleCount = {}
        Ho = 0

        for index in indexMap[population]:
            genotype = l[index + 1]
            allele = genotype[0]
            for x in range(2):
                if allele in alleleCount:
                    alleleCount[allele] +=1
                else:
                    alleleCount[allele] = 1
                allele = genotype[1]
            if str(genotype[0]) == str(genotype[1]):
                Ho = Ho + 1

       # mission is now to find max allele
        maxAllele = ''
        maxValue = 0
        for allele in alleleCount:
            if(alleleCount[allele] > maxValue):
                maxValue = alleleCount[allele]
                maxAllele = allele
        # max is in maxAllele

        n = len(indexMap[population])
        maxFrequency = float(maxValue) / n / 2
        He = 1 - ( float(maxFrequency) * 2 * ( 1 - float(maxFrequency)))
        Ho = float(Ho) / n
        populationDataMap[population] = (maxFrequency, n, He, Ho, maxAllele)


    # add chr:position here
    outfile.write(currentSnp + '\t' + SnpPositionMap[currentSnp][0] + '\t' + SnpPositionMap[currentSnp][1])

    Alleles = []
    MAF = []
    n = []
    for population in populationDataMap:
        outfile.write('\t' + str(populationDataMap[population][4])+'\t' + str(populationDataMap[population][0]) + '\t' + str(populationDataMap[population][1]) + '\t' + str(populationDataMap[population][2]) + '\t' + str(populationDataMap[population][3]))

        #Fst calculations
        Alleles[1:1] = [str(populationDataMap[population][4])]
        MAF[1:1] = [str(populationDataMap[population][0])]
        n[1:1] = [str(populationDataMap[population][1])]

    if Alleles[0] == Alleles[1]:
        MAF1 = float(MAF[0])
        MAF2 = float(MAF[1])
    else:
        MAF1 = float(MAF[0])
        MAF2 = 1 - float(MAF[1])

    #MSP
    n1 = float(n[0])
    n2 = float(n[1])

    MSP = n1*(MAF1-(MAF1+MAF2)/2)**2 + n2*(MAF2-(MAF1+MAF2)/2)**2

    #MSG
    MSG = 1/(n1+n2-2)*((n1*MAF1*(1-MAF1))+(n2*MAF2*(1-MAF2)))

    #nc
    nc = (n1+n2)-(((n1)**2+(n2)**2)/(n1+n2))

    #Fst
    if MSP == 0:
        Fst = 0
    elif MSG == 0:
        Fst = 0
    else:
        Fst = (MSP-MSG)/(MSP+(nc-1)*MSG)

    outfile.write('\t'+str(Fst))
    outfile.write('\n')

outfile.close()

# ya!
