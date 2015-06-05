import sys
import random
import numpy
import scipy.stats
import scipy.special
from collections import deque
import time
import math
from operator import itemgetter, attrgetter


class BinarySearch:
    NOT_FOUND = False

    def binarySearch(self, list, value, left, right):
        if right < left:
            return self.NOT_FOUND

        mid = (left + right) / 2

        if right - left == 1:
            return right

        if value > list[mid]:
            return self.binarySearch(list, value, mid, right)

        elif value <= list[mid]:
            return self.binarySearch(list, value, left, mid)

    def search(self, list, value):
        left = 0
        right = len(list) - 1
        return self.binarySearch(list, value, left, right)

start = time.clock()

# New
input = str(sys.argv[1])
# print 'What is the complete name of the input file ?'
# input = raw_input('-->')
# input = input#[:-1]
inputFH = open(input, 'r')

reportFH = str(sys.argv[2])
# print 'Choose the name of the report file?'
# reportFH = raw_input('-->')
# reportFH = reportFH#[:-1]
report = open(reportFH, 'w')

output = str(sys.argv[3])
outputFH = open(output, 'w')
outputFH.write('SNP\tChromosome\tPosition\t'+ 'Max window He1'+ '\t' + 'He1 Percentile'+ '\t '+ 'Maxwindow He2'+
               '\t'+ 'He2 Percentile' + '\t'+ 'Max window Fst' + '\t'+ 'Fst Percentile' + '\n')

numberOfWindows = int(sys.argv[4])
# numberOfWindows = int(numberOfWindows)
# numberOfWindows = int(raw_input('-->'))

currentSize = str(sys.argv[5])
startSize = currentSize = int(currentSize)
# print "What size for the shortest window ?"
# startSize = currentSize = int(raw_input('-->'))
windowIncrement = 2
maxWindowSize = startSize + windowIncrement * (numberOfWindows - 1)

resampling = int(sys.argv[6])
# print "How much resampling ?"
# resampling = int(raw_input('-->'))

sensitivity = str(sys.argv[7])
# print "Choose the sensitiviy"
# sensitivity = raw_input('-->')
# sensitivity = float(sensitivity[:-1])
sens = float(sensitivity)

####

print "Info:"
print "Number of windows: " + str(numberOfWindows)
print "Start size (min win size): " + str(startSize)
print "Window increment: " + str(windowIncrement)
print "Maximum window size: " + str(maxWindowSize)
print "Resampling value: " + str(resampling)
print "Sensitivity: " + sensitivity

# Creating windows
windows = []
while len(windows) < numberOfWindows:
    windows.append(currentSize)
    currentSize += windowIncrement

# skip the header of input file and read it to memory
inputFH.readline()

# Getting the SNPs and values
allSNPs = []
for line in inputFH:
    snpData = line.split('\t')#character used to split might need to be changed according to the format of the input file
    if snpData[7] == "1.0" and snpData[12] == "1.0":
        continue
    else:
        allSNPs.append([str(snpData[0])+','+str(snpData[1])+','+str(snpData[2]), snpData[7], snpData[12], snpData[13]]) #6,13,17 for Hexpected, 7,14,17 for Hobserved
SIZE = float(len(allSNPs))

# Start Sampling
currentWindow = {}
sampled = {}
counter = 0

for item in allSNPs:

    for window in windows:
        avgSnp = [0,0,0]

        if window in currentWindow:
            currentWindow[window][0].append(item[0])
            currentWindow[window][1].append(float(item[1]))
            currentWindow[window][2].append(float(item[2]))
            currentWindow[window][3].append(float(item[3]))
        else:
            currentWindow[window] = [deque([item[0]]), deque([float(item[1])]),
                                     deque([float(item[2])]),
                                     deque([float(item[3])])]
        if len(currentWindow[window][0]) == window:
            avgSnp[0] = sum(currentWindow[window][1]) / window
            avgSnp[1] = sum(currentWindow[window][2]) / window
            mean_avgSnp_2 = sum(currentWindow[window][3]) /window
            sq_diff = 0
            for elem in currentWindow[window][3]:
                sq_diff += (elem - mean_avgSnp_2)**2 / (window -1)
            avgSnp[2] = round(sq_diff, 5)

            TargetIndex = window /2
            Target = currentWindow[window][0][TargetIndex]

            rm = currentWindow[window][0].popleft()
            rm = currentWindow[window][1].popleft()
            rm = currentWindow[window][2].popleft()
            rm = currentWindow[window][3].popleft()

            if Target in sampled:
                sampled[Target].append([window, avgSnp[0], avgSnp[1], numpy.log(avgSnp[2])])
            else:
                sampled[Target] = [[window, avgSnp[0], avgSnp[1], numpy.log(avgSnp[2])]]

    counter += 1
    if counter % 1000 == 0:
        percentDone = (counter/ SIZE) * 100
        print "Sampling Progress: "+ str(percentDone) + "%"

print 'Done in ' + str(time.clock() - start)


del currentWindow

# Start re-sampling
currentWindow = {}
max = len(allSNPs) - 1
counter = 0

resampled = {}
for i in xrange(0, resampling):

    for window in windows:
        rdm_avgSnp = [0,0,0]
        currentWindow[window] = [[], [], []]
        while len(currentWindow[window][0]) < window:
            rdm = numpy.random.randint(0, max)
            rdmSnp = allSNPs[rdm]
            currentWindow[window][0].append(float(rdmSnp[1]))
            currentWindow[window][1].append(float(rdmSnp[2]))
            currentWindow[window][2].append(float(rdmSnp[3]))

        rdm_avgSnp[0] = sum(currentWindow[window][0]) / window
        rdm_avgSnp[1] = sum(currentWindow[window][1]) / window
        mean_rdm_avgSnp_2 = sum(currentWindow[window][2]) /window
        sq_diff = 0
        for elem in currentWindow[window][2]:
            sq_diff += (elem - mean_avgSnp_2)**2 / (window -1)
        rdm_avgSnp[2] = round(sq_diff, 5)

        if window in resampled:
            resampled[window][0].append(rdm_avgSnp[0])
            resampled[window][1].append(rdm_avgSnp[1])
            resampled[window][2].append( numpy.log(rdm_avgSnp[2]))
        else:
            resampled[window] = [[rdm_avgSnp[0]], [rdm_avgSnp[1]], [ numpy.log(rdm_avgSnp[2])]]

    counter += 1
    if counter % 1000 == 0:
        percentDone = (counter / float(resampling)) *100
        print "Resampling Progress: "+ str(percentDone) + "%"

print 'Done in ' + str(time.clock() - start)


print 'Sorting and tie ranking the resampled distributions'
tied_ranks = {}
tied_ranks2 = {}
tied_ranks3 = {}
for window in resampled:

    resampled[window][0].sort()
    value = 0
    tie = 0
    for n in resampled[window][0]:  # resampled he1
        if n > value:
            value = n
            tie += 1
        if window in tied_ranks:
            tied_ranks[window].append(tie)
        else:
            tied_ranks[window] = [tie]

    resampled[window][1].sort()
    value = 0
    tie = 0
    for n in resampled[window][1]:  # resampled he2
        if n > value:
            value = n
            tie += 1
        if window in tied_ranks2:
            tied_ranks2[window].append(tie)
        else:
            tied_ranks2[window] = [tie]

    resampled[window][2].sort()
    value = 0
    tie = 0
    for n in resampled[window][2]:  # resampled Fst
        if n > value:
            value = n
            tie += 1
        if window in tied_ranks3:
            tied_ranks3[window].append(tie)
        else:
            tied_ranks3[window] = [tie]

print 'Done in ' + str(time.clock() - start)

# Here we calculate the mean and standard deviation of the resampled distribution
distribution = {}

for window in resampled:
    distribution[window] = [numpy.mean(resampled[window][0]), numpy.mean(resampled[window][1]), numpy.mean(resampled[window][2]), numpy.std(resampled[window][0]), numpy.std(resampled[window][1]), numpy.std(resampled[window][2])]


allSNPsOut = []
counter = 0
for snp in sampled:
    snpPerc = [-1000,-1000,-1000]
    SizeHe1 = 0
    SizeHe2 = 0
    SizeFst = 0

    for window in sampled[snp]:
        PercentileIndex = BinarySearch()

    PercHe1 = (window[1] - float(distribution[window[0]][0])) / float(distribution[window[0]][3])
    PercHe1 = ( 1 + scipy.special.erf(PercHe1 / numpy.sqrt(2))) / 2

    PercHe2 = (window[2] - float(distribution[window[0]][1])) / float(distribution[window[0]][4])
    PercHe2 = ( 1 + scipy.special.erf(PercHe2 / numpy.sqrt(2))) / 2

    PercFst = (window[3] - float(distribution[window[0]][2])) / float(distribution[window[0]][5])
    PercFst = ( 1 + scipy.special.erf(PercFst / numpy.sqrt(2))) / 2

    if PercHe1 > snpPerc[0]:
        snpPerc[0] = PercHe1
        SizeHe1 = window[0]
        if PercHe2 > snpPerc[1]:
            snpPerc[1] = PercHe2
            SizeHe2 = window[0]
        if PercFst > snpPerc[2]:
            snpPerc[2] = PercFst
            SizeFst = window[0]

    snpArray = snp.split(',')
    allSNPsOut.append([ snpArray[0] , snpArray[1] , int(snpArray[2]) ,[ '\t' + str(SizeHe1) + '\t'+  str(snpPerc[0]) + '\t'+ str(SizeHe2) + '\t' + str(snpPerc[1]) + '\t' + str(SizeFst) + '\t' + str(snpPerc[2]) + '\n']])
    counter += 1
    if counter % 100 == 0:
        percentDone = (counter / SIZE) *100
        print "Percentile Calculation Progress: "+ str(percentDone) + "%"

allSNPsOut.sort(key=itemgetter(1,2))
for item in allSNPsOut:
    outputFH.write(str(item[0])+'\t'+str(item[1])+'\t'+str(item[2])+item[3][0])


outputFH.flush()
outputFH.close()
print 'Done in ' + str(time.clock() - start)

# Now let's build the report of identified sweeps
print "Building the Report"
inputFH = open(output, 'r')
# inputFH = open('SmileFinderCompleteTable.csv', 'r')

report.write('Type of Selection\t' + 'SNP'+'\t'+ 'Max window He1'+ '\t' + 'He1 Percentile'+ '\t '+ 'Maxwindow He2'+
               '\t'+ 'He2 Percentile' + '\t'+ 'Max window Fst' + '\t'+ 'Fst Percentile' + '\n')

inputFH.readline()
for line in inputFH:
    line = line.strip()
    array = line.split('\t')
    if float(array[4]) > (1 - sens) and float(array[6]) > (1 - sens) and float(array[8]) < (1 - sens):
        report.write('Ancestral\t' + line + '\n')
    if float(array[4]) > (1 - sens) and float(array[6]) > (1 - sens) and float(array[8]) > (1 - sens):
        report.write('Recent in both Populations\t' + line+'\n')
    if float(array[4]) < (1 - sens) and float(array[6]) > (1 - sens) and float(array[8]) > (1 - sens):
        report.write('Recent in Population 2\t' + line + '\n')
    if float(array[4]) >(1 - sens) and float(array[6]) < (1 - sens) and float(array[8]) > (1 - sens):
        report.write('Recent in Population 1\t' + line + '\n')

inputFH.close()
report.close