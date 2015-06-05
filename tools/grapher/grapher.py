#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import sys
import os

plt.figure()
plt.ylim([0, 7])

position = []
He1 = []
He2 = []
Fst = []

infile = str(sys.argv[1])
inputFH = open(infile, 'r')
inputFH.readline()
Gene = str(sys.argv[2])
chromosome = str(sys.argv[3])
strand = str(sys.argv[4])
geneStart = int(sys.argv[5])
geneEnd = int(sys.argv[6])
outfile = (sys.argv[7])
size = geneEnd - geneStart
start = geneStart - size*10
end = geneEnd + size*10
if strand == "+":
    geneStart = geneStart - 100
if strand == "-":
    geneEnd = geneEnd + 100

gene = [geneStart, geneEnd]
YAxis = [5, 5]

for line in inputFH:
    line = line.rstrip()
    array = line.split('\t')
    if array[1] == chromosome and start < int(array[2]) < end:
        position.append(float(array[2]))
        if float(array[4]) != 1:
            He1.append( - numpy.log10(1 - float(array[4])))
        else:
            He1.append( - numpy.log10(0.0000001))
        if float(array[6]) != 1:
            He2.append( - numpy.log10(1 - float(array[6])))
        else:
            He2.append( - numpy.log10(0.0000001))
        if float(array[8]) != 1:
            Fst.append( - numpy.log10(1 - float(array[8])))
        else:
            Fst.append( - numpy.log10(0.0000001))
plt.axis([start, end, 0, 10])
plt.plot(position, He1, label='He Pop1')
plt.plot(position, He2, label='He Pop2')
plt.plot(position, Fst, label='S2Fst')
plt.plot(gene, YAxis, label=Gene)
plt.title(Gene)
plt.legend(loc="upper left")
plt.xlabel("Position")
plt.ylabel("-log Percentile")

# the plotter detects types by file extension
png_out = outfile + '.png'  # force it to png
plt.savefig(png_out)

# shuffle it back and clean up
data = file(png_out, 'rb').read() 
fp = open(outfile, 'wb')
fp.write(data)
fp.close()
os.remove(png_out)

