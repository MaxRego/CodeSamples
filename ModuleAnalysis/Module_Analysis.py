#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Module_Analysis.py
#  
#  Copyright 2014 max <mr255509@ohio.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#	How To Run:
#		python Module_Analysis outputDirectoryOfPipeline positiveFastaFile
#
#
#  
#  

#imports
import gzip
import math
import random
from Bio import SeqIO
import shutil
import time
import sys, os, re
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes

# calls FIMO tool on positive fasta file
def callFimo(pwmFile, fastaFile, jid):
	fimoCmd = '/home/max/meme/bin/fimo'
	handle = open(fastaFile, "rU")
	numSeqs = 0
	for record in SeqIO.parse(handle, "fasta"):
		seqId = str(record.id)
		numSeqs += 1
	pVal = 0.0001#default p-value
	fimoDirName = jid + 'Fimo_Output'
	#surrounf the pwm filw with " cause it has spaces
	pwmFile = '"'+pwmFile+'"'
	cmdList = [fimoCmd, '--oc', fimoDirName,'--verbosity', str(1), '--thresh', str(pVal), pwmFile, fastaFile]
	cmdStr = ' '.join(cmdList)
	#print 'fimoCmd:', cmdStr
	os.system(cmdStr)
	
	#dict to store for each motif list of seqs that it occurs in
	tmpDict = {}
	#a dict between the seq name and the motifs that occur in it
	seqDict = {}
	#read the fimo.txt file
	with open(fimoDirName + '/fimo.txt', 'rb') as handler:
		for line in handler:
			line = line.strip()
			if re.search(r'#', line):
				continue
			lineSplit = line.split()
			motifName = lineSplit[0]
			seqName = lineSplit[1]
			start = int(lineSplit[2])
			stop = int(lineSplit[3])
			pval = float(lineSplit[5])
			#skip reverse strand if you want
			if stop < start:
				continue
			#fill the tmp dict
			if motifName not in tmpDict:
				tmpDict[motifName] = []
			if seqName not in tmpDict[motifName]:
				tmpDict[motifName].append(seqName)
			#fill the seq dict
			if seqName not in seqDict:
				seqDict[seqName] = []
			if motifName not in seqDict[seqName]:
				seqDict[seqName].append(motifName)
			
	
	motifDict = {}
	#go thru the motifs and do a coverage
	for motifName, seqList in tmpDict.items():
		cov = len(seqList)/numSeqs
		motifDict[motifName] = cov
	
	
	#check dict
	#for motifName, cov in motifDict.items():
		#print motifName, '\t', cov
	#return dict between motif name and covergae
	return motifDict, tmpDict, seqDict, numSeqs
	
# Parse the output to the MD pipeline and get important information
def parseRamiFile(ramiFile, outputDir, motifLogos, motifLogoDir):
	
	#get Depth
	myFile = open(ramiFile, "r")
	depth = 0
	for line in myFile:
		match = re.search("#Name", line)
		if match:
			depth = depth + 1
	myFile.close()		

	#make directories
	for i in range(depth):
		i = i+1
		if not os.path.exists(outputDir+"/Depth_"+str(i)):
			os.makedirs(outputDir+"/Depth_"+ str(i))
	
	#get lines in file
	myFile = open(ramiFile, "r")	
	lineList = myFile.readlines()
	myFile.close()
	
	motifsInCov = []
	index = -1
	currDepth = 0
	for line in lineList:
		index = index + 1
		match = re.search("#Name", line)
		if match:
			currDepth = currDepth + 1
			depthFile = open(outputDir+"/Depth_"+str(currDepth)+"/"+"Depth_"+str(currDepth)+"_Stats.csv", "w")
			depthFile.write(line)
			currIndex = index + 1
			switch = 1
			motif = 1
			while switch == 1:
				if currIndex < len(lineList):
					if not "*" in lineList[currIndex]:
						currLine = lineList[currIndex]
						lineSplit = currLine.split()
						if lineSplit:
							motifName = lineSplit[0].split("_",1)[1]
							motifsInCov.append(motifName)
							motifName = "logo"+motifName+".png"
							motifName = motifName.replace("-", "_")
							shutil.copy(motifLogoDir+"/"+motifName, outputDir+"/Depth_"+str(currDepth)+"/")
							os.rename(outputDir+"/Depth_"+str(currDepth)+"/"+motifName, outputDir+"/Depth_"+str(currDepth)+"/"+ str(motif)+ "_" + motifName)
							newLine = " " + "\t" + lineSplit[1] + "\t" + str(float(lineSplit[2])*100) + "%\t" + lineSplit[3] + "\t" + str(float(lineSplit[4])*100)
							newLine = newLine + "%\t"+ lineSplit[5] + "\t"+ str(float(lineSplit[6])*100) + "%\t"+ lineSplit[7] + "\t"+ str(float(lineSplit[8])*100)+"%\n"
							depthFile.write(newLine)
							motif = motif + 1
						currIndex = currIndex + 1
					else:
						switch = 0
				else:
					return motifsInCov
					break

# generate distances of each motif occurance and generate box plots					
def generateModuleStats(motif1, motif2, outputDir, fimoDir):
	
	#print "Motifs = " + motif1 + " : " + motif2 + " IN : " + outputDir + " FIMO : " + fimoDir 
	
	motif1 = motif1[4:-4]
	motif2 = motif2[4:-4]
	motif1 = motif1.replace("_", "-")
	motif2 = motif2.replace("_", "-")
	
	motifOccurs = {}
	motifOccurs[motif1] = []
	motifOccurs[motif2] = []
	
	fimoFile = fimoDir + "/fimo.txt"
	myFile = open(fimoFile, "r")
	for line in myFile:
		match = re.search("#", line)
		if not match:
			lineSplit = line.split()
			matchMotif1 = re.search(motif1, lineSplit[0])
			matchMotif2 = re.search(motif2, lineSplit[0])
			if matchMotif1:	
				matchInfo = [lineSplit[1],lineSplit[2], lineSplit[3]]
				motifOccurs[motif1].append(matchInfo) 
			elif matchMotif2:
				matchInfo = [lineSplit[1],lineSplit[2], lineSplit[3]]
				motifOccurs[motif2].append(matchInfo) 
	
	myFile.close()
	motif1Info = motifOccurs[motif1]
	motif2Info = motifOccurs[motif2]
	
	unionSeqs = {}
	
	for value in motif1Info:
		found = 0
		for innerValue in motif2Info:
			if not found:
				if innerValue[0] == value[0]:
					unionSeqs[value[0]] = 1
					found = 1
	
	#	Seq = [ [seq, start, end], [seq,start,end]] 			
	finalMotif1Dict = {}
	finalMotif2Dict = {}
	
	for key, value in unionSeqs.items():
		finalMotif1Dict[key] = []
		finalMotif2Dict[key] = []

	for key, value in unionSeqs.items():
		for innerValue in motif1Info:
			if innerValue[0] == key:
				finalMotif1Dict[key].append(innerValue)
		for innerValue in motif2Info:
			if innerValue[0] == key:
				finalMotif2Dict[key].append(innerValue)
	
	

	minimum = 99999
	maximum = -999999
	count = 0
	total = 0
	distances = []
	for key, value in finalMotif1Dict.items():
		
		tmpMotif2Occur = finalMotif2Dict[key]
		
		for i in value:
			for j in tmpMotif2Occur:
				count += 1
				tmp = abs(int(i[1])-int(j[1]))
				distances.append(tmp)
				if tmp > maximum:
					maximum = tmp
				if tmp < minimum:
					minimum = tmp
				total = total + tmp
	
	if count:
		average = float(total)/float(count)
		infoFile = open(outputDir+"/Module_Stats.csv", "w")
		infoFile.write("Motif_Name	Sequence_Name	Start\n")
		for key, value in finalMotif1Dict.items():
			for i in value:
				infoFile.write(motif1 + "	" + key + "	" + str(i[1]) + "\n")
		for key, value in finalMotif2Dict.items():
			for i in value:
				infoFile.write(motif2 + "	" + key + "	" + str(i[1]) + "\n")	 
		infoFile.write("#	#	#" + "\n")
		infoFile.write("Number_Pairs	Minimum	Maximum	Average\n")
		infoFile.write(str(len(distances))+ "	" + str(minimum) + "	" + str(maximum) + "	" + str(average) + "\n")
		infoFile.close()
		makeBoxPlot(outputDir, distances)

#call to make the box plots	
def makeBoxPlot(outputDir, distances):

	plt = figure()
	ax = axes()
	hold(True)
	bp = boxplot(distances)
	savefig(outputDir + "/boxplot.png")
		
#fills the matrix of modules and creates folders for all output	
def getModules(outputDir, motifsInCov, allMotifs, numSeqs, motifLogoDir):
	
	fimoDir = outputDir + "/Fimo_Output"
	
	#make output directory for module analysis
	outputDir = outputDir + "/Module_Analysis"	
	if not os.path.exists(outputDir):
			os.makedirs(outputDir)
	
	#create Matrix
	moduleMatrix = [[0 for x in xrange(len(motifsInCov))] for x in xrange(len(motifsInCov))] 
			
	#Fill Matrix		
	i = -1
	j = -1		
	for targetMotif in motifsInCov:
		i = i + 1
		j = 0
		targetMotifCov = allMotifs["_" + targetMotif]
		for currMotif in motifsInCov:
			total = 0
			currMotifCov = allMotifs["_" + currMotif]
			for gene in targetMotifCov:
				if gene in currMotifCov:
					total = total + 1
			moduleMatrix[i][j] = total
			j = j + 1
		
	#write matrix to file		
	moduleFile = open(outputDir+"/Module_Stats.csv", "w")		
	header = "\t"
	for motif in motifsInCov:
		header = header + motif + "\t"
	moduleFile.write(header+"\n")
	
	i = 0
	for motif in motifsInCov:
		line = motif + "\t"
		for val in moduleMatrix[i]:
			line = line + str(val) + "\t"
		moduleFile.write(line+"\n")
		i = i + 1

	moduleFile.close()		

	
	#make motif dictionary to get coverage
	topModules = {}
	k = -1
	v = 0
	for i in motifsInCov:
		v = 0
		k = k + 1
		for j in motifsInCov:
			if i != j:
				if v < k:
					topModules[i+":"+j] = moduleMatrix[k][v]
			v = v + 1	
	
	ID = 1
	moduleFile = open(outputDir+"/Module_Coverage.csv", "w")
	moduleFile.write("ID"+"\t"+"Motif_1"+"\t"+"Motif_2"+"\t"+"NumSeqs"+"\t"+"Coverage"+"\n")
	for key, value in topModules.items():
		cov = (float(value)/float(numSeqs))*100
		keySplit = key.split(":")
		line = str(ID) + "\t" + keySplit[0] + "\t" +keySplit[1] + "\t" + str(value) + "\t" + str(cov) + "\n"
		moduleFile.write(line)
		if not os.path.exists(outputDir+"/ID_"+str(ID)):
			os.makedirs(outputDir+"/ID_"+ str(ID))
		motif1 = "logo"+keySplit[0]+".png"
		motif2 = "logo"+keySplit[1]+".png"
		motif1 = motif1.replace("-", "_")
		motif2 = motif2.replace("-", "_")
		shutil.copy(motifLogoDir+"/"+motif1, outputDir+"/ID_"+str(ID)+"/")
		shutil.copy(motifLogoDir+"/"+motif2, outputDir+"/ID_"+str(ID)+"/")

		generateModuleStats(motif1, motif2, outputDir+"/ID_"+str(ID), fimoDir)
		
		ID = ID + 1
		
						
	moduleFile.close()
	
	
	
def main():
	
	startTime = time.time()
	print "******************\n" + "Start Report Generation : Total Time : " + str(0) + " Seconds"
	
	#check arguements 
	if len(sys.argv) != 3:
		 print "Incorrect amount of arguements"
		 exit(0)
	
	#get arguements
	outputFiles = [f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1],f)) ]
	seqFile = sys.argv[2]
	pwmFile = sys.argv[1]
	ramiFile = sys.argv[1]
	
	#generate output file
	outputDir = sys.argv[1]+"/Report"
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	
	#get motif logos
	motifLogoDir = sys.argv[1] + "/" + "logos/Combination of All Motifs"
	motifLogos = [f for f in listdir(motifLogoDir) if isfile(join(motifLogoDir,f)) ]
	
	#get pwmFile
	for files in outputFiles:
			match = re.search("Combination of All Motifs__meme.txt", files)
			if match:
				pwmFile = pwmFile +"/"+ files
	
	#get ramiFile
	for files in outputFiles:
			match = re.search("rami_format2.txt", files)
			if match:
				ramiFile = ramiFile +"/"+ files
	
	#generate report, and get the motifs that are in the coverage report
	motifsInCov = parseRamiFile(ramiFile, outputDir, motifLogos, motifLogoDir)
	
	#call FIMO
	print "Starting FIMO : Toatl Time : " + str(time.time() - startTime)
	motifDict, tmpDict, seqDict, numSeqs = callFimo(pwmFile, seqFile, outputDir+"/")
	print "Finished FIMO : Total Time : " + str(time.time() - startTime)
	
	#get Modules
	getModules(outputDir, motifsInCov, tmpDict, numSeqs, motifLogoDir)
	
	print "Finished Module Analysis : Total Time : " + str(time.time() - startTime)
	print "**************************************"+"\n"
	
	return 0

if __name__ == '__main__':
	main()

