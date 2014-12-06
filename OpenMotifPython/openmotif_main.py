#!/usr/bin/python
########################################################################################################################
#
#	Program:	openmotif_main.py
#	Author:		Max Rego
#	Email:		mr255509@ohio.edu
#
#	Description:
#
#	Date:		5-27-2013
#
#	notes:
#		if min and max word length are equal, we get errors
#
#
########################################################################################################################

#imports
import matplotlib
matplotlib.use('Agg')
import argparse, os, sys, subprocess, shutil, math, Levenshtein, time
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio import Motif
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import scipy.cluster.hierarchy
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.pyplot import show, savefig
import visualize_motifs

######## CLASSES ########################################################################################################

class Word:
	'Words that are worth storing after they are scored by wordseeker'
	wordCount = 0
	#Word,S,S/allS,E_s,O,E,O/S,O-S,O/E,O*ln(O/E),S*ln(S/E_s),Pval

	def __init__(self, word, markovOrder, S, SdivAllS, E_S, O, E, OdivS, OminusS, OdivE, OtimesNatLog, StimesNatLog, Pval):
		self.word = word
		self.order = int(markovOrder)
		self.length = len(word)
		self.S = float(S)
		self.SdivAllS = float(SdivAllS)
		self.E_S = float(E_S)
		self.O = int(O)
		self.E = float(E)
		self.OdivS = float(OdivS)
		self.OminusS = float(OminusS)
		self.OdivE = float(OdivE)
		self.OtimesNatLog = float(OtimesNatLog)
		self.StimesNatLog = float(StimesNatLog)
		self.Pval = Pval
		self.zScore = 0
		self.cluster = 0		
		self.seqDict = {}
		Word.wordCount += 1
		
   
	#accessors
	def getWord(self):
		return self.word
	def getOrder(self):
		return self.order	
	def getLength(self):
		return self.length
	def getS(self):
		return self.S
	def getSdivAllS(self):
		return self.SdivAllS
	def getE_S(self):
		return self.E_S
	def getO(self):
		return self.O
	def getE(self):
		return self.E
	def getOdivS(self):
		return self.OdivS
	def getOminusS(self):
		return self.OminusS
	def getOdivE(self):
		return self.OdivE
	def getOtimesNatLog(self):
		return self.OtimesNatLog
	def getStimesNatLog(self):
		return self.StimesNatLog
	def getPval(self):
		return self.Pval
	def getZScore(self):
		return self.zScore
	def getCluster(self):
		return self.cluster
#end Word Class

class MyMotif:
        """ MyMotif class """
        def __init__(self, seedWord):
                self.seedWord = seedWord
                self.wordList = list()	#list of simialr words to seed word
                self.index = 0
                self.size = 0 		#how many words for this motif
                self.seqDict = {} 	#the sequences that this motif hits
                self.seqCount = 0 	#number of seqeunces
		self.seqCov = 0
#end Motif Class

#########################################################################################################################

####### FUNCTIONS  ######################################################################################################

def callWordSeeker(fastaFile, jid, minWLen, maxWLen):
	"""Call wordseeker with scoring off for min - max word lengths and move the csv to an output directory.
	Parameters
    	----------
    	fastaFile : string
        	Name of the fasta file in the cwd.
    	jid : string
        	The job ID assigned from the command line to use in naming output.
	minWLen : int
        	The minimum word length to pass to wordseeker to enumerate words.
	maxWLen: int
		The maximum word length to pass to wordseeker to enumerate words. 
	Return
	----------
	outputDir : string
		The location of the directory containing the csv files outputed from wordseeker with scoring turned off.
	"""

	#setting up output directory to store output in
	os.mkdir("".join([jid,'_seekerOut_',str(minWLen),'_',str(maxWLen)]))
	cwd = os.getcwd()
	outputDir = "".join([cwd,'/',jid,'_seekerOut_',str(minWLen),'_',str(maxWLen)])	

	#primes the control for the while loop, sets up the needed variable names that dont change on each iteration
	i = minWLen
	execPath = '/bin/OWEFexec'
	fastaFile = "".join([cwd,'/',fastaFile])
	cmd = "".join([cwd,execPath])
	prefix = "".join([jid,'_seeker'])
	
	while (i <= maxWLen): 	
		#sets the new word length and sets the args to call wordseeker
		wordLength = i	
		args = [cmd, '--count', '-i', fastaFile,'-l', str(wordLength), '-n', '--prefix', prefix]
	
		#calls wordseeker and waits for it to finish
		process = subprocess.Popen(args,stdout=subprocess.PIPE)
		process.wait()

		#sets up the csv file name and the directory to delete upon finishing
		countFileName = ''.join([cwd, '/', prefix, '_', str(wordLength), '/', prefix, '_', str(wordLength), '.csv'])
        	deleteDirName = ''.join([cwd, '/', prefix, '_', str(wordLength)])	
		
		#executing move of csv and deleting directory 
		shutil.move(countFileName, outputDir)
		shutil.rmtree(deleteDirName)
		
		#iterate to the next word length
		i += 1
	return outputDir
#end callWordSeeker

def callWordSeekerScoring(fastaFile, jid, markovOrder, savedLengths):
	"""Call wordseeker with scoring on for all lengths in savedLengths and move the csv to a directory.
	Parameters
    	----------
    	fastaFile : string
        	Name of the fasta file in the cwd.
    	jid : string
        	The job ID assigned from the command line to use in naming output.
	markovOrder : int
        	The markov order to be passed to wordseeker. 
	savedLengths: dictionary
		The lengths that are significant(numWords > minWordCount): key is the length, value is the number of words that are > minSeqCoverage.  
	Return
	----------
	outputDir : string
		The location of the directory containing the csv files outputed from wordseeker with scoring on.
	"""

	#setting up output directory to store output in
	os.mkdir("".join([jid,'_seekerScoring_',str(markovOrder)]))
	cwd = os.getcwd()
	outputDir = "".join([cwd,'/',jid,'_seekerScoring_',str(markovOrder)])	

	#setting up needed variable for the args to call wordseeker	
	execPath = '/bin/OWEFexec'
	fastaFile = "".join([cwd,'/',fastaFile])
	cmd = "".join([cwd,execPath])
	prefix = "".join([jid,'_seeker'])
	

	for key, value in savedLengths.items():
		wordLength = key	
		args = [cmd, '--count', '-i', fastaFile, '-l', str(wordLength), '-n', '--score', '-o', str(markovOrder), '-p' ,'--prefix', prefix]			    		
		
		#calls wordseeker and waits for it to finish
		process = subprocess.Popen(args,stdout=subprocess.PIPE)
		process.wait()

		#sets up the csv file name and the directory to delete upon finishing
		scoreFileName = ''.join([cwd, '/', prefix, '_', str(wordLength), '_', str(markovOrder)])
		seekerScoreFiles = os.listdir(scoreFileName)
		
		for outFile in seekerScoreFiles:
			target = ''.join([prefix, '_', str(wordLength),'_'])
			if outFile.find(target) != -1:
				targetFileName = ''.join([cwd, '/', prefix, '_', str(wordLength), '_', str(markovOrder), '/', outFile])			
				#executing move of csv and deleting directory 
				shutil.move(targetFileName, outputDir)
				shutil.rmtree(scoreFileName)
	return outputDir
#end callWordSeekerScoring

def getFastaNumSeqs(fastaFile):
	"""Parses the fastaFile and finds how many sequences there are in the file.
	Parameters
    	----------
    	fastaFile : string
        	Name of the fasta file in the cwd.
	Returns
	----------
	numSeqs : int
		The number of sequences in the file.		
	"""
	numSeqs = 0
	handle = open(fastaFile, "rU")
	
	for record in SeqIO.parse(handle, "fasta"):
		seqId = str(record.id)
		#print 'SEQ:', record.seq
		numSeqs += 1
	
	handle.close()
	return numSeqs
#end getFastaNumSeqs

def getSeqDict(fastaFile):
	seqDict = {}
	i = 1
	handle = open(fastaFile, "rU")
	for record in SeqIO.parse(handle, "fasta"):
		seqId = str(record.id)
		seqName = "seq_"+ str(i)
		seqDict[seqId] = seqName
		i += 1
	return seqDict
#end getSeqDict

def getWordCount(seekerOutputDir, numSeqs, minWCount):
	"""Finds the word count for each length that are > minSeqCoverage and have words > minWCount.
	Parameters
    	----------
    	seekerOutputDir : string
        	Location of the output from wordseeker with scoring off.
    	numSeqs : float (will always be ceil(value))
        	A word must have sequences > numSeqs to be considered significant.
	minWCount:
		A length must have words found to be significant from above parameter per word length > minWCount. 
	Return
	----------
	savedLengths : Dictionary
		A dictionary with the lenths that pass the parameters above, key is the length, value is the significant word count.
	"""
	
	#returns a list with each output file in the directory
	seekerFiles = os.listdir(seekerOutputDir)

	#holds the final dictionary with word lenths that pass the minWCount	
	savedLengths = {}

	for outputFile in seekerFiles:	
		IN = open("".join([seekerOutputDir,'/',outputFile]), 'r')
		words = 0
		length = 0		
		for line in IN:
			if line.find('#') != -1:
				#do nothing its the header of the csv file
				words = 0
				length = 0	
			else:
				#lineSplit[0] = Word; lineSplit[1] = Sequences; lineSplit[2] = Count
				lineSplit = line.split(',')
				length = len(lineSplit[0])
				if float(lineSplit[1]) >= numSeqs:
					words += 1
		#if there are enough words >= sequence coverage then add it to the Dict 		
		if words >= minWCount:	
			savedLengths[length] = words
		
		IN.close()

	#for key, value in savedLengths.items():
	#    	print key, value	

	return savedLengths
#end getWordCount

def setWordClasses(scoringOutputDir, markovOrder):
	"""Takes each word and creates a class "word" object for it.
	Parameters
    	----------
    	scoringOutput : string
        	Location of the directory that holds the output from wordseeker with scoring turned on.
    	markovOrder : int
        	The markov order that was passed to wordseeker. 
	Return
	----------
	wordList : Dictionary
		A dictionary with the words that were scored by wordseeker, key is the word, value is the class "word" object.
	"""

	scoringFiles = os.listdir(scoringOutputDir)	
	wordList = {}
	i = 0

	for File in scoringFiles:
		IN = open("".join([scoringOutputDir,'/',File]), 'r')
		word = ""
		S = 0
		SdivAllS = 0
		E_S = 0
		O = 0
		E = 0
		OdivS = 0
		OminusS = 0
		OdivE = 0
		OtimesNatLog = 0
		StimesNatLog = 0
		Pval = 0		
		for line in IN:
			if line.find('#') != -1:
				#do nothing its the header of the csv file
				word = ""
				S = 0
				SdivAllS = 0
				E_S = 0
				O = 0
				E = 0
				OdivS = 0
				OminusS = 0
				OdivE = 0
				OtimesNatLog = 0
				StimesNatLog = 0
				Pval = 0	
			else:	
				#Word,S,S/allS,E_s,O,E,O/S,O-S,O/E,O*ln(O/E),S*ln(S/E_s),Pval
				#assigning all variables for this line in the csv file				
				lineSplit = line.split(',')
				word = lineSplit[0]
				S = lineSplit[1]
				SdivAllS = lineSplit[2]
				E_S = lineSplit[3]
				O = lineSplit[4]
				E = lineSplit[5]
				OdivS = lineSplit[6]
				OminusS = lineSplit[7]
				OdivE = lineSplit[8]
				OtimesNatLog = lineSplit[9]
				StimesNatLog = lineSplit[10]
				Pval = lineSplit[11]
				
				#wordList.append(word)
				#wordList[i] = Word(word,markovOrder,S,SdivAllS,E_S,O,E,OdivS,OminusS,OdivE,OtimesNatLog,StimesNatLog,Pval)
				wordList[word] = Word(word,markovOrder,S,SdivAllS,E_S,O,E,OdivS,OminusS,OdivE,OtimesNatLog,StimesNatLog,Pval)				
				i += 1
		IN.close()
	return wordList							
#end setWordClasses

def printWordClass(word):
	"""A function to print out values from the word class object.
	Parameters
    	----------
    	word : class object
        	Must be of class word to print out its values.
	Return
	---------
	None :
		Just prints out values to the screen.
    	"""

	print "".join(["Word = ", word.getWord()])
	print "".join(["MarkovOrder = ", str(word.getOrder())])
	print "".join(["Length = ", str(word.getLength())])
	print "".join(["S = ", str(word.getS())])
	print "".join(["S/allS = ", str(word.getSdivAllS())])
	print "".join(["E_s = ", str(word.getE_S())])
	print "".join(["O = ", str(word.getO())])
	print "".join(["E = ", str(word.getE())])
	print "".join(["O/S = ", str(word.getOdivS())])
	print "".join(["O-S = ", str(word.getOminusS())])
	print "".join(["O/E = ", str(word.getOdivE())])
	print "".join(["O*ln(O/E) = ", str(word.getOtimesNatLog())])
	print "".join(["S*ln(S/E_s = ", str(word.getStimesNatLog())])
	print "".join(["Pval = ", str(word.getPval())])
	print "".join(["Z-Score = ", str(word.getZScore())])
	print "".join(["Cluster = ", str(word.getCluster())])	

#end printWordClasses

def getSin_zscore(wordList):
	"""Gets the sin z-score for each word length in wordList.
	Parameters
    	----------
    	wordList : dictionary
		Holds all the words that were scored from wordseeker, key is the word, value is a class "word" object.
	Return
	----------
	None:
		Z-Socres are saved in each words class object.
	"""
	
	sinScore = []
	avg = 0
	stdDev = 0

	#adding Sin values to a list that belongs to its length
	for key, value in wordList.items():
		sinScore.append(value.getStimesNatLog())

	#getting the avg and stdDev
	avg = getAvg(sinScore)
	stdDev = getStdDev(sinScore)

	#appending each word so it has a z_score
	for word, value in wordList.items():

		#appending word with the Z-Score
		value.zScore = ((float(value.getStimesNatLog())-avg)/stdDev)
		
#end getZ_Score

def getAvg(values):
	"""Calcultes the average from a list of numbers.
	Parameters
    	----------
    	values : list
		A list with all the numbers to calculate an average from.
	Return
	----------
	avg : float
		The calculated average.
	"""

	if not values:
		return 0

	total = 0
	count = 0

	for num in values:
		total += float(num)
		count += 1

	avg = float(total)/float(count)

	return avg
#end getAvg

def getStdDev(values):
	"""Calcultes the standard deviation from a list of numbers.
	Parameters
    	----------
    	values : list
		A list with all the numbers to calculate an stdDev from.
	Return
	----------
	stdDev : float
		The calculated standard deviation.
	"""

	avg = getAvg(values)
	count = len(values)
	totalSqrs = 0

	i = 0
	while i < count:
		num = float(values[i]) 
		num = num - avg
		num = num**2
		values[i] = num
		i += 1

	for num in values:
		totalSqrs += num

	stdDev = math.sqrt((totalSqrs/(count-1)))
	
	return stdDev
#end getStdDev	

def filterSin_zscores(wordList, threshold):
	"""filters out words that have a z-score less than that of the threshold.
	Parameters
    	----------
    	wordList : dictionary
		Holds all the words that were scored from wordseeker, key is the word, value is a class "word" object.
	threshold : float
		Each word must have a z-score > threshold to be added to the filterList.
	Return
	----------
	filteredList : dictionary
		A dictionary containing the words that pass the z-score threshold, key is the word, value is the class "word" object.	
	"""

	filteredList = {}

	#finds all words with z-score >= threshold from input and saves them to filteredList
	for key,value in wordList.items():
		if value.getZScore() >= threshold:
			filteredList[key] = value

	return filteredList
#end filterZScores

#from http://matplotlib.org/examples/api/radar_chart.html
def radar_factory(num_vars, frame='circle'):
	"""Create a radar chart with `num_vars` axes. This function creates a RadarAxes projection and registers it.
	Parameters
    	----------
    	num_vars : int
        	Number of variables for radar chart.
    	frame : {'circle' | 'polygon'}
        	Shape of frame surrounding axes.
	"""
    	# calculate evenly-spaced axis angles
    	theta = 2*np.pi * np.linspace(0, 1-1./num_vars, num_vars)
    	# rotate theta such that the first axis is at the top
    	theta += np.pi/2

    	def draw_poly_patch(self):
       		verts = unit_poly_verts(theta)
        	return plt.Polygon(verts, closed=True, edgecolor='k')

    	def draw_circle_patch(self):
       		# unit circle centered on (0.5, 0.5)
        	return plt.Circle((0.5, 0.5), 0.5)

    	patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
    	if frame not in patch_dict:
        	raise ValueError('unknown value for `frame`: %s' % frame)

    	class RadarAxes(PolarAxes):

        	name = 'radar'
        	# use 1 line segment to connect specified points
        	RESOLUTION = 1
        	# define draw_frame method
        	draw_patch = patch_dict[frame]

        	def fill(self, *args, **kwargs):
            		"""Override fill so that line is closed by default"""
            		closed = kwargs.pop('closed', True)
            		return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

        	def plot(self, *args, **kwargs):
            		"""Override plot so that line is closed by default"""
            		lines = super(RadarAxes, self).plot(*args, **kwargs)
            		for line in lines:
                		self._close_line(line)

        	def _close_line(self, line):
            		x, y = line.get_data()
            		# FIXME: markers at x[0], y[0] get doubled-up
            		if x[0] != x[-1]:
                		x = np.concatenate((x, [x[0]]))
                		y = np.concatenate((y, [y[0]]))
                		line.set_data(x, y)

        	def set_varlabels(self, labels):
            		self.set_thetagrids(theta * 180/np.pi, labels)

        	def _gen_axes_patch(self):
            		return self.draw_patch()

        	def _gen_axes_spines(self):
            		if frame == 'circle':
                		return PolarAxes._gen_axes_spines(self)
            		# The following is a hack to get the spines (i.e. the axes frame)
            		# to draw correctly for a polygon frame.

            		# spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            		spine_type = 'circle'
            		verts = unit_poly_verts(theta)
            		# close off polygon by repeating first vertex
            		verts.append(verts[0])
            		path = Path(verts)

            		spine = Spine(self, spine_type, path)
            		spine.set_transform(self.transAxes)
            		return {'polar': spine}

    	register_projection(RadarAxes)
    	return theta
#end radarFactory()

def unit_poly_verts(theta):
    	"""Return vertices of polygon for subplot axes.
		This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    	"""
    	x0, y0, r = [0.5] * 3
    	verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    	return verts
#end unit_poly_verts(theta)

def makeSpiderChart(jobID, savedLengths, outputDir):
    	
	#setting up the directory to save the chart in
	cwd = os.getcwd()
	prefix = "".join([cwd,'/',jobID,'_SignificantWords.png'])

	#local variables
	data = []
	labels = []
	largestWords = 0
	
	#assiging values to lists for the spider chart
	for key, value in savedLengths.items():
		data.append(value)
		labels.append(key)
		if value > largestWords:
			largestWords = value

	N = len(labels)
    	theta = radar_factory(N, frame='circle')
    	spoke_labels = labels

    	#size of the file
    	fig = plt.figure(figsize=(13, 13))
    
    	fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    	color = 'b'
    	# Plot the four cases from the example data on separate axes
    	ax = fig.add_subplot(2, 2, 1, projection='radar')
    	plt.rgrids([1, (largestWords*0.25), (largestWords*0.5), (largestWords*0.75), largestWords])
    	ax.set_title(jobID, weight='bold', size='medium', position=(0.5, 1.1), horizontalalignment='center', verticalalignment='center')
    	ax.plot(theta, data, color)
    	ax.fill(theta, data, color, alpha=0.25)
    	ax.set_varlabels(spoke_labels)

    	# add legend relative to top-left plot
    	plt.subplot(2, 2, 1)
    	labels = ('Significant Words', 'Fill')
    	legend = plt.legend(labels, loc=(0.9, .95), labelspacing=0.1)
    	plt.setp(legend.get_texts(), fontsize='small')

    	plt.figtext(0.5, 0.965, 'Significant Words for Interesting Word Lengths',
        	ha='center', color='black', weight='bold', size='large')
    	
	savefig(prefix)
	shutil.move(prefix, outputDir)

#end makeSpiderChart

def editDist(str1, str2):
	"""Calculates the edit distance between two strings.
	Parameters
    	----------
    	str1, str2 : string
		Both are strings containing words to calculate an edit distance between.
	Return
	----------
	editDistance : int
		The edit distance.
	"""
	return Levenshtein.distance(str1, str2)
#end editDistance

def createDendrogram(wordCluster, filteredList):
	wordLabelList = []

	for key, value in filteredList.items():
		wordLabelList.append(key)

	#creates a dendrogram and saves it to the above file
	wordDendrogram = scipy.cluster.hierarchy.dendrogram(wordCluster, labels=wordLabelList)	
	
	return wordDendrogram
#end createDendrogram

def makeMotifsFromSeed(fastaFile, seedList, wordObjDict, outputDir, jid, numMotifs, hamDist):
	maxDist = hamDist
        motifObjList = []
	motifSeeds = []
	fileList = []
	zScoreDict = {}
	numSeqs = getFastaNumSeqs(fastaFile)

	#creating dictionary with seedword as key and zscore as value
	for seedWord in seedList:
		zScoreDict[seedWord] = wordObjDict[seedWord].zScore
	
	#create a list with top seedWords, checks for length of seedList for safety, then removes the used seedWord
	for i in range(numMotifs):
		if i < len(seedList):
			largest = 0.0
			largeWord = ''
			for seedWord, zScore in zScoreDict.items():
				if zScore > largest:
					largeWord = seedWord
					largest = zScore
			motifSeeds.append(largeWord)		
			zScoreDict.pop(largeWord, None)			

	counter = 0
	for seedWord in motifSeeds:
		newMotif = MyMotif(seedWord)			
		m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
		
		for nextWord, wordObj in wordObjDict.items():
			if len(nextWord) != len(seedWord):
				continue
			dist = editDist(nextWord, seedWord)
			if dist > maxDist:
                        	continue
			newMotif.wordList.append(nextWord)			
			for seq in wordObj.seqDict:
                                newMotif.seqDict[seq] = 1
                        #logo
                        wordCount = int(wordObj.O)
                        for i in range(wordCount):
                                m.add_instance(Seq(nextWord,m.alphabet))
                
		newMotif.seqCount = len(newMotif.seqDict)
                motifObjList.append(newMotif)
		
                #make a logo
                if counter < int(numMotifs):
                        logoName = ''.join([jid, '_', seedWord, '_seqCount_',str(round(float(newMotif.seqCount)/float(numSeqs),3)),'.png'])
                        m.weblogo(logoName)
                        fileList.append(logoName)
                counter += 1

	os.mkdir("".join([jid,'_OpenMotif_MotifLogos']))
	cwd = os.getcwd()
	motifDir = "".join([cwd,'/',jid,'_OpenMotif_MotifLogos'])

	for motif in fileList:
		motifLoc = ''.join([cwd,'/',motif])
		shutil.move(motifLoc, motifDir)
	shutil.move(motifDir, outputDir)
       
        return motifObjList
#end makeMotifsFromSeed

def occurrences(string, sub):
	count =	0 
	start = 0
	string = str(string)
	string = string.upper()
	while True:
		start = string.find(sub, start) + 1
        	if start > 0:
            		count+=1
       		else:
        		return count
#end occurrences	

def getSeqNamesForWord(word, fastaFile):
        """ given a word and a fasta file return all the seqs it occurs in """
        seqDict = {}
        wordCount = 0
       
        with open(fastaFile, "rU") as handler:
                for record in SeqIO.parse(handler, "fasta"):
                        seqId = str(record.id)
                        count = occurrences(record.seq, word)
                        if count > 0:
                                seqDict[seqId] = 1
                        wordCount += count
                        
        return seqDict
#end getSeqNamesForWord	

def addSeqDictToWords(wordList, fastaFile):
	
	for word, wordObj  in wordList.items():
		wordObj.seqDict = getSeqNamesForWord(word,fastaFile)

	return
#end addSeqDictToWords

def clusterWordsHier(filteredList, clusterMethod, outputDir, jobID, wordObjList):
	"""Performs hierarchical clustering. Check scipy clustering and http://www.cs.swarthmore.edu/~turnbull/cs67/s09/labs/lab05.pdf """
	wordDistMatrix = np.zeros( (len(filteredList), len(filteredList)) )
	wordDistList = list()	
	wordLabelList = list()
	#make a dict between the word and an index and between an inex and the word
	wordIndexDict = {}
	indexWordDict = {}
	indexCounter  = 0
	#go thru the filtered list and get the mapping between the words and their indices and vice versa
	for key, value in filteredList.items():
		wordIndexDict[key] = indexCounter
		indexWordDict[indexCounter] = key
		indexCounter += 1
		if key not in wordLabelList:
			wordLabelList.append(key)
	
	#print 'filtList len:', len(filteredList)
	
	keyCounter = 0
	#runs through the filtered list and makes a 2D dictionary
	for key, value in filteredList.items():
		#print 'key:', key
		innerKeyCounter = 0
		for innerKey, innerValue in filteredList.items():
			#print '\tinnerkey:',innerKey
			dist = editDist(key,innerKey)
			wordDistMatrix[keyCounter][innerKeyCounter] = dist
			innerKeyCounter += 1
		keyCounter += 1
		
	sqfrm = scipy.spatial.distance.squareform(wordDistMatrix)
	linkageMatrix = scipy.cluster.hierarchy.linkage(sqfrm, method=clusterMethod)
	
	wordDendrogram = scipy.cluster.hierarchy.dendrogram(linkageMatrix, labels=wordLabelList, leaf_font_size = 7)
	cwd = os.getcwd()
	prefix = "".join([cwd,'/',jobID,'_dendrogram.png'])
	savefig(prefix)
	shutil.move(prefix, outputDir)

	#set the max number of clusters to be no more than half the number of words in the filtered list. 
	tThreshold =  math.ceil(len(wordLabelList))/2	
	clusterLabels = scipy.cluster.hierarchy.fcluster(linkageMatrix, t = tThreshold, criterion = 'maxclust')
	#print 'clusterLabels:', clusterLabels
	#print 'clusterWords :',wordLabelList 
	
	#make the cluster dict. Each cluster id will hold a list of words that belong to it 
	clusterDict = {}
	wordIndex = 0
	for clusterId in clusterLabels:
		#initialize the dict of lists if not already
		if clusterId not in clusterDict:
			clusterDict[clusterId] = list()
		#get the word form its index
		word = indexWordDict[wordIndex]
		clusterDict[clusterId].append(word)
		wordIndex += 1
		
	#check the cluster and pick the word with highest seq sln score
	seedList = list()
	for clusterId in clusterDict:
		#print 'clusterId:', clusterId
		maxScore = 0
		chosenWord = ''
		for word in clusterDict[clusterId]:
			score = wordObjList[word].zScore
			#print '\tword:',word,' score:',score
			if score > maxScore:
				maxScore = score
				chosenWord = word
			#print '\t',word
		seedList.append(chosenWord)
	
	#print 'seedList:', seedList
	
	return seedList
#end clusterWordsHier

def visMotifs(fastaFile, jobID, reverseOption, numMotifs, motifObjList, outputDir):
        """ visualize the motifs """
        #make a motif file in the format for the visualizer
	cwd = os.getcwd()
        motifFileName= ''.join([cwd,'/',jobID,'_visMotifs'])
        motifFile = open(motifFileName, 'wb')
	numMotifs = int(numMotifs)
        counter = 0
        for motifObj in motifObjList:
                seedWord = motifObj.seedWord
                motifFile.write('>seedWord:' + seedWord + '\n' )
                for word in motifObj.wordList:
                        motifFile.write(word + '\n')
                if counter == numMotifs - 1:
                        break
                counter += 1

        motifFile.close()
     
        htmlFileName = visualize_motifs.callVis(fastaFile, motifFileName, jobID, reverseOption)
 
        #move the the other xslt files too
        cismlFile = ''.join([jobID,'_cisml.xml'])
        xslFile = ''.join([jobID,'_mcast-to-html.xsl'])
        xmlFile = ''.join([jobID,'_motifsXML.xml'])
   
	#moving all files to output dir and deleting unused ones
	os.remove(motifFileName)
	os.remove(cismlFile)
	os.remove(xslFile)
	os.remove(xmlFile)
	shutil.move(htmlFileName, outputDir)
	return
#end visMotifs

def makeMotifFiles(fastaFile, motifObjList, outputDir, jobID):

	#adding seqCov to motifs object
	numSeqs = getFastaNumSeqs(fastaFile)
	for motifObj in motifObjList:
		motifObj.seqCov = float(motifObj.seqCount)/float(numSeqs)

	#making csv which has motifs and sequences that it occurs in
	seqDict = getSeqDict(fastaFile)
	cwd = os.getcwd()
        motifSeqsFile = ''.join([cwd,'/',jobID,'_motifSeqs.csv'])
	OUT = open(motifSeqsFile, 'wb')
	for motifObj in motifObjList:
		OUT.write('>' + motifObj.seedWord+'\n')
		for key, value in motifObj.seqDict.items():
			OUT.write(seqDict[key]+'\n')
	OUT.close()
	shutil.move(motifSeqsFile, outputDir)
	
	#making csv for motifs and words and which seqs that they occur in
	motifDetFile = ''.join([cwd,'/',jobID,'_motifDetSeqs.csv'])
	OUT = open(motifDetFile, 'wb')
	for motifObj in motifObjList:
		OUT.write('>' + motifObj.seedWord+'\n')
		for word in motifObj.wordList:
			OUT.write(word + '\n')
			seqs = getSeqNamesForWord(word, fastaFile)
			for key, value in seqs.items():
				OUT.write(seqDict[key]+'\n')
	OUT.close()
	shutil.move(motifDetFile, outputDir)
	return
#end makeMotifFiles

def findPFM(jobID, motifObjList, wordObjDict, numMotifs, outputDir):
    	""" find the PFM(Position Frequency Matrix) for the top motifs """
    	#write the words of motif in Jaspar site format like here: https://github.com/biopython/biopython/blob/master/Doc/cookbook/motif/Arnt.sites
    	alphaList = ['A', 'C', 'G', 'T']
    	siteFileName = ''.join([jobID,'_jasparWordFile'])
    	pfmFileName = ''.join([jobID,'_PFM'])
    	pfmFile = open(pfmFileName, 'wb')
    	counter = 1
    	#write the words
    	for motifObj in motifObjList:
    		seedWord = motifObj.seedWord
        	siteFile = open(siteFileName, 'wb')
        	for word in motifObj.wordList:
            		wordCount = wordObjDict[word].O
            		for i in range(int(wordCount)):
                		siteFile.write('>site ' + str(counter) + '\n' + word + '\n')
                		counter += 1
        	siteFile.close()
        
        	srf = Motif.read(open(siteFileName),'jaspar-sites')
        	srf.make_counts_from_instances()
        	pfmFile.write('\n>' + seedWord + '\n')
        	for alpha in alphaList:
        		pfmFile.write(alpha + ' ' + str(srf.counts[alpha]) + '\n')
	     
	shutil.move(pfmFileName, outputDir) 
    	os.remove(siteFileName)
    	pfmFile.close()
	return
#end findPFM

######################################################################################################################

######## MAIN  #######################################################################################################
def main():

	start = time.time()
	#parsing command line args
	parser = argparse.ArgumentParser(description='Enter command line arguments.')
	parser.add_argument("-P", "--pSeq", help="Enter Positive Input Sequence (Fasta)")
	#parser.add_argument("-N", "--nSeq", help="Enter Negative Input Sequence (Fasta)")
	parser.add_argument("-maxLen", "--maxWordLen", help="Enter Maximum Word Length", type=int)
	parser.add_argument("-minWLen", "--minWordLen", help="Enter Minimum Word Length", type=int)
	parser.add_argument("-minSCov", "--minSeqCov", help="Enter Minimum Sequence Coverage % (10 = 10% Sequence Coverage)", type=int)
	parser.add_argument("-minWCnt", "--minWordCnt", help="Enter Minimum Word Count", type=int)
	parser.add_argument("-zScrTh", "--zScoreThreshold", help="Enter Z-Score Threshold", type=float)
	#parser.add_argument("-hamDist", "--hamDist", help="Enter Hamming Distance")
	#parser.add_argument("-di-Mismatch", "--diMisValue", help="Enter Di-mismatch value")
	parser.add_argument("-numMotifs", "--numMotifs", help="Number of Motifs to Report", type=int)
	parser.add_argument("-jid", "--jobID", help="Enter a Job ID")
	args = parser.parse_args()
	
	markovOrder = 2
	hamDist = 1

	#creating a directory to store all output files
	os.mkdir("".join([args.jobID,'_OpenMotif_OutputFiles']))
	cwd = os.getcwd()
	outputDir = "".join([cwd,'/',args.jobID,'_OpenMotif_OutputFiles'])
	
	#calling wordseeker, returns location of the count file
	seekerOutput = callWordSeeker(args.pSeq, args.jobID, args.minWordLen, args.maxWordLen)
	
	end = time.time()
	print "non scoring" + str((float(end) - float(start)))
	#getting number of sequences in the fasta file
	numSeqs = getFastaNumSeqs(args.pSeq)
	#doing the math on how many sequences a word must cover to be interesting 
	numSeqs = math.ceil(numSeqs * (float(args.minSeqCov)/100))
	
	#getting a dictionary that maps each word length with how many words cover the min count
	savedLengths = getWordCount(seekerOutput, numSeqs, args.minWordCnt)
	
	#calling wordseeker again with only lenths we found interesting and with scoing on, returns the directory with the output
	scoringOutput = callWordSeekerScoring(args.pSeq, args.jobID, markovOrder, savedLengths)		
	end = time.time()
	print "scoring " + str((float(end) - float(start)))	

	#returns a list with all the words in class form, see class above for form and use print function to show a word
	wordList = setWordClasses(scoringOutput, markovOrder)
	shutil.move(seekerOutput, outputDir)
	shutil.move(scoringOutput, outputDir)
	end = time.time()
	print "classes " + str((float(end) - float(start)))

	#adds seqDictionary to each word that gives which sequences a word occurs in
	addSeqDictToWords(wordList, args.pSeq)
	end = time.time()
	print "add seqDict " + str((float(end) - float(start)))
	
	#getting z-score for each word in the wordList, z-scores are added as an attribute to each word object
	getSin_zscore(wordList)
	end = time.time()
	print "z_score " + str((float(end) - float(start)))
	
	#getting the words with z-Scores >= Z-ScoreThreshold, returns a dictionary with only those words
	filteredList = filterSin_zscores(wordList, args.zScoreThreshold)
		
	#perform clustering, make dendrogram and return seedwords
	clusterMethod = "average"
	seedWords = clusterWordsHier(filteredList, clusterMethod, outputDir, args.jobID, wordList)	
	end = time.time()
	print "cluster " + str((float(end) - float(start)))

	#making motifs
	motifObjList = makeMotifsFromSeed(args.pSeq, seedWords, wordList, outputDir, args.jobID, args.numMotifs, hamDist)
	end = time.time()
	print "motifs " + str((float(end) - float(start)))

	#calling visualizer
	reverseOption = 1
	visMotifs(args.pSeq, args.jobID, reverseOption, args.numMotifs, motifObjList, outputDir)
	end = time.time()
	print "visualize " +str((float(end) - float(start)))

	#adding seqCov for each motif to its Object and making csv files motifSeqs and motifDetSeqs
	makeMotifFiles(args.pSeq, motifObjList, outputDir, args.jobID)
	end = time.time()
	print "motifCSVfiles " +str((float(end) - float(start)))	

	#makes PMF file to use on jaspar site
	findPFM(args.jobID, motifObjList, wordList, args.numMotifs, outputDir)
	end = time.time()
	print "PMF " + str((float(end) - float(start)))

	#making spider chart for significant words
	#makeSpiderChart(args.jobID, savedLengths, outputDir)
	
#end main

if __name__ == '__main__':
     """ MAIN """
     main()
