#!/usr/bin/python
#####################################################################################################################
#
#	Program:	memeVisualizer-v2.py
#	Author:		Max Rego
#	Email:		mr255509@ohio.edu
#
#	Description:	This program takes in a seedWord file, finds the occurances in a fasta file
#			and visualizes the data in an html file.
#
#	Date:		4-25-2013
#
#####################################################################################################################

#imports
import sys
import libxml2
import libxslt
from Bio import SeqIO

### FUNCTIONS  #######################################################################################################

def processWords(inputWordsFile):
   	
	#opening file
	IN = open(inputWordsFile, 'r')
	
	#variables	
	motifHash = {}		#holds the final dictionary, called it a hash to keep old variables names
	wordsArr = []		#holds all words under seedWord
	seedWord = ""		#holds the seedWord
	wordCount = 0		#if 0 then its the first run, if not then add to hash

	for line in IN:
		#remove new line		
		line = line.rstrip()
		#if line is a seedword
		if line.find('>seedWord:') != -1:
			if (wordCount != 0):
				if wordsArr:
					motifHash[seedWord] = wordsArr
				#end
			#end
			#reset wordsArr and get new seedWord 	
			wordsArr = []	
			seedWord = line.split(':')
			seedWord = seedWord[1]
			#increase wordcount
			wordCount += 1

		#if line is not a seedword add it to wordsArr
		else:
			wordsArr.append(line)
		#end
	#end for
	#add last seedWord if it exists
	if wordsArr:	
		motifHash[seedWord] = wordsArr				
	#end			

	return motifHash
#end processWords

def memeVisual(jid, inputSeqFile, motifHash):

	
	numSeqs = 0	#total number of headers >seq
	numBases = 0	#total number of letters in all sequences

	#getting fasta data into a list
	IN = open(inputSeqFile, "rU")
	seqList = list(SeqIO.parse(IN, "fasta"))
	IN.close()
	#print records[0].id     #seq header
	#print records[0].seq	 #sequence

	#setting numSeqs and numBases
	numSeqs = len(seqList)
	for i in seqList:
		numBases = numBases + len(i.seq)
	#end 

	motifFile = jid + "_motifsXML.xml"

    	MOTIFOUT = open(motifFile, 'w')
   	
	temp = ""
    	temp = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + '<' + jid + '_motifsXML version="4.8.1" release="0 EST 20">\n'
      	temp +='xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n' + 'xsi:schemaLocation=  xmlns:mhmmscan="http://noble.gs.washington.edu/schema/fimo"\n'
    	temp += '>\n' + '<command-line>./openMotif</command-line>\n' + '<sequence-data num-sequences="' + str(numSeqs) + '" num-residues="' + str(numBases) + '" />\n'
    	temp += '<alphabet>nucleotide</alphabet>\n'
	MOTIFOUT.write(temp)	

	for k, v in motifHash.iteritems():
     		#print k, v
		length = len(k)
        	wordString=""
		i = length
        	while(i > 0):
        		wordString += " "
			i -= 1
        	#end
        	MOTIFOUT.write('<motif name="' + k + '" width="' + str(length) + '"  best-possible-match="' + wordString + '"/>\n')
    	#end

    	MOTIFOUT.write('<cisml-file>' + jid + '_cisml.xml</cisml-file>\n')
    	MOTIFOUT.write('</' + jid + '_motifsXML>\n')
    	MOTIFOUT.close()
	#done with _motifsXML.xml

	#start _cisml.xml file
	cismlFile = jid + "_cisml.xml"
	CISMLOUT = open(cismlFile, 'w')
	
	temp = ""
	temp = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + '<?xml-stylesheet type="text/xsl" href="mhmmscan-to-html.xsl"?>\n'
	temp += '<cis-element-search\n' + 'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n' + 'xsi:schemaLocation="http://zlab.bu.edu/schema/cisml cisml.xsd"\n'
      	temp += 'xmlns="http://zlab.bu.edu/schema/cisml"\n' + 'xmlns:mem="http://noble.gs.washington.edu/meme">\n' + '<program-name>MEME-visualizer</program-name>\n'
	CISMLOUT.write(temp)
	
	patternFile = inputSeqFile
    	inputFile = inputSeqFile
	temp = ""
    	temp = '<parameters>\n' + '<pattern-file>' + patternFile + '</pattern-file>\n' + '<sequence-file>' + inputFile + '</sequence-file>\n' 
        temp += '<pattern-pvalue-cutoff>0.0005</pattern-pvalue-cutoff>\n' + '<sequence-pvalue-cutoff>1</sequence-pvalue-cutoff>\n' + '</parameters>\n'
	CISMLOUT.write(temp)

	clusterIndex = 1
	#for each header in fasta file -> i
	for i in seqList:
		seqID = i.id
		seqString = i.seq
		metric = "motifSeqCov"
		counter = 0
		pValue = "nan"
		seqCov = 16.1782
		CISMLOUT.write('<multi-pattern-scan score="' + str(seqCov) + '" pvalue="' + pValue + '">\n')
		#iterate through motifHash, seedword = k, words list = v		
		for seedWord, wordList in motifHash.iteritems():
			#iterate through each word in word list
			for word in wordList:
				
				offset = 0
				result = seqString.find(word, offset)
								
				while (result != -1):				
					if seqID.find('REVERSE') != -1:
                    				
						#whats the point of revResult?
						revResult = len(seqString) - result - len(word)
                        			CISMLOUT.write('<pattern accession="'+seedWord+'" name="'+seedWord+'">\n' + '\t<scanned-sequence accession="'+seqID+'" name="'+seqID+'">\n')

                        			start = result + 1
                       				stop = result + len(word)
                        			score = 0
                        			
						temp = ""	
                        			temp = '\t\t<matched-element start="'+str(stop)+'" stop="'+str(start)+'" pvalue="'+str(score)+'">\n'
                        			temp += '\t\t\t<sequence>'+seedWord+'</sequence>\n'+'\t\t</matched-element>\n'+'\t</scanned-sequence>\n'+'</pattern>\n\n'
						CISMLOUT.write(temp)

					else:
						CISMLOUT.write('<pattern accession="'+seedWord+'" name="'+seedWord+'">\n' + '\t<scanned-sequence accession="'+seqID+'" name="'+seqID+'">\n')

                        			start = result + 1
                       				stop = result + len(word)
                        			score = 0
                        			
						temp = ""	
                        			temp = '\t\t<matched-element start="'+str(start)+'" stop="'+str(stop)+'" pvalue="'+str(score)+'">\n'
                        			temp += '\t\t\t<sequence>'+seedWord+'</sequence>\n'+'\t\t</matched-element>\n'+'\t</scanned-sequence>\n'+'</pattern>\n\n'
						CISMLOUT.write(temp)						
					#end if reverse or not
					offset = result + 1
					result = seqString.find(word, offset)
				#end while result != -1
			#end each word
		#end each seedword
		
		#write the sequnce itself
        	clusterID = "cluster-" + str(clusterIndex)
        	seqStart = 1
        	seqEnd = len(seqString)
        	CISMLOUT.write('<mem:match cluster-id="'+clusterID+'" seq-name="'+seqID+'" start="'+str(seqStart)+'" stop="'+str(seqEnd)+'" evalue="nan" qvalue="nan">'+str(seqString)+'\n')
        	CISMLOUT.write('</mem:match>\n' + '</multi-pattern-scan>\n\n')
		clusterIndex += 1
	#end each sequence	
	
	CISMLOUT.write('</cis-element-search>\n')
    	CISMLOUT.close()

	#Converting to html
	htmlOut = jid + "_visual.html"
	str(htmlOut)	
	xsltFile = writeXSLTFile(jid, motifFile, cismlFile)
	style = getXSLT(xsltFile)
	doc = libxml2.parseFile(motifFile)
	result = style.applyStylesheet(doc, None)

	OUT = open(htmlOut, 'w')
	OUT.write(str(result))
	OUT.close()

	return
#end memeVisual

def writeXSLTFile(jid, motifFile, cismlFile):

	outFile = jid + "_mcast-to-html.xsl"
	motifFileName = motifFile.replace(".xml", "")
	cismlFileName = cismlFile.replace(".xml", "")	
	OUT = open(outFile, 'w')
	IN = open("mcast-to-html.xsl", 'r')
	for line in IN:
		if line.find('xmlns:cis="http://zlab.bu.edu/schema/cisml') != -1:	
			
			OUT.write(line)			
		
		elif line.find('cisml') != -1:
			
			line = line.replace('cisml', cismlFileName)
			OUT.write(line)
			
		elif line.find('mcast') != -1:
			if line.find('(mcast)') == -1:
				line = line.replace('mcast', motifFileName)
			OUT.write(line)
			
		else:
			OUT.write(line)
		#end if
	#end for
	OUT.close()
	IN.close()

	return outFile

#end writeXSLTFile

def getXSLT(xsl_filename):
    
	# parse the stylesheet xml file into doc object
	styledoc = libxml2.parseFile(xsl_filename)

	# process the doc object as xslt
	style = libxslt.parseStylesheetDoc(styledoc)

	return style

#end getXSLT

####  MAIN  ####################################################################################################################

if __name__ == '__main__':

	#command line inputs
	inputSeqFile = sys.argv[1]
	inputWordsFile = sys.argv[2]
	jid = sys.argv[3]
	numArgs = len(sys.argv)

	#variables
	motifHash = {}	#holds the seedWords hash

	#make sure correct num args, if under 3 then traceback will occur, over 3 and this statement will exit program
	if numArgs != 4:
		print "Incorrect number of arguements"
		sys.exit(1) # Or something that calls sys.exit()
	#end

	#get motifHash from function
	motifHash = processWords(inputWordsFile)
	
	memeVisual(jid, inputSeqFile, motifHash)






