#!/usr/bin/env python

'''
Script created by Christopher Laumer to select longest ORF from TransDecoder output. Worked with older versions of Trinity output, such as 20140413.
This script has been updated by Tauana Cunha (http://github.com/tauanajc) on August 2016 to work with the current version of Trinity (2.2.0).
The headers of the transcripts are different, and now include DN[0-9]+ as part of the gene identity.
The older version of the script would greatly reduce the number of final genes beacuse it was lumping as one gene many different things.
New headers are: gene (DN\d*_c\d*_g\d*:), isoform (DN\d*_c\d*_g\d*\_i\d*:)
'''


from Bio import SeqIO
import re
import argparse

def parser():
	args = argparse.ArgumentParser(description='This program takes an assembly of peptides and outputs the longest ORF for each gene.')
	args.add_argument('-i', '--input', required=True, help='The name of the fasta input file.')
	args.add_argument('-o', '--output', required=True, help='The name of the fasta output file you wish to create.')
	args.add_argument('-t', '--test', help='Choose whether to run the test suite.', action='store_true')
	args.add_argument('-l', '--longest', help='A flag that allows you to retain only one protein per subcomponent', action='store_true')
	args = args.parse_args()
	return args

if parser().test:
	print 'Test suite turned on. This could take a while. If "-l" is on, check 1 should fail (-l chooses only one peptide per transcript).'

infile = open(parser().input, 'rb')
outfile = open(parser().output,'wb')

def write_output(dict_final, out):
	for each in dict_final:
		out.write('>' + str(each) + '\n')
		out.write(str(dict_final[each]) + '\n\n')
	out.close()

transcripts = SeqIO.to_dict(SeqIO.parse(infile, "fasta"), key_function = lambda rec : rec.description)

def simplify(trans_dict):
	for i in trans_dict:
		trans_dict[i] = str(trans_dict[i].seq)
	return trans_dict

def strip_stop_codons(dict):
	'''
	In a dictionary representation of transdecoder output, removes all '*' characters (stop codons)
	'''
	for i in dict:
		dict[i] = dict[i].strip('*')
	return dict

def sorted_headers(trans_dict):
	'''
	Takes as input an unsorted sequence dictionary of predicted ORFs from trinity transcripts.
	Outputs a list of the keys in that dictionary, sorted lexicographically by their component ID.
	'''
	headers = trans_dict.keys()
	transcript_id = re.compile('DN\d*_c\d*_g\d*\_i\d*:') # Added DN\d*_ for new Trinity gene hierarchy
	for each in range(len(headers)):
		id = transcript_id.search(headers[each]).group()
		headers[each] = (id,headers[each])
	headers.sort()
	for i in range(len(headers)): #strips 
		headers[i] = headers[i][1] 
	return headers

def split_by_isoform(dict, L): 
	'''
	Takes as input a.) dict, a dictionary representation of all transcripts in a trinity assembly
	and b.) L, a list of the keys in dict, sorted lexicographically by comp_c_seq
	and c.) a compiled regular expression object to match comp_c (subcomponent ID only)
	Returns a dict composed of all predicted ORFs from the seq ID that contains the longest ORF in its subcomponent
	'''
	subcomp_id = re.compile('DN\d*_c\d*_g\d*') # Added DN\d*_ for new Trinity gene hierarchy
	transcript_id = re.compile('DN\d*_c\d*_g\d*\_i\d*:') # Added DN\d*_ for new Trinity gene hierarchy
	subcomponent = {}
	longest = {}
	P = []
	for i in range(len(L)):
		q = subcomp_id.search(L[i]).group()
		P.append(q)
	for j in range(len(P)):
		if P.count(P[j]) == 1:
			longest[L[j]] = dict[L[j]]
	for i in range(len(L)-1):
		if subcomp_id.search(L[i]).group() == subcomp_id.search(L[i+1]).group(): #checks if the ith and i+1th headers have the same subcomp ID
			subcomponent[L[i]] = dict[L[i]]
			subcomponent[L[i+1]] = dict[L[i+1]] 	
			try:
				if subcomp_id.search(L[i+1]).group() != subcomp_id.search(L[i+2]).group():
					string = id_longest_seq(subcomponent, transcript_id) #stores the trinity ID of the longest ORF in subcomponent
					add_to_final(string, subcomponent, longest) #adds all seqs containing str to longest
					subcomponent = {}
			except:
				string = id_longest_seq(subcomponent, transcript_id)
				add_to_final(string, subcomponent, longest)
				subcomponent = {} #clears subcomponent dictionary
	return longest
		
def id_longest_seq(tempdict, re):
	'''
	Input is a.) a sequence dictionary composed of all elements with the same subcomponent ID
	and b.) a compiled regular expression object to match comp\d*_c\d*_seq\d*
	Returns a string with the comp_c_seq ID of the longest sequence in the tempdict
	'''
	longest = 0	
	for i in tempdict:
		if len(tempdict[i]) > longest:
			longest = len(tempdict[i])
			string = i
	match = re.search(string)
	return match.group()

def add_to_final(component_id,tempdict,longestdict):
	'''
	Takes as input a.) component_id, the trinity ID (of form comp_c_seq) of the longest sequence in tempdict
	b.) tempdict, a dictionary containing all peptides with the same subcomponent ID
	c.) longestdict, a dictionary used to accumulate all predicted ORFs from the transcript that has the longest such ORF in a subcomponent
	Mutates longestdict to contain all ORFs with str as a substring of their key
	'''
	if parser().longest:
		for i in tempdict:
			if tempdict[i] == max(tempdict.values(), key=len) and component_id in i:
				longestdict[i] = tempdict[i]
				pass
	else:
		for i in tempdict:
			if component_id in i:
				longestdict[i] = tempdict[i]

def test(in_dict, out_dict):
	subcomp_id = re.compile('DN\d*_c\d*_g\d*') # Added DN\d*_ for new Trinity gene hierarchy
	transcript_id = re.compile('DN\d*_c\d*_g\d*\_i\d*:') # Added DN\d*_ for new Trinity gene hierarchy
	def check1(in_dict, out_dict):
		'''
		For each peptide in the filtered dictionary, check to see if the transcript that this peptide comes from has the same number of peptides in the input dictionary as in the filtered dictionary.
		Failure indicates that some of the peptides from the selected transcript are missing (or possibly, were duplicated).
		'''
		for k in out_dict:
			m = transcript_id.search(k)
			count1 = 0
			count2 = 0 
			for i in in_dict.keys():
				if m.group() in i:
					count1 += 1
			for j in out_dict.keys():
				if m.group() in j:
					count2 += 1
			if count1 != count2:
				print 'Check 1 failed! Found at least one transcript in the filtered fasta that has an unequal number of peptides than present for that transcript in the input file.'
				return False
		print 'Check 1 passed! Each component ID present in the filtered fasta has the same number of peptides as in the input fasta.'
	def check2(in_dict, out_dict):
		'''
		Checks to see if there are any singleton subcomponents present in the input dictionary that aren't in the output dictionary.
		'''
		for k in in_dict:
			count = 0
			m = subcomp_id.search(k)
			for j in in_dict.keys():
				if m.group() in j:
					count += 1
					p = j
			if count == 1:
				if p not in out_dict:
					print 'Check 2 failed! Found a singleton subcomponent in input fasta that didn\'t make it into the filtered fasta.'
					return False
		print 'Check 2 passed! No singleton subcomponents present in the input file were missing from the filtered fasta.'
	def check3(in_dict, out_dict):
		'''
		Checks to make sure that the peptides in the filtered dictionary all come from the transcript which contains the longest predicted peptide in its subcomponent.
		'''
		for k in in_dict:
			longest1 = 0
			longest2 = 0
			if k not in out_dict:
				m = transcript_id.search(k)
				p = subcomp_id.search(k)
				if m and len(in_dict[k]) > longest1:
					longest1 = len(in_dict[k])
				for i in out_dict.keys():
					if p.group() in i and len(out_dict[i]) > longest2:
						longest2 = len(out_dict[i])
			if longest1 > longest2:
				print 'Check 3 failed! Found a transcript in the input fasta that makes a longer peptide than the one made by the transcript from the same subcomponent in the filtered fasta.'
				return False
		print 'Check 3 passed! The peptides in the filtered fasta all come from the transcript that makes the longest peptide in its subcomponent.'
	a = check1(in_dict, out_dict)
	b = check2(in_dict, out_dict)
	c = check3(in_dict, out_dict)
	if False in [a, b, c]:
		return False
	else:
		return True	
	
#def pathology1(in_dict, out_dict):
#	for i in range(100):
#		if in_dict.keys()[i] not in out_dict:
#			out_dict[in_dict.keys()[i]] = in_dict[in_dict.keys()[i]]

#def pathology2(out_dict):
#	for i in range(100):
#			del out_dict[out_dict.keys()[i]]

transcripts = simplify(transcripts)
transcripts = strip_stop_codons(transcripts)
print 'The process starts with', len(transcripts), 'translated transcripts'
sorted = sorted_headers(transcripts)
out_transcripts = split_by_isoform(transcripts, sorted)
print 'The process ends with', len(out_transcripts), 'peptides (longest isoform of each)'
write_output(out_transcripts, outfile)
if parser().test:
	test(transcripts, out_transcripts)

