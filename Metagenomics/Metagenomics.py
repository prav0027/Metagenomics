from matplotlib import pyplot as plt
import numpy as np


def create_histogram(dataList):
    npData = np.array(dataList)
    fig, ax = plt.subplots(figsize =(10, 7))
    ax.hist(npData, bins = 40)
    plt.show()

def read_genome(ref_fn):
    with open(ref_fn, 'r') as f:
        next(f)
        output_reference = ''
        for line in f:
            if '>' in line:
                continue
            line = line.strip()
            output_reference += line  # We append each line to the output reference string.
    return output_reference

def read_reads(read_fn):
    all_reads = []
    #readCount = 0
    with open(read_fn, 'r') as f:
        for line in f:
            if '>' in line: #and readCount < maxReads:
                continue
            line = line.strip()
            all_reads.append(line)
            #readCount+=1
    return all_reads

def listGenomes(read):
    outputList = list()
    for genomeInt in range(0, 1000):
        genomeString = 'project4a_10000_genome_'
        genomeString = genomeString + str(genomeInt) + '.fasta'
        if read in read_genome(genomeString):
            outputList.append(genomeInt)
    return outputList

def getSmallestKMer(read, k):
    currStr = read[:k]
    lexMin = currStr
    for i in range(k, len(read)):
        currStr = currStr[1 : k] + read[i]
        if (lexMin >currStr):
            lexMin = currStr
    return lexMin

readsListTemp = read_reads('project4a_10000_reads.fasta')
readsList = readsListTemp[:1000]
print("Length: ", len(readsList))
genomesList = list()
for read in readsList:
    smallKmer = getSmallestKMer(read, 10)
    for i in listGenomes(smallKmer):
        genomesList.append(i)

genomeslisttemp = list()

count=0
setGenomes = set(genomesList)
#print("134: ", setGenomes)
while count<10:
    maxGenome = max(setGenomes, key=genomesList.count)
    genomeslisttemp.append(maxGenome)
    setGenomes.remove(maxGenome)
    count+=1

#create_histogram(genomesList)

#358
#377
#622
#146
#489
#631
#785
#550
#671
#118

genomesList = genomeslisttemp

#copied from 1b
class BWT:
	def __init__(self, reference):
		suffixRotation = list()
		reverseSuffixRotation = list()
		suffix_array = list()
		bwt = list()
		C = dict()
		Occ = dict()
		Occ_reverse = dict()
		alphabet = set()
		referenceReverse = reference[::-1]
		
		alphabet.add('C')
		alphabet.add('G')
		alphabet.add('T')
		alphabet.add('A')
		
		for letter in alphabet:
			C[letter] = 0
			Occ[letter] = list()# from website: in Occ, each character has an associated list of integer values
			Occ_reverse[letter] = list()
	
		#add ending character to reference
		reference = reference+"$"
		referenceReverse = referenceReverse+"$"

		#create suffix combinations of the reference and reverse reference
		#don't store everything??? how do you not store everything and then still be able to sort it
		for i in range(len(reference)):
			if i+200 > len(reference):
				new_rotation = "%s%s" % (reference[i:],reference[0:200-len(reference)+i])
				new_rotation_reverse = "%s%s" % (referenceReverse[i:],referenceReverse[0:200-len(referenceReverse)+i])
			else:
				new_rotation = reference[i:i+200] #"%s%s" % (reference[i:i+200],reference[0:i])
				new_rotation_reverse = referenceReverse[i:i+200]
			suffixObj = Suffix(new_rotation,i)
			suffixRotation.append(suffixObj)
			
			
			#new_rotation_reverse = "%s%s" % (referenceReverse[i:],referenceReverse[0:i])
			suffixObj_rev = Suffix(new_rotation_reverse,i)
			reverseSuffixRotation.append(suffixObj_rev)
		
			#number of characters that are smaller than each char
			if reference[i]!='$':
				for char in alphabet:
					if reference[i] < char:
						C[char] = C[char] + 1	
		
#		print("sorting")
		suffixRotation.sort(key=textKey)
		reverseSuffixRotation.sort(key=textKey)

		#print("finished sorting")
					
		for i in suffixRotation:
			suffix_array.append(i.pos)
		
			#create the Occ
			for letter in alphabet:
				if len(Occ[letter]) == 0:
					prev = 0
				else:
					prev = Occ[letter][-1]
				text = "%s%s" % (reference[i.pos:],reference[0:i.pos])
				if text[-1:]==letter:
					Occ[letter].append(prev+1)
				else:
					Occ[letter].append(prev)

		#print("created suffix array")			
		#record results into the suffix array and the BWT array and calculate Occ
		for i in reverseSuffixRotation:
			#create the Occ
			for letter in alphabet:
				if len(Occ_reverse[letter]) == 0:
					prev = 0
				else:
					prev = Occ_reverse[letter][-1]
				text = "%s%s" % (referenceReverse[i.pos:],referenceReverse[0:i.pos])
				if text[-1:]==letter:
					Occ_reverse[letter].append(prev+1)
				else:
					Occ_reverse[letter].append(prev)					
		
		self.SA = suffix_array
		self.C = C
		self.Occ = Occ
		self.Occ_reverse = Occ_reverse
		self.n = len(reference)
		self.newArray = list()
		self.alphabet = alphabet

	def getMatchIndex(self,read,num_differences):
		if num_differences == 0:
			return self.matchNoDifferences(read)
		else:
			return self.matchWithDifferences(read,num_differences)


	def matchNoDifferences(self, read):
		#read = read.lower()
		i = 0
		j = self.n - 1
		
		for x in range(len(read)):
			newChar = read[-x-1]
			newI = self.C[newChar] + self.OCC(newChar,i-1) + 1
			newJ = self.C[newChar] + self.OCC(newChar,j)
			i = newI
			j = newJ
		matches = self.SA[i:j+1]
		return matches

	#threshold is number of allowed mutations
	def matchWithDifferences(self,read,threshold):
		self.calculate_distance(read)
		SA_indeces = self.recursionMatchingWithDifferences(read, len(read)-1, threshold, 0, self.n-1)
		return [self.SA[x] for x in SA_indeces]

	def recursionMatchingWithDifferences(self,read,i,threshold,m,n):
		tempset = set()
		resultZ = dict()
		resultDict = dict()
		if (threshold < self.get_D(i)):
			return set()
		if i < 0:
			for j in range(m,n+1):
				tempset.add(j)
				resultZ[j] = threshold
				resultDict[j] = resultZ
			return tempset
			
		result = set()
		#insertion
		resultAdd = self.recursionMatchingWithDifferences(read,i-1,threshold-insertion_penalty,m,n)
		result = result.union(resultAdd)
		for char in self.alphabet:
			newM = self.C[char] + self.OCC(char,m-1) + 1 
			newN = self.C[char] + self.OCC(char,n)
			if newM <= newN:#if the substring was found
				resultAdd = self.recursionMatchingWithDifferences(read,i,threshold-deletion_penalty,newM,newN)
				result = result.union(resultAdd)
				if char == read[i]:#char aligned correctly
					resultAdd = self.recursionMatchingWithDifferences(read,i-1,threshold,newM,newN)
					result = result.union(resultAdd)
				else:#unaligned
					resultAdd = self.recursionMatchingWithDifferences(read,i-1,threshold-mismatch_penalty,newM,newN)
					result = result.union(resultAdd)
		return result

	def calculate_distance(self,read):
		m = 0
		l = self.n-1
		distance = 0
		self.newArray = list()
		for i in range(len(read)):
			m = self.C[read[i]] + self.OCC(read[i],m-1,reverse=True) + 1
			l = self.C[read[i]] + self.OCC(read[i],l,reverse=True)
			if m > l:
				m = 0
				l = self.n - 1
				distance = distance + 1
			self.newArray.append(distance)

	def OCC(self,char,index,reverse=False):
		if index < 0:
			return 0
		else:
			if reverse:
				return self.Occ_reverse[char][index]
			else:
				return self.Occ[char][index]
	
	def get_D(self,index):
		if index < 0:
			return 0
		else:
			return self.newArray[index]

class Suffix:
	def __init__(self, text, position):
		self.text = text
		self.pos = position

def textKey( a ): return a.text

resultsDict = dict()

for i in genomesList:
	referenceStr = 'project4a_10000_genome_'
	referenceStr = referenceStr + str(i) + '.fasta'

	reference = read_genome(referenceStr)
	#reference = reference.lower()
	bw = BWT(reference)

	reads = readsListTemp

	insertion_penalty = 1
	deletion_penalty = 1
	mismatch_penalty = 1

	matchIndices = []

	count = 0
	for query in reads:
		difference_threshold = 0
		#query = query.lower()
		matches = []
		while len(matches) == 0 and difference_threshold <= 4:
			matches = bw.getMatchIndex(query,difference_threshold)
			difference_threshold+=1
		if len(matches) > 0:
			if count in list(resultsDict.keys()):
				if difference_threshold < resultsDict[count][1]:
					tempTuple = (i, difference_threshold)
					resultsDict[count] = tempTuple
			else:
					resultsDict[count] = (i, difference_threshold)
		else:
			#print("difference threshold: ", difference_threshold)
			#print("resultsdict: ", resultsDict[count], resultsDict[count][0])
			if count not in list(resultsDict.keys()):
				resultsDict[count] = (i, 6)

		#findMutations(query, matches[0], difference_threshold-1, 0, data)
		print("count: ", count)
		count+=1

print(resultsDict)
returnStr = ''
for key, value in resultsDict.items():
	if len(returnStr)==0:
		#genomeNumber = value[0]
		returnStr = '>read_'+str(key)+' Genome_Number_'+str(value[0])
	else:
		returnStr = returnStr + '\n'+'>read_'+str(key)+' Genome_Number_'+str(value[0])

with open('predictions.csv', 'w') as f:
	f.write(returnStr)



