"""
CS 309 - Computational Biology
Lab 1
Loc Tran, Phu Nguyen, Hang Ngo
"""

from enum import unique
from re import sub
from Bio import SeqIO
import matplotlib.pyplot as plt
import random

# 2
"""
This function reads the Fasta file, coverts the content into a Python string, and returns said string

parameters filename name of the Fasta file
ptype .Fasta file
return string representing the genome sequence
rtype string
"""
def readGenome(filename):
    record = SeqIO.read(filename, 'fasta')
    return str(record.seq)

# 3
"""
This function returns the size of the genome sequence (the length of the string representing the genome sequence)

parameters seqeuence the string represents the genome sequence
ptype string
return the length of the string represents the genome sequence
rtype integer
"""
def sequenceSize(sequence):
    return len(sequence)


"""
This function returns the dictionary containing all the base compositions appearing in the sequence and their frequencies

parameters sequence the string represents the genome sequence
ptype string
return dictionary containing all the base compositions appearing in the sequence and their frequencies
rtype dictionary
"""
def baseComposition(sequence):
    baseContent = {}
    for i in sequence:
        if i not in baseContent:
            baseContent[i] = 1
        else:
            baseContent[i] += 1
    return baseContent


#4
"""
This function returns the list of values on the x axis that is used to plot a graph

parameters sequence the string represents the genome sequence
           windowSize the dimension of the window
ptype string
      integer
return the list containing values on the x axis that is used to plot a graph
rtype list
"""
def x_axis(sequence, windowSize):
    xList = []
    step = windowSize//5
    for i in range(0, len(sequence)-windowSize+1, step):
        xList.append((i+i+windowSize)/2)
    return xList


"""
This function returns the list containing all the unique base compositions in the sequence

parameters sequence the string represents the genome sequence
ptype string
return the list containing all the unique base compositions in the sequence
rtype list
"""
def compositionList(sequence):
    return list(set(sequence))


"""
This function returns a dictionary containing the base compositions as keys and their percentages based on their frequency in each window as values

parameters subDict the template dictionary 
           subSequence the substring that has the length of a window size
           windowSize the dimension of the window
ptype dictionary
      string
      integer
return dictionary containing the base compositions as keys and their percentages based on their frequency in each window as values
rtype dictionary
"""
def compositionPercentPerWindow(subDict,subSequence, windowSize):
    for i in subSequence:
        if i in subDict:
            subDict[i] += 1
    countDictionary = baseComposition(subSequence)
    for key in subDict:
        if key in countDictionary:
            subDict[key] = (countDictionary[key]/windowSize) * 100
    return subDict


"""
This function returns the dictionary that stores the initial frequency values of all the base compositions in the sequence

parameters uniqueList list containing all the unique base compositions in the sequence
ptype list
return  dictionary that stores the initial frequency values of all the base compositions in the sequence
rtype dictionary
"""
def dictionaryTemplate(uniqueList):
    result = {}
    for item in uniqueList:
        if item not in result:
            result[item] = 0
    return result


"""
This function plots a graph that shows the distribution of all the base compositions in all windows with different window sizes

parameters windowSize the dimension of the window
           sequence the string represents the genome sequence
ptype integer
      string
return none
"""
def localized(windowSize, sequence):
    step = windowSize//5
    uniqueList = compositionList(sequence)
    baseCompositionDict = {list: [] for list in uniqueList}
    subDict = dictionaryTemplate(uniqueList)

    for i in range(0, len(sequence)-windowSize+1, step):
        subSequence = sequence[i:i+windowSize]
        subDict = compositionPercentPerWindow(subDict, subSequence, windowSize)
        for key in baseCompositionDict:
            if key in subDict:
                baseCompositionDict[key].append(subDict[key])

    plt.plot(x_axis(sequence, windowSize), baseCompositionDict[uniqueList[0]], label = uniqueList[0], color = 'orange')
    plt.plot(x_axis(sequence, windowSize), baseCompositionDict[uniqueList[1]], label = uniqueList[1], color = 'red')
    plt.plot(x_axis(sequence, windowSize), baseCompositionDict[uniqueList[2]], label = uniqueList[2], color = 'green')
    plt.plot(x_axis(sequence, windowSize), baseCompositionDict[uniqueList[3]], label = uniqueList[3], color = 'blue')
    plt.xlabel('Position')
    plt.ylabel('Percentage')
    plt.legend()
    plt.show()


# 5
"""
This function draws the skew diagram of the sequence and returns the minimum skew position

parameters sequence the string represent the genome sequence
ptype string
return the minimum skew position
rtype integer
"""
def skewDiagram(sequence):
    """Draw the skew diagram (#G - #C) for dna."""

    skew = 0
    skewList = []
    for nt in sequence:
        if nt == 'C':
            skew -= 1
        elif nt == 'G':
            skew += 1
        skewList.append(skew)


    plt.plot(range(len(sequence)), skewList)
    plt.xlabel('Position')
    plt.ylabel('Skew')
    plt.show()
    return skewList.index(min(skewList)) + 1


# 6
"""
This function finds the most frequent patterns in the genome sequence

parameters sequence the string represents the genome sequence
           k the length of the pattern
ptype string
      integer
return dictionary containing the most frequent patterns in the genome sequence and their frequencies
rtype dictionary
"""
def frequentPatternsInSequence(sequence, k):
    freq = {}
    for index in range(0, len(sequence) - k + 1):
        pattern = sequence[index:index+k]
        if pattern not in freq:
            freq[pattern] = 0
        freq[pattern] += 1

    maxFreq = max(freq.values())
    mostFrequentPatterns = {}
    for pattern in freq:
        if freq[pattern] == maxFreq:
            mostFrequentPatterns[pattern] = freq[pattern]

    return mostFrequentPatterns


"""
This function finds the most frequent patterns at the minimum skew position of the genome sequence

parameters sequence the string represents the genome sequence
           k the length of the pattern
ptype string
      integer
return dictionary containing the most frequent patterns at the minimum skew position of the genome sequence and their frequencies
rtype dictionary
"""
#   """Return the set of the most frequent k-mers at the position of the minimum skew."""
def frequentPatternsAtMinimumSkew(sequence, k):
    freq = {}
    step = 500//5
    begin = 885900
    end = 893600
    for index in range(begin, end-k+1, step):
        window = sequence[index:index+500]
        freqOfWindow = frequentPatternsInSequence(window, k)
        freq.update(freqOfWindow)

    maxFreq = max(freq.values())
    mostFrequentPatterns = {}
    for pattern in freq:
        if freq[pattern] == maxFreq:
            mostFrequentPatterns[pattern] = freq[pattern]

    return mostFrequentPatterns


# 7
"""
This function calculates the hamming distance between 2 patterns

parameters pattern1 string that represents the first pattern
           pattern2 string that represents the second pattern
ptype string
      string
return the hamming distance between the 2 patterns
rtype integer
"""
def hamming_distance(pattern1, pattern2):
    distance = 0
    for i in range(len(pattern1)):
        if pattern1[i] != pattern2[i]:
            distance += 1
    return distance


"""
This function generates a set of possible patterns with mismatches of a pattern

parameters pattern string represents the pattern
           d the number of mismatches
ptype string
      integer
return a set of possible patterns with mismatches of a pattern
rtype set
"""
def neighborhood(pattern, d):
    if len(pattern) == 0:
        return set([''])
    neighbors = set()
    suffixNeighbors = neighborhood(pattern[1:], d)
    for n in suffixNeighbors:
        if hamming_distance(n, pattern[1:]) < d:
            for nt in 'ACGT':
                neighbors.add(nt+n)
        else:
            neighbors.add(pattern[0]+n)
    return neighbors


"""
This function returns the reverse complement of a given pattern

parameters pattern string represents the pattern
ptype string
return the reverse complement of a given pattern
rtype string
"""
def revereseComplement(pattern):
    reversePattern = pattern[::-1]
    complement = ""
    for i in reversePattern:
        if i == "A":
            complement += "T"
        elif i == "T":
            complement += "A"
        elif i == "G":
            complement += "C"
        else:
            complement += "G"
    return complement


"""
This function finds the most frequent patterns in the genome sequence, taking their mismatches and their reverse complement into account

parameters sequence the string represents the genome sequence
           k the length of the pattern
           d the number of mismatches
ptype string
      integer
      integer
return dictionary containing the most frequent patterns in the genome sequence, taking their mismatches and their reverse complement into account
rtype dictionary
"""
def frequentPatternWithMismatches(sequence, k, d):
    freq = {}
    for index in range(0, len(sequence)-k+1, 1):
        pattern = sequence[index:index+k]
        neighbors = neighborhood(pattern, d)
        neighbors.add(revereseComplement(pattern))
        for neighbor in neighbors:
            if neighbor in freq:
                freq[neighbor] += 1
            else:
                freq[neighbor] = 1

    maxFreq = max(freq.values())
    mostFrequentPatterns = {}
    for pattern in freq:
        if freq[pattern] == maxFreq:
            mostFrequentPatterns[pattern] = freq[pattern]

    return mostFrequentPatterns


"""
This function finds the most frequent patterns at the minimum skew position of the genome sequence, taking their mismatches and their reverse complement into account

parameters sequence the string represents the genome sequence
           k the length of the pattern
           d the number of mismatches
ptype string
      integer
      integer
return dictionary containing the most frequent patterns at the minimum skew position of the genome sequence, taking their mismatches and their reverse complement into account
rtype dictionary
"""
def frequentPatternWithMismatchesAtMinimumSkew(sequence, k, d):
    freq = {}
    begin = 885900
    end = 893600
    for index in range(begin, end-k+1, 1):
        pattern = sequence[index:index+k]
        neighbors = neighborhood(pattern, d)
        neighbors.add(revereseComplement(pattern))
        for neighbor in neighbors:
            if neighbor in freq:
                freq[neighbor] += 1
            else:
                freq[neighbor] = 1

    maxFreq = max(freq.values())
    mostFrequentPatterns = {}
    for pattern in freq:
        if freq[pattern] == maxFreq:
            mostFrequentPatterns[pattern] = freq[pattern]

    return mostFrequentPatterns

def main():
    sequence = readGenome('GCA_013267415.1.fasta')

    length = sequenceSize(sequence)
    print(length)
    print(baseComposition(sequence))

    localized(20, sequence)
    localized(90, sequence)
    localized(20000, sequence)

    minimumPos = skewDiagram(sequence)
    print(minimumPos)

    print(frequentPatternsInSequence(sequence, 9))
    print(frequentPatternsAtMinimumSkew(sequence, 9))

    print(frequentPatternWithMismatches(sequence, 9, 2))
    print(frequentPatternWithMismatchesAtMinimumSkew(sequence, 9, 2))




main()
