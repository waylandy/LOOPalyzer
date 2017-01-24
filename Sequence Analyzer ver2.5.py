'''
####################################################################
#
#  LOOPER SEQUENCE ANALYZER
#  created by: Wayland Yeung
#
#  Shoutout to the Hili Lab, best lab at UGA!
#
#  This program will analyze sequencing data from performing LOOPER
#  polymerization. Input a merged FASTQ sequencing file (we used PEAR
#  to merge) and let the program do all the work! I made sure to do
#  plenty of commenting. Easy to read code makes life less stressful!
#
####################################################################
'''

import time

file = input('Input name of FASTQ:')

pr1 = 'GGATCCGAGCTCCACGTG'          # 5' primer
pr2 = 'TGCGACGGCAGGCGAATC'          # 3' primer

s = ','                             # seperator (delimiter)
e = '-'                             # error (sequence did not fulfill condition)
p = 'NOTICE MEEEEEE'                # easily spotted test phrase

start = time.time()                 # start the process timer

###########################################################

'''
# stick raw sequences from the FASTQ into the main pipeline
#
# contains nested functions for each step of the analysis
# each nested function requires specific inputs/outputs
# the specifics are listed at the bottom of the pipeline
'''

def pipeline(raw):
    # input raw sequence; output reading region
    # output error if primer regions not present
    # output error if reading region not length 40
    def noprimer(dna):
        if (pr1 and pr2) in dna:
            read = dna[dna.find(pr1)+len(pr1):dna.find(pr2)]
            if len(read) == 40:
                return read
            else:
                return 'ERROR: length ' + str(len(read))
        elif (pr1 in dna) and (pr2 not in dna):
            return 'ERROR: missing pr2'
        elif (pr2 in dna) and (pr1 not in dna):
            return 'ERROR: missing pr1'
        elif (pr1 and pr2) not in dna:
            return 'ERROR: missing both primers'
        else:
            return 'ERROR: unexpected'                      # should NEVER happen
    # input reading region; output 8 codons (5-mers)
    def grouper(dna):
        if 'ERROR' not in dna:
            n = 5 
            cod = [dna[i:i+n] for i in range(0, len(dna), n)]
            f = ''
            for num in range(8):
                f += cod[num] +s
            return f.rstrip(s)
        else:
            return e+s+e+s+e+s+e+s+e+s+e+s+e+s+e
    # input reading region; output codon fidelity
    # if all codons start with A, returns NOICE
    # if any codon not start with A, returns first mutation occurance
    # should checked in direction of polymerization (ensure this is correct)
    def mistakes(dna):
        if e not in dna:
            c = dna.split(s)
            f = ''
            for num in range(8):
                f += c[num][0]
            n = 0
            for letter in f:
                n +=1
                if letter is not 'A':
                    return 'POS' + str(n)
                elif (n == len(f)) and letter is 'A':
                    return 'GOOD'
        else:
            return e
    # input codons, reading region, codon fidelity; output rGroups
    # if passed the codon fidelity check, returns rGroups
    def rgroup(dna,read,check):
        if check is 'GOOD':
            rg = ''
            for num in range(8):
                rg += read[num*5+3:num*5+5]+s
            return rg.rstrip(s)
        else:
            return e+s+e+s+e+s+e+s+e+s+e+s+e+s+e

    noprimer = noprimer(raw)
    grouper = grouper(noprimer)
    mistakes = mistakes(grouper)
    rgroup = rgroup(grouper,noprimer,mistakes)
    
    return (noprimer+s+grouper+s+mistakes+s+rgroup)      # add .split(s) if desired
    
'''
#   If you want the pipeline to return specific attributes
#   instead of everything, simply add .split(s) to the return
#   statement above. The list below shows the individual
#   attributes available.

#   attr[1]  = reading region or error message if not available
#
#   attr[2]  = codon 1          attr[11] = rGroup sequence 1
#   attr[3]  = codon 2          attr[12] = rGroup sequence 2
#   attr[4]  = codon 3          attr[13] = rGroup sequence 3
#   attr[5]  = codon 4          attr[14] = rGroup sequence 4
#   attr[6]  = codon 5          attr[15] = rGroup sequence 5
#   attr[7]  = codon 6          attr[16] = rGroup sequence 6
#   attr[8]  = codon 7          attr[17] = rGroup sequence 7
#   attr[9]  = codon 8          attr[18] = rGroup sequence 8
#
#   attr[10] = LOOPER mutation checker
'''

# count number of lines in the input file
with open(file) as f:
    for i, l in enumerate(f):
        pass
    c = i+1

'''
#   The loop below iterates through the FASTQ file in chunks of four
#
#   label   = starts with '@'
#             followed by sequence identifier & optional description
#   seq     = raw nucleotide sequence
#   labely  = starts with '+'
#             optionally followed by sequence identifier & description
#   quality = quality score for raw sequence
#             must contain same number of characters as sequence
'''

w = open(file + '_ATTR', 'w')       # create new file in write mode
n = 0

for l in open(file):
    n += 1
    l = l.rstrip('\n')
    if (n+4)%4 == 1:
        label = l
    elif (n+4)%4 == 2:
        seq = l
    elif (n+4)%4 == 3:
        labely = l
    elif (n+4)%4 == 0:
        quality = l
        w.write(seq+s+pipeline(seq) + '\n')

w.close()
timer = time.time() - start         # end the process timer

print('\n' + '\n' + str(int(c/4)) + ' sequences processed')
print('Elapsed time: ' + str(round(timer, 2)) + ' seconds')
print('Output saved as: ' + file + '_ATTR')
