from __future__ import division
import time


####################################################################
#
#  LOOPalyzer: LOOPER SEQUENCE ANALYZER
#  created by: Wayland Yeung
#
#  Shoutout to the Hili Lab, best lab at UGA!
#
#
#  This program will analyze sequencing data from performing LOOPER
#  polymerization. Input a merged FASTQ sequencing file (we used PEAR
#  to merge) and let the program do all the work! I made sure to do
#  plenty of commenting. Easy to read code makes the world go round.
#
#  This version of LOOPalyzer is for Python 2.
#
####################################################################



# NOTE: when inputing positions, start counting at 1 NOT 0

pr1 = 'GGATCCGAGCTCCACGTG'          # 5' primer sequence
pr2 = 'TGCGACGGCAGGCGAATC'          # 3' primer sequence
lread = 40                          # length of reading region (assumed to be directly in between the primers)
kmer = 5                            # length of each codon
kper = int(lread/kmer)              # codons per reading region (DO NOT CHANGE)
modnt = 'A'                         # identity of the modified nucleotide (for error checking)
modpos = 1                          # position of modified nucleotide on codon (for error checking)
idstart = 4                         # first position of codon ID sequence
idend = 5                           # final position of codon ID sequence

s = ','                             # seperator (delimiter)
e = '-'                             # error (sequence did not fulfill condition)
p = 'NOTICE MEEEEEE'                # easily spotted test phrase

####################################################################
####################################################################

# stick raw sequences from the FASTQ into the main pipeline
#
# contains nested functions for each step of the analysis
# each nested function requires specific inputs/outputs
# the specifics are listed at the bottom of the pipeline

def pipeline(raw):
    # input raw sequence; output reading region
    # output error if primer regions not present
    # output error if reading region not as specified (default: 40)
    def noprimer(dna):
        if (pr1 in dna) and (pr2 in dna):
            read = dna[dna.find(pr1)+len(pr1):dna.find(pr2)]
            if len(read) == lread:
                return read
            else:
                return 'ERROR: length ' + str(len(read))
        elif (pr1 in dna) and (pr2 not in dna):
            return 'ERROR: missing pr2'
        elif (pr2 in dna) and (pr1 not in dna):
            return 'ERROR: missing pr1'
        elif (pr1 and pr2) not in dna:
            return 'ERROR: missing both primers'
    # input reading region; output codons (default: 8 pentamers)
    def grouper(dna):
        if 'ERROR' not in dna:
            cod = [dna[i:i+kmer] for i in range(0, len(dna), kmer)]
            f = ''
            for num in range(kper):
                f += cod[num] +s
            return f.rstrip(s)
        else:
            return ((e+s)*int(kper)).rstrip(s)
    # input reading region; output codon fidelity
    # default: if all codons start with A, returns GOOD
    # if any codon is incorrect, returns first mutation occurance
    # should checked in direction of polymerization (ensure this is correct)
    def mistakes(dna):
        if e not in dna:
            c = dna.split(s)
            f = ''
            for num in range(kper):
                f += c[num][modpos - 1]
            n = 0
            for letter in f:
                n +=1
                if letter is not modnt:
                    return 'POS' + str(n)
                elif (n == len(f)) and letter is modnt:
                    return 'GOOD'
        else:
            return e
    # input codons, reading region, codon fidelity; output rGroups
    # if passed the codon fidelity check, returns rGroups
    def rgroup(dna,read,check):
        if check is 'GOOD':
            rg = ''
            for num in range(kper):
                rg += read[num*kmer+(idstart-1):num*kmer+(idend)]+s
            return rg.rstrip(s)
        else:
            return ((e+s)*int(kper)).rstrip(s)
    # this directs workflow and saves intermediate results for output later
    noprimer = noprimer(raw)                        # reading region
    grouper = grouper(noprimer)                     # codons
    mistakes = mistakes(grouper)                    # mutations
    rgroup = rgroup(grouper,noprimer,mistakes)      # rgroups
    
    return (noprimer+s+grouper+s+mistakes+s+rgroup)

####################################################################
####################################################################
#
#   The list below shows the individual attributes available with
#   default settings. Keep in mind that the 's' variable is the
#   delimiter used to make the csv. Change the return statement to
#   any other you want and it should work. 
#
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
#
####################################################################
####################################################################
#
#   The loop below iterates through the FASTQ file in chunks of four
#
#   label   = starts with '@'
#             followed by sequence identifier & optional description
#   seq     = raw nucleotide sequence
#   labely  = starts with '+'
#             optionally followed by sequence identifier & description
#   quality = quality score for raw sequence
#             must contain same number of characters as sequence
#
####################################################################
####################################################################

def main():
    file = raw_input('Input name of FASTQ:')
    start = time.time()
    # loop that iterates through the FASTQ file
    with open(file + '_ATTR', 'w') as w:
        n = 0
        with open(file) as f:
            for l in f:
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
    # count number of lines in file
    with open(file) as f:
        for i, l in enumerate(f):
            pass
        c = i+1
    # results summary
    print '\n' + str(int(c/4)) + ' sequences processed'
    print 'Elapsed time: ' + str(round(time.time() - start, 2)) + ' seconds'
    print 'Output saved as: ' + file + '_ATTR'
if __name__ == "__main__": main()
