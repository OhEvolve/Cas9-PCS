
"""
Build Plates
"""

# standard libraries

# nonstandard libraries
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# homegrown libraries


""" TCR Sequence declaration """

# Exons for TCRA-Constant (w/ UTR on last segment)

TCRA_exons = [
        'ATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAG',
        'AAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAG',
        'ATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGCTGAG',
        'ATCTGCAAGATTGTAAGACAGCCTGTGCTCCCTCGCTCCTTCCTCTGCATTGCCCCTCTTCTCCCTCTCCAAACAGAGGGAACTCTCCTACCCCCAAGGAGGTGAAAGCTGCTACCACCTCTGTGCCCCCCCGGCAATGCCACCAACTGGATCCTACCCGAATTTATGATTAAGATTGCTGAAGAGCTGCCAAACACTGCTGCCACCCCCTCTGTTCCCTTATTGCTGCTTGTCACTGCCTGACATTCACGGCAGAGGCAAGGCTGCTGCAGCCTCCCCTGGCTGTGCACATTCCCTCCTGCTCCCCAGAGACTGCCTCCGCCATCCCACAGATGATGGATCTTCAGTGGGTTCTCTTGGGCTCTAGGTCCTGCAGAATGTTGTGAGGGGTTTATTTTTTTTTAATAGTGTTCATAAAGAAATACATAGTATTCTTCTTCTCAAGACGTGGGGGGAAATTATCTCATTATCGAGGCCCTGCTATGCTGTGTATCTGGGCGTGTTGTATGTCCTGCTGCCGATGCCTTC']

TCRB1_exons = [
        'GACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTTCCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAG',
        'ACTGTGGCTTTACCTCGG',
        'TGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATG',
        'GTCAAGAGAAAGGATTTCTGAAGGCAGCCCTGGAAGTGGAGTTAGGAGCTTCTAACCCGTCATGGTTTCAATACACATTCTTCTTTTGCCAGCGCTTCTGAAGAGCTGCTCTCACCTCTCTGCATCCCAATAGATATCCCCCTATGTGCATGCACACCTGCACACTCACGGCTGAAATCTCCCTAACCCAGGGGGACCTTAGCATGCCTAAGTGACTAAACCAATAAAAATGTTCTGGTCTGGCCTGA']

# note: assumed UTR is same as B1...
TCRB2_exons = [
        'GAGGACCTGAAAAACGTGTTCCCACCCAAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTACCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCTGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCA',
        'GACTGTGGCTTCACCTCC',
        'GAGTCTTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCTTGCTAGGGAAGGCCACCTTGTATGCCGTGCTGGTCAGTGCCCTCGTGCTGATGGCCATG',
        'GTCAAGAGAAAGGATTCCAG']
        #'GTCAAGAGAAAGGATTCCAGAGGGCAGCCCTGGAAGTGGAGTTAGGAGCTTCTAACCCGTCATGGTTTCAATACACATTCTTCTTTTGCCAGCGCTTCTGAAGAGCTGCTCTCACCTCTCTGCATCCCAATAGATATCCCCCTATGTGCATGCACACCTGCACACTCACGGCTGAAATCTCCCTAACCCAGGGGGACCTTAGCATGCCTAAGTGACTAAACCAATAAAAATGTTCTGGTCTGGCCTGA']

# Introns for TCR-Constant (flanking on 5',3' sides)
# NOTE: these allow the use of PAMS close to the 5'/3' edges, 
#       even if these sequences are not in mRNA transcript

TCRA_introns = [('','GTAAGGGCAGCTTTGGTGCC'),
                ('ATGTCTGTTTTTCCTTTTAG','GTAAGACAGGGGTCTAGCCT'),
                ('GTGGCCTCTTGGTTTTACAG','GTGAGGGGCCTTGAAGCTGG'),
                ('TCTGTTCTTCCTCATTCCAG','ATTAAAATGATTTGGAAGAG')]

TCRB1_introns = [('','GTGAGTGGGGCCTGGGGAGA'),
                ('CTCTCCCTGCTTTCTTTCAG','GTAAGTAAGCCCTTCCTTTT'),
                ('TCCTTCCTCCGTGCCAACAG','GTAAGCAGGAGGGCAGGATG'),
                ('CAAAACTTTCTCTTCTGCAG','')]

TCRB2_introns = [('','GGTGAGTGGGGCCTGGGGAG'),
                ('TCTTCCCCTGTTTTCTTTCA','GGTAAGTGAGTCTCTCCTTT'),
                ('ATGCTCTGTTCTTGTCAACA','GTAAGGAGGAGGGTGGGATA'),
                ('AAATGTTCTCTCTTCCACAG','')]

# Compile all of these sequences into one place 

TCR_sequences = {
        'A':{
            'exons':TCRA_exons,
            'introns':TCRA_introns
            },
        'B1':{
            'exons':TCRB1_exons,
            'introns':TCRB1_introns
            },
        'B2':{
            'exons':TCRB2_exons,
            'introns':TCRB2_introns
            }
        }


""" Start main """

def main(sequences):


    # container for good oligos
    my_oligos = {}
    my_sites = {}
    my_names = {}

    labels = sequences.keys()

    # annealing overlaps
    #overlap5 = 'cttgtggaaaggacgaaacacc'.upper() # TM: 58-59 C
    #overlap3 = 'gttttagagctagaaatagcaagttaaaataaggc'.upper() # TM: 58-59 C
    
    # pcr overlaps
    overlap5 = 'TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC'.upper() # TM: 58-59 C
    overlap3 = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTC'.upper() # TM: 58-59 C

    mt_threshold = 55.0
    gRNA_minimum_spacing = 5 
    gRNA_5prime_spacing =10

    print 'Starting oligo generation...'

    for label in labels:

        print 'Starting generation for {}...'.format(label)
        my_oligos[label] = []
        my_sites[label] = []
        my_names[label] = []

        name_index = 1 

        for exon,intron in zip(sequences[label]['exons'],sequences[label]['introns']):
            
            full_seq = intron[0] + exon + intron[1]
            full_seq_rev = Seq(full_seq).reverse_complement()
            shift = len(intron[0]) # shift for indexing correct regions
            shift_rev = len(intron[1]) # shift for indexing correct regions

            # Forward read

            ind = gRNA_5prime_spacing 
            while ind < len(exon):
                if full_seq[shift+ind:shift+ind+2] == 'GG':
                    site = full_seq[ind+shift-20:ind+shift-1]
                    oligo = Seq(overlap5 + 'G' + site)
                    oligo_rev = (Seq('G') + site + overlap3).reverse_complement()
                    if site == '':
                        ind += 1
                        continue 
                    if mt.Tm_Wallace('G' + site) <= mt_threshold:
                        ind += 1
                        continue
                    
                    my_oligos[label].append(oligo)
                    my_sites[label].append(site)
                    my_names[label].append('gRNA-TCR{}C-{}F'.format(label,name_index))
                    my_oligos[label].append(oligo_rev)
                    my_sites[label].append(site)
                    my_names[label].append('gRNA-TCR{}C-{}F-RC'.format(label,name_index))
                    name_index += 1
                    ind += gRNA_minimum_spacing
                else:
                    ind += 1

            # Reverse read

            ind = gRNA_5prime_spacing 
            while ind < len(exon):
                if full_seq_rev[shift_rev+ind:shift_rev+ind+2] == 'GG':
                    site = full_seq_rev[shift_rev+ind-20:shift_rev+ind-1]
                    #oligo = overlap5 + 'G' + site + overlap3 
                    oligo = overlap5 + 'G' + site
                    oligo_rev = (Seq('G') + site + overlap3).reverse_complement()
                    if site == '':
                        ind += 1
                        continue 
                    if mt.Tm_Wallace('G' + site) <= mt_threshold:
                        ind += 1
                        continue

                    my_oligos[label].append(oligo)
                    my_sites[label].append(site)
                    my_names[label].append('gRNA-TCR{}C-{}R'.format(label,name_index))
                    my_oligos[label].append(oligo_rev)
                    my_sites[label].append(site)
                    my_names[label].append('gRNA-TCR{}C-{}F-RC'.format(label,name_index))
                    name_index += 1
                    ind += gRNA_minimum_spacing
                else:
                    ind += 1


            #raw_input()

        print 'Gene {}: {} sites'.format(label,len(my_oligos[label])/2)

    wells = ['{}{}'.format(a,b) for a in 'ABCDEFGH' for b in xrange(1,13)]


    #'''#
    for label in ['A','B1','B2']:
        print '\nPlate; Well; Name; Sequence; MT (Wallace)'
        i = 0
        for name,site,oligo in zip(my_names[label],my_sites[label],my_oligos[label]):
            #print '{}; {}; {}'.format(wells[i%96],name,oligo)
            print 'TCR{}-{}; {}; {}; {}; {}'.format(label,1+(i/96),wells[i%96],name,oligo,mt.Tm_Wallace('G'+site))
            i += 1
    #'''# 
    
    print 'HERRREE'

    B_oligos = my_oligos['B1'] + my_oligos['B2']
    print (len(B_oligos) - len(list(set(B_oligos))))/2

    #'''#
    for label in ['A','B1','B2']:
        print '\nNew plate:'
        i = 0
        for name,site,oligo in zip(my_names[label],my_sites[label],my_oligos[label]):
            if '-RC' in name:
                continue
            print '{}; {}; TM:{}'.format(name,site,mt.Tm_Wallace('G'+site))
    #'''# 

""" Namespace catch """

if __name__ == "__main__":
    main(TCR_sequences)





















