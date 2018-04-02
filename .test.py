
def get_dna_code():
    return dict([
        ('A','A'),
        ('T','T'),
        ('C','C'),
        ('G','G'),
        ('B','CGT'),
        ('D','AGT'),
        ('H','ACT'),
        ('K','GT'),
        ('M','AC'),
        ('N','ACGT'),
        ('R','AG'),
        ('S','CG'),
        ('V','ACG'),
        ('W','AT'),
        ('Y','CT')])
        

def check_pam(pam,seq,dna_code = {}):

    # load DNA translation code
    dna_code = get_dna_code()

    # check length of sequences is the same
    if len(pam) != len(seq):
        raise ValueError('Lengths not equal!')
    
    # iterate through letters
    for p,s in zip(pam,seq):
        if not s in dna_code[p]:
            return False

    # if no template breaking matches
    return True


pam = 'NGG'

seqs = ['ATG','GGG','TGG']

for seq in seqs:
    print check_pam(pam,seq)


