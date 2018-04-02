
# This script is to get some initial estimations for how well CRISPR-Cas9 project is going to work

# standard libraries

# nonstandard libraries
import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt

# homegrown libraries
from methods import match_probability


### MAIN ###

class Analysis:

    def __init__(self,*args,**kwargs):
        """ Initialze experiment class """
        # get default settings
        settings = _get_default_analysis_settings()
        
        # update with user specification
        for arg in args: settings.update(arg)
        settings.update(kwargs)

        # update attributes with dictionary
        self.update(settings)

    def update(self,mysettings):
        """ Update self attributes with dictionary """
        for k,v in mysettings.items():
            setattr(self, k, v)

    def process_sample(self,sample):
        
        mark_occurences = {}
        mark_count = {}
        
        for label in ('alpha','beta'): 
            
            marks = np.concatenate([s[1] for s in sample[label]])
            unique,counts = np.unique(marks, return_counts=True)
            mark_count[label] = len(marks)
            mark_occurences[label] = dict(zip(unique, counts))

        marks = np.unique([m.keys() for m in mark_occurences.values()])

        total_cells = len(sample['alpha']) # arbitrary pick
        total_marks = len(marks)

        alpha = 1
        cpw = mark_count['alpha']*mark_count['beta']/(2*total_marks)

        print 'Effective CPW: {}'.format(cpw)

        for mark in marks:

            for alpha in sample['alpha']:
                for beta in sample['beta']:

                    pair_data = _data_intersect((alpha[1],),(beta[1],),(total_marks,))
                    pair_data['alpha'] = 0
                    pair_data['cpw'] = (cpw,)
                    pair_data['label'] = ((alpha[0],),(beta[0],))

                    prob = match_probability(pair_data)
                    if alpha[0] == beta[0]:
                        print pair_data
                        print '{}: {}'.format(pair_data['label'],prob)
            


class Experiment:

    def __init__(self,*args,**kwargs):
        """ Initialze experiment class """
        # get default settings
        settings = _get_default_experiment_settings()
        
        # update with user specification
        for arg in args: settings.update(arg)
        settings.update(kwargs)

        # update attributes with dictionary
        self.update(settings)


    def update(self,mysettings):
        """ Update self attributes with dictionary """
        for k,v in mysettings.items():
            setattr(self, k, v)

    def generate_cells(self):

        if self.cell_freq_distro == 'constant':
            self.freqs = [1./self.num_clonotypes for _ in xrange(self.num_clonotypes)]
        elif self.cell_freq_distro == 'power-law':
            self.freqs =_power_law_distribution(self.num_clonotypes,self.cell_freq_max,self.cell_freq_constant)

        # cell index indices
        a_inds,b_inds = np.arange(self.num_clonotypes),np.arange(self.num_clonotypes)

        # randomize indices
        if self.randomize_indices == True:
            np.random.shuffle(a_inds)
            np.random.shuffle(b_inds)

        self.cells,self.cell_freqs = [],{} # initialize containers 
        
        # create cell frequency dictionary
        for a,b,f in zip(a_inds,b_inds,self.freqs):
            self.cells.append((a,b))
            self.cell_freqs[(a,b)] = f

    def generate_sample(self):

        sample = {
                'alpha':[],
                'beta': [],
                }

        indices = np.random.choice(self.num_clonotypes,size=(self.num_cells,),p=self.freqs)
        
        for ind in indices:
            a,b = self.cells[ind]
            # figure out how many marking guides got into each cell
            if self.moi_distribution == 'poisson':
                marks = self.a_sites,size=(np.random.poisson(self.moi),)
            elif self.moi_distribution == 'constant':
                marks = np.random.choice(self.a_sites,size=(self.moi,))
            else:
                raise TypeError('Unkown moi_distribution type!')

            # figure out how many actually landed 
            if self.edit_efficiency_distribution == 'custom':
                raise NotImplementedError('Working on it...!')    
            elif self.edit_efficiency_distribution == 'constant':
                obs_a_marks = [a for a,r in zip(marks,np.random.random(len(marks))) 
                        if r < self.a_efficiency]
                obs_b_marks = [b for b,r in zip(marks,np.random.random(len(marks))) 
                        if r < self.b_efficiency]

            sample['alpha'].append((self.cells[ind][0],np.unique(obs_a_marks)))
            sample['beta'].append((self.cells[ind][1],np.unique(obs_b_marks)))
    
        # return dictionary
        return sample
        
#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _get_default_analysis_settings():
    return {
            # sample properties
            'p':0.05
            }

#------------------------------------------------------------------------------# 

def _get_default_experiment_settings():
    return {
            # sample properties
            'num_cells':100,
            # clonotype properties
            'num_clonotypes':1000,
            'cell_freq_distro':'power-law',
            'cell_freq_constant':1,
            'cell_freq_max':0.05,
            'randomize_indices':False,
            # number of targeted sites
            'a_sites':50,
            'b_sites':50,
            # tranfection characteristics
            'moi_distribution':'constant', # constant,poisson
            'moi':5,
            # edit efficiency
            'edit_efficiency_distribution':'constant', # constant,custom
            'a_efficiency':0.5,
            'b_efficiency':0.5,
            }

#------------------------------------------------------------------------------# 

def _power_law_distribution(num_cells,max_freq,alpha):
    """ Returns power law distribution using given parameters """ 
    # Lower bound
    if max_freq <= 1./num_cells:
        print 'Max. freq too low! Returning uniform distribution...'
        return [1./num_cells for _ in xrange(num_cells)]
    
    # Upper bound
    if max_freq >= 1.:
        print 'Max. freq too high! Returning delta distribution...'
        return [1] + [0 for _ in xrange(num_cells-1)]
 
    # Find a shift
    shift = root(_get_max_freq_diff, 1.0, args=(num_cells,max_freq,alpha)).x
    
    # Find best
    return _get_freq_distribution(shift,num_cells,max_freq,alpha)

#------------------------------------------------------------------------------# 

def _get_max_freq_diff(shift,num_cells,max_freq,alpha):
    """ Function for finding diff. b/w max_freq and current distribution """
    freqs = _get_freq_distribution(shift,num_cells,max_freq,alpha)
    return max_freq - freqs[0]

#------------------------------------------------------------------------------# 

def _get_freq_distribution(shift,num_cells,max_freq,alpha):
    """ Generate a normalized power-law distribution """
    freqs = np.arange(shift,num_cells+shift) ** -alpha
    return freqs/sum(freqs)

#------------------------------------------------------------------------------# 

def _data_intersect(d1,d2,num_wells):
    """ Inputs two lists of sets of indices, formats for match_probability """
    w_ij = tuple(len(set(s1).intersection(set(s2))) for s1,s2 in zip(d1,d2))
    w_i  = tuple(len(s1) - w for s1,w in zip(d1,w_ij))
    w_j  = tuple(len(s2) - w for s2,w in zip(d2,w_ij))
    w_o  = tuple(w4 - w2 - w3 - w1 for w1,w2,w3,w4 in zip(w_ij,w_i,w_j,num_wells))
    return {
            'w_i':w_i,
            'w_j':w_j,
            'w_ij':w_ij,
            'w_o':w_o,
            'w_tot':num_wells
           }

#------------------------------------------------------------------------------# 

if __name__ == '__main__':

    """ Unit tests """

    experiment = Experiment()
    experiment.generate_cells()
    sample = experiment.generate_sample()

    analysis = Analysis()
    analysis.process_sample(sample)


