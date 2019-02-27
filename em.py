import random
import numpy as np
import math
import matplotlib.pyplot as plt

def poisson_random(lmbd, numsamples=1):
    '''
    Draws samples from a Poisson distribution
    '''
    cdf = np.append(np.cumsum([math.exp(-lmbd)*lmbd**x/math.factorial(x) \
                    for x in range(24)]),[1])
    sample = [random.random() for _ in range(numsamples)]
    return [np.where(cdf > sample[i])[0][0] for i in range(numsamples)]

def zip_random(lmbd, pi, numsamples=1):
    '''
    Draws samples from a zero-inflated Poisson (ZIP) distribution
    '''
    return [0 if random.random()<pi else poisson_random(lmbd)[0] \
            for _ in range(numsamples)]

#sample = (zip_random(3,0.3,10000))
#bins = np.arange(0, max(sample) + 1.5) - 0.5
#plt.hist(sample,bins)
#plt.show()

def poisson_pdf(x,lmbd):
    '''
    Returns the Poisson probability of each given sample
    '''
    return math.exp(-lmbd)*lmbd**x/math.factorial(x)

class ZIP_EM(object):
    '''
    This class implements the EM algorithm for finding the MLE of parameters of
    the zero-inflated Poisson distribution
    '''
    def estep(self,pi,lmbd):
        E_zi = np.divide(pi*(x>0),
                         np.multiply(pi*(x>0),
                                     (1-pi)*poisson_pdf(x,lmbd)
                                    )
                        )
        return E_zi

    def mstep(self, E_zi):
        pi = np.sum(E_zi)/np.size(E_zi)
        lmbd = np.sum(np.multiply(1-E_zi,x))/np.sum(1-E_zi)
        return pi, lmbd

    def stop(self):
        '''
        Implementation of the AAC stopping criterion
        '''

    def run(self, sample, init_pi=0.5, init_lambda=None, rho=1e-4):
        '''
        Runs the EM algorithm
        :param sample: A vector of samples from ZIP distribution.
        :param init_pi: Initial value for pi. If not provided the default value
        0.5 will be used
        :param init_lambda: Initial value for Poisson lambda parameter. If not
        provided the sample mean will be used
        :param rho: Stopping criterion tolerance.
        '''
        while not self.stop():
            E_zi = self.estep(self.pi, self.lambda)
            self.pi, self.lambda = self.mstep(E_zi)
