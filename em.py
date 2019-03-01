import random
import numpy as np
import math
import matplotlib.pyplot as plt

def poisson_pdf(x,lmbd):
    '''
    Returns the Poisson probability of each given sample
    '''
    return np.array([math.exp(-lmbd)*lmbd**_x/math.factorial(_x) for _x in x])

def poisson_random(lmbd, numsamples=1):
    '''
    Draws samples from a Poisson distribution
    '''
    cdf = np.append(np.cumsum([math.exp(-lmbd)*lmbd**x/math.factorial(x) \
                    for x in range(24)]),[1])
    sample = [random.random() for _ in range(numsamples)]
    return [np.where(cdf > sample[i])[0][0] for i in range(numsamples)]

def zip_pdf(x,pi,lmbd):
    '''
    Returns the ZIP probability of each given sample
    '''
    return poisson_pdf(x,lmbd)+pi*(x>0)

def zip_random(lmbd, pi, numsamples=1):
    '''
    Draws samples from a zero-inflated Poisson (ZIP) distribution
    '''
    return np.array([0 if random.random()<pi else poisson_random(lmbd)[0] \
            for _ in range(numsamples)])

class ZIP_EM(object):
    '''
    This class implements the EM algorithm for finding the MLE of parameters of
    the zero-inflated Poisson distribution
    '''
    def estep(self, x, pi, lmbd):
        '''
        E-step
        We don't need to compute the actual Q function here, as we use the
        exact estimators for ZIP distribution. Instead the purpose of this
        function is to pre-compute all the values associated with this function
        which are needed in the next step.
        '''
        E_zi = np.divide(pi*(x>0),
                         pi*(x>0) + (1-pi)*poisson_pdf(x,lmbd)
                        )
        return E_zi

    def mstep(self, x, E_zi):
        '''
        M-step
        Use the estimators to compute the estimates for the given sample.
        '''
        pi = np.sum(E_zi)/np.size(E_zi)
        lmbd = np.sum(np.multiply(1-E_zi,x))/np.sum(1-E_zi)
        return pi, lmbd

    def stop(self, x, pi, lmbd):
        '''
        Implementation of the AAC stopping criterion
        '''
        next_likelihood = np.sum(np.log(zip_pdf(x, pi, lmbd)))
        self.error = np.abs(next_likelihood - self.likelihood)
        self.likelihood = next_likelihood
        return self.error < self.rho

    def run(self, x, init_pi=0.5, init_lmbd=None, rho=1e-4):
        '''
        Runs the EM algorithm
        :param x: A vector of samples from ZIP distribution.
        :param init_pi: Initial value for pi. If not provided the default value
        0.5 will be used
        :param init_lmbd: Initial value for Poisson lmbd parameter. If not
        provided the sample mean will be used
        :param rho: Stopping criterion tolerance.
        '''
        self.lmbd = init_lmbd if init_lmbd is not None else np.mean(x)
        self.pi = init_pi
        self.rho = rho
        self.likelihood = np.inf
        headers = ["Iter", "Pi", "lmbd", "l", "Error"]
        print(("{:>15}"*len(headers)).format(*headers))
        i = 1
        while not self.stop(x, self.pi, self.lmbd):
            E_zi = self.estep(x, self.pi, self.lmbd)
            self.pi, self.lmbd = self.mstep(x,E_zi)
            print(("{:>15}"+"{:>15.5f}"*(len(headers)-1)).format(i,self.pi,self.lmbd,self.likelihood,self.error))
            i+=1
        return self.pi, self.lmbd

if __name__=="__main__":
    solver = ZIP_EM()
    sample = zip_random(3,0.3,1000)
    pi, lmbd = solver.run(sample)

#sample = (zip_random(3,0.3,10000))
#bins = np.arange(0, max(sample) + 1.5) - 0.5
#plt.hist(sample,bins)
#   plt.show()
