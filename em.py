import random
import numpy as np
import math
import matplotlib.pyplot as plt

def poisson_pdf(x,lmbd):
    '''
    Returns the Poisson probability of each given sample
    '''
    return math.exp(-lmbd)*lmbd**x/math.factorial(x)

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
    return [poisson_pdf(x)+pi*(x>0)]

def zip_random(lmbd, pi, numsamples=1):
    '''
    Draws samples from a zero-inflated Poisson (ZIP) distribution
    '''
    return [0 if random.random()<pi else poisson_random(lmbd)[0] \
            for _ in range(numsamples)]

#sample = (zip_random(3,0.3,10000))
#bins = np.arange(0, max(sample) + 1.5) - 0.5
#plt.hist(sample,bins)
#   plt.show()


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
                         np.multiply(pi*(x>0),
                                     (1-pi)*poisson_pdf(x,lmbd)
                                    )
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

    def stop(self,x):
        '''
        Implementation of the AAC stopping criterion
        '''
        next_likelihood = np.sum(np.log(zip_pdf(x)))
        self.error = np.abs(next_likelihood - self.likelihood)
        self.likelihood = next_likelihood
        return self.error < self.rho

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
        headers = ["Iter", "Pi", "Lambda", "l", "Error"]
        print("{:>5}"*len(headers).format(headers))
        i = 1
        while not self.stop():
            E_zi = self.estep(self.pi, self.lambda)
            self.pi, self.lambda = self.mstep(E_zi)
            print("{:>5}"*len(headers).format(i,self.pi,self.lambda,self.likelihood,self.error))
            i+=1
        return self.pi, self.lambda

if __name__=="__main__":
    solver = EM()
    pi, lmbd = solver.run()


