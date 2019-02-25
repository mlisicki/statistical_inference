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

sample = (zip_random(3,0.3,10000))
bins = np.arange(0, max(sample) + 1.5) - 0.5
plt.hist(sample,bins)
plt.show()
