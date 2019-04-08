# STAT 6841 - Assignment 4

*Michal Lisicki*



## (1.)

Suppose that we wish to calculate $P(X<1, Y<1)​$, where $(X,Y)​$ are bivariate normally distributed with $corr(X,Y)=0.25​$. Use Monte Carlo integration to approximate this probability.

**Solution**

$(X,Y)\sim \text{BVN}(\mathbf{0}, \boldsymbol\Sigma)$

$\mu_X = 0$

$\mu_Y = 0$$

$\sigma_X^2 = 1$

$\sigma_Y^2 = 1​$

$\sigma_{XY} = \rho\sigma_X \sigma_Y =0.25​$

$$
\boldsymbol\Sigma = 
\begin{vmatrix}
\sigma_X^2 & \rho \sigma_X \sigma_Y \\
\rho \sigma_X \sigma_Y & \sigma_Y^2 \\
\end{vmatrix}
=
\begin{vmatrix}
1 & 0.5 \\
0.5 & 1 
\end{vmatrix}
= Cov(X,Y)
$$


$$
corr(X,Y) = \rho = \frac{Cov(X,Y)}{\sigma_X \ \sigma_Y} = 0.25
$$

$$
\begin{align}
f_\text{MVN}(\mathbf{z};\boldsymbol\mu,\boldsymbol\Sigma) &=(2\pi)^{-\frac{k}{2}} \det (\boldsymbol\Sigma)^{-\frac{1}{2}} e^{-\frac{1}{2}(\mathbf{z}-\boldsymbol\mu)^T\boldsymbol\Sigma(\mathbf{z}-\boldsymbol\mu)} \\
f_\text{BVN}(x,y;\boldsymbol\mu,\boldsymbol\Sigma) &= \frac{1}{2\pi \sigma_X \sigma_Y \sqrt{1-\rho^2}} \exp\bigg\{-\frac{1}{2(1-\rho^2)}\bigg[\frac{(x-\mu_X)^2}{\sigma_X^2} + \frac{(y-\mu_Y)^2}{\sigma_Y^2} - \frac{2\rho(x-\mu_X)(y-\mu_Y)}{\sigma_X\sigma_Y} \bigg]\bigg\} \\
&= \frac{1}{2\pi \sqrt{1-0.0625}} \exp\bigg\{-\frac{1}{2(1-0.0625)} \bigg[x^2 + y^2 - 2*0.25xy \bigg] \bigg\} \\
&= \frac{1}{2\pi\sqrt{0.9375}} \exp\bigg\{-\frac{1}{2(0.9375)} \bigg[x^2+y^2-0.5xy\bigg] \bigg\}
\end{align}
$$


$$
\begin{align}
P(X<1,Y<1) &= \int_{-\infty}^1 \int_{-\infty}^1 f(x,y;\boldsymbol\Sigma) dx dy \\
&= \int\int_A f(x,y;\boldsymbol\Sigma) dx dy \\
&= \int\int I_A f(x,y;\boldsymbol\Sigma) dx dy \\
&= E[I_A(x,y)]
\end{align}
$$

where

$$
I_A(x,y) =
\begin{cases}
1 & \text{if } (x,y) \in A \\
0 & \text{o.w.}
\end{cases}.
$$

Then

$$
P(X<1,Y<1) = \frac{1}{n} \sum_{i=1}^n I_A (x_i,y_i)
$$

where $X_i,Y_i \sim \text{BVN}(\boldsymbol0,\boldsymbol\Sigma)$.

To sample from BVN we need to iterate over two steps:

```python
def bvn(mu_x = 0, sigma_x = 1, mu_y = 0, sigma_y =1, rho = 0.25, n=100):
    samples = []
    for _ in range(n):
        # Step 1: generate x
        _x = random.normalvariate(mu_x,sigma_x)
        
        # Step 2: generate y
        _mu_y = mu_y + rho*(sigma_y/sigma_x)*(_x-mu_x)
        _sigma_y = sigma_y**2 * (1-rho**2)
        _y = random.normalvariate(_mu_y,_sigma_y)
        
        samples += [(_x,_y)]
    return samples
```

```python
z = bvn(n=1000)
x,y = zip(*z)
```

```python
import matplotlib.pyplot as plt
%matplotlib inline

plt.scatter(x,y)
```

![monte_carlo1](assets/monte_carlo1.png)

To make the code more efficient we could also pre-calculate $\mu_Y^* = 0.25 x​$ and $\sigma_Y^{*2} = 0.9375​$.





## (2.)

Suppose $X \sim\text{Cauchy}(2,1)​$, with density:

$$
f(x) = \frac{1}{\pi (1+(x-2)^2)},\quad -\infty < x < \infty
$$

Use importance sampling method to estimate $P(X>3)​$.

**Definitions**

Suppose we are interested in obtaining an estimate of the parameters of some function of $\phi$ a distribution $f(x)$. Assuming $f(x)$, which in Bayesian statistics can be thought of as a product of a likelihood $L(\theta)$ and a prior $\pi(\theta)$, is hard to sample from, we can use importance sampling to obtain an estimate by sampling from an easier distribution. The following notation applies:

$\theta = E_f[\phi(X)] = E_g[\chi(X)]​$

$w(x) = \frac{f(x)}{g(x)}$

$\chi(x) = \phi(x) w(x)$

The variance of the estimate is minimized when

$$
g\equiv g_0 = \frac{|\phi f|}{\int |\phi f| dx}
$$

We should choose $g(x)$ to be similar in shape to $f(x)$ but more heavy tailed.

**Solution**

In our case the function $\phi$ we want to calculate our estimate over is implicit in the constraint $P(X>3)$.

$\theta = P(X>3) = \int_3^\infty f(x) dx = \int_{-\infty}^\infty I_A(x) f(x) dx = E_f[I_{(x>3)}]$

where we set $\phi(x) = I_A(x)$ and $A = \{x: x>3\}$.

We need a function $g(x)​$ that is easier to sample from and similar in behavior to $f(x)​$. To keep things simple, we can try $g(x) = \frac{c}{(x-2)^2}​$. The constant $c​$ is a scalar, which ensures that $g(x)​$ represents a proper density, i.e. $\int_{-\infty}^\infty g(x) dx =1​$. In our case as $\int_3^\infty \frac{1}{(x-2)^2} dx = -\frac{1}{x-2} \bigg|_3^\infty​$ already equals 1 we can set $c=1​$.

Now, instead of using Monte Carlo integration to approximate $\theta$ as

$$
\hat{\theta}_f = E_f[\phi(x)] = \frac{1}{n} \sum_{i=1}^n \phi(x_i),\ x_i\sim f(x)
$$

we can use the importance sampling estimate

$$
\begin{align}
\hat{\theta}_g = E_g[\chi(x)] &= \frac{1}{n} \sum_{i=1}^n \chi (x_i) \\
&= \frac{1}{n} \sum_{i=1}^n \phi(x_i) w(x_i) \\
&= \frac{1}{n} \sum_{i=1}^n I_A(x_i) w(x_i) \\
&= \frac{1}{n} \sum_{i=1}^n  w(x_i),\ x_i\sim g(x),\ x>3
\end{align}
$$

In order to sample from any distribution numerically we can sample uniformly from the output domain of its cdf: $u_i\sim \text{Unif}(0,1)​$, and then calculate its inverse $G^{-1}(u)​$ to obtain $x_i​$.

$G(x) =P_g(X<x) = \int_3^x g(t)dt = \int_3^x \frac{1}{(t-2)^2} dt = -\frac{1}{x-2}+1 = 1-\frac{1}{x-2},\ x>3​$

To make the calculation simpler we can also sample from the inverse of the complementary cdf $\overline{G}^{-1}(u)​$ as $U_i​$ is uniformly distributed, so it doesn't if we flip the domain. The resulting distribution of $X​$ will be still the same. Therefore we have:

$\overline{G}(x) = 1-G(x) = 1- (1-\frac{1}{x-2}) = \frac{1}{x-2}​$

and the inverse can be calculated as:

$u_i = \frac{1}{x_i-2}$

$x_i=\frac{1}{u_i}+2 = \overline{G}^{-1}(u) ​$

This allows us to generate samples from $g(x)$ computationally

```python
import numpy as np
x = 1/np.random.uniform(0,1,1000)+2
```

```python
import seaborn as sns
ax = sns.distplot(x[x<10])
ax.set_xlabel('$x$')
ax.set_ylabel('$g(x)$')
```


![monte_carlo2](assets/monte_carlo2.png)

 and estimate the parameters by plugging them into:
 
$$
\begin{align}
\hat{\theta}_g &= \frac{1}{n} \sum_{i=1}^n  w(x_i) \\
&= \frac{1}{n} \sum_{i=1}^n \frac{f(x_i)}{g(x_i)} \\
&= \frac{1}{n} \sum_{i=1}^n \frac{(x_i-2)^2}{\pi (1+(x_i-2)^2)} \\
\end{align}
$$

```python
w = (x-2)**2/(np.pi+1+(x-2)**2)
```

```python
plt.hist(theta_star,50)
plt.xlabel('weight')
plt.ylabel('count')
```

![monte_carlo3](assets/monte_carlo3.png)



```python
theta = np.mean(w)
print("Theta = {}.".format(theta))
```

​    Theta = 0.5434807864325402.

which is the estimator of the expected value of $P(X>3)$.





**Additional notes**

$$
\begin{align}
\theta_f = E_f[\phi(x)] &= \int_{-\infty}^\infty \phi(x) f(x) dx \\
&= \int_{-\infty}^\infty \phi(x) \frac{f(x)}{g(x)} g(x) dx \\
&= \int_{-\infty}^\infty \phi(x) w(x) g(x) dx \\
&= \int_{-\infty}^\infty I_A(x) w(x) g(x) dx \\
&= \int_{3}^\infty w(x) g(x) dx \\
\end{align}
$$





## (3.) 

Suppose we are interested in sampling from a normal distribution $N(\mu,1)​$ where $\mu​$ is generated from $N(0,1)​$. We consider the following proposal distribution

$$
q(x'|x)\sim N(x,\tau^2)
$$

to generate a Metropolis-Hastings MCMC to estimate $\mu$. Try different values of the proposal variance $\tau^2$, e.g. $\tau^2=0.1,1,10$. From MCMC chains, report the acceptance probability, the estimate of $\mu$ and plot the trace plot for each choice of $\tau^2$. 

**Definitions**

*Proposal distribution* 

$q(x'|x)\sim N(x,\tau^2)$

