---
layout: default
title: STAT6841 Assignment 2
description: Michal Lisicki
---


## (1.)

Let $$Y\sim \text{Uniform}[0,1]$$, and let $$X_n=Y^n$$. Prove that $$X_n\rightarrow 0$$ almost surely.

**Definitions**

- Almost sure convergence:

  $$\forall \epsilon$$ we have $$P(\lim_\limits{n\rightarrow \infty} |X_n-X|<\epsilon)=1$$, which is stronger than a regular convergence in probability defined as $$\lim_\limits{n\rightarrow \infty} P(|X_n-X|<\epsilon)=1$$. 

  We denote it as: $$X_n \xrightarrow{a.s.} X$$.

- Uniform distr:

$$
f(x)={\begin{cases}{\frac {1}{b-a}}&\mathrm {for} \ a\leq x\leq b,\\[8pt]0&\mathrm {for} \ x<a\ \mathrm {or} \ x>b\end{cases}}
$$

**Solution**

$$
f_Y(y)=
\begin{cases}
1 &\mathrm {for} \ 0\leq y\leq 1,\\
0 &\mathrm {for} \ y<0\ \mathrm {or} \ y>0
\end{cases}
$$

$$
F_Y(y) = p(Y\leq y) = \int_\limits{0}^y f_Y(y)dy = \int_\limits{0}^y 1dy
$$

Then:

$$
\begin{align}
P(\lim_\limits{n\rightarrow \infty} |Y^n-0|<\epsilon) &= P(\lim_\limits{n\rightarrow \infty} |Y^n|<\epsilon) \\
&= P(\lim_\limits{n\rightarrow \infty} Y^n<\epsilon) & \text{($$Y$$ always positive)} \\
&= P(\lim_\limits{n\rightarrow \infty} Y<\epsilon^{\frac{1}{n}}) \\
&= \lim_\limits{n\rightarrow \infty} \int_\limits{0}^{\epsilon^{\frac{1}{n}}} 1 dy \\
&= \lim_\limits{n\rightarrow \infty} y \bigg|_0^{\epsilon^{\frac{1}{n}}} \\
&= \lim_\limits{n\rightarrow \infty} \epsilon^{\frac{1}{n}} \\
&= 1
\end{align}
$$

$$\hspace{400pt}\blacksquare$$



##(2.)

Consider $$X_1,...X_n \sim Bernoulli(p),\ 0 < p < 1$$ unknown. Write the pdf in terms of
natural parameters and from that, find the MLEs of $$p$$.

**Definitions**

Bernoulli distribution

$$
f(x;p) = p^x(1-p)^{1-x}, \quad x\in\{0,1\}
$$

Exponential family

$$
f_\theta(x) = \exp\{\sum_{k=1}^K a_k(\theta) b_k(x)+c(x)+d(\theta) \}
$$

**Solution**

Bernoulli distribution can be written in terms of natural parameters, i.e. in exponential form, can be derived as:

$$
\begin{align}
f(x;p) &=\exp\{\log(p^x(1-p)^{1-x}\} \\
&= \exp\{x\log p +(1-x)\log(1-p) \} \\
&= \exp\{x(\log p-\log(1-p)) + \log(1-p) \} \\
&= \exp\{x\log \frac{p}{1-p} + \log(1-p) \}
\end{align}
$$

which is in the exponential family parametrized by $$p$$ if we assume $$a(p) = \log\frac{p}{1-p}$$ , $$b(x)=x$$, $$c(x)=0$$ and $$d(p) = \log(1+e^{a(p)})$$. The regularity conditions of MLE are satisfied for all exponential family distributions, and so we can proceed with a solution for a regular case, i.e. finding the parameters for which the score function is equal to zero:

$$
S(p) = \frac{\part}{\part p} l(p) = 0
$$

In our case the likelihood is just a product of distributions of independent samples:

$$
l(p) = \prod_{i=1}^n f(x_i;p) = \exp\{\log \frac{p}{1-p}\sum_{i=1}^n x_i + n\log(1-p) \}
$$

We can now use a property of exponential family distributions that:

$$
 \frac{\part}{\part \theta} l(\theta) = 0 \iff t_k = E(T_k)
$$

where $$t_k = \sum_\limits{i=0}^n b_k(x_i)$$ and $$T=\sum_\limits{i=0}^n b_k(X_i)$$. In our case $$K=1$$, so looking at the likelihood equation above we have $$t=\sum_\limits{i=1}^n x_i$$ and $$T = \sum_\limits{i=1}^n X_i$$. Then:

$$
\begin{align}
E(T) &= t \\
E(\sum_\limits{i=1}^n X_i) &= \sum_\limits{i=1}^n x_i \\
n E(X) &= \sum_\limits{i=1}^n x_i  & \text{(iid samples)} \\
\end{align}
$$

Because $$X  \sim Bernoulli(p)$$ , we have $$E(X)=p$$. However $$E(T) = t$$ rule applies only to exponential family, so we need to solve for $$a(p)$$  first. 

$$
\begin{align}
a(p) &= \log \frac{p}{1-p} \\
e^{a(p)} &= \frac{p}{1-p} \\
e^{a(p)}-pe^{a(p)} &= p \\
e^{a(p)} &=(1+e^{a(p)})p \\
p&=\frac{e^{a(p)}}{1+e^{a(p)}}
\end{align}
$$

The we can substitute for $$p$$ to obtain $$E(X)=\frac{e^{a(p)}}{1+e^{a(p)}}$$ and continue:

$$
\begin{align}
\frac{e^{a(p)}}{1+e^{a(p)}} &= \frac{\sum_\limits{i=1}^n x_i}{n} \\
\frac{e^{a(p)}}{1+e^{a(p)}} &= \overline{X} \\
e^{a(p)} &= \overline{X}+\overline{X}e^{a(p)} \\
e^{a(p)} (1-\overline{X}) &= \overline{X} \\
e^{a(p)} &= \frac{\overline{X}}{1-\overline{X}} \\
\hat{a}(p) &= \log \overline{X} - \log(1- \overline{X})
\end{align}
$$

This gives us the estimate  $$\hat{a}(p)$$. Now we can use the invariance principle which says that if $$\hat{\theta}$$ s the MLE of $$\theta$$, then the MLE of $$g(\theta)$$ is $$g(\hat{\theta})$$ and vice versa. By expanding $$a(p)$$ we then have:

$$
\begin{align}
\log \frac{p}{1-p} &= \log \overline{X} - \log(1- \overline{X}) \\
\frac{p}{1-p} &= \frac{\overline{X}}{1-\overline{X}} \\
p-\overline{X}p &= \overline{X} - \overline{X}p \\
\hat{p} &= \overline{X}
\end{align}
$$

$$\hspace{400pt}\blacksquare$$









## (3.)

Suppose $$X_1,...,X_n$$ is a vector of count response variables that each variable follows
a zero-inflated Poisson (ZIP) distribution with pmf as follows:

$$
P(X=x) = 
\begin{cases}
\pi +(1-\pi)e^{-\lambda},& x=0 \\
(1-\pi)\frac{e^{-\lambda}\lambda^x}{x!},& x>0
\end{cases}
$$

where $$0 < \pi < 1$$ determines the mixing proportion of excess zero that is not from the Poisson ($$\lambda$$) distribution. Formulate an EM algorithm to obtain the maximum likelihood estimates (MLEs) of $$\pi$$ and $$\lambda$$. Fully specify the observed and unobserved data, the complete-data log-likelihood function, the incomplete-data log-likelihood, the E-step (including the solution in the E-step), the M-step (including the general
solution of the MLEs), and convergence criterion.

**Definitions**

*Conditional expectation*

For example the mean of X among those times where Y=y:

$$E(X|Y=y) = \int_X x f_{X|Y}(x|y) dx$$

With respect to a function $$r(X,Y)$$:

(1) Single variable case. We need a subscript. Otherwise we wouldn't know over which density to integrate.

$$E_X(r(X,Y)) = \int_X r(X,Y)f_X(x) dx$$

(2) Conditional case

$$E(r(X,Y)|Y) = E_{X|Y}(r(X,Y)) = \int r(x,y) f_{X|Y}(x|y)dx$$

*EM algorithm*

MLE of the parameters $$\theta$$ is determined by maximizing the marginal likelihood of the observed data [Wikipedia: EM algorithm]:

$$
{\displaystyle L({\boldsymbol {\theta }};\mathbf {X} )=p(\mathbf {X} |{\boldsymbol {\theta }})=\int p(\mathbf {X} ,\mathbf {Z} |{\boldsymbol {\theta }})d\mathbf {Z} }
$$

*E-step*

$$Q(\theta|\theta^{(t)}) = E_{Z|X,\theta^{(t)}}[log L(\theta;X,Z)]$$ 

*M-step*

$$\boldsymbol\theta^{(t+1)} = \underset{\boldsymbol\theta}{\operatorname{arg\,max}} \ Q(\boldsymbol\theta|\boldsymbol\theta^{(t)})$$

We can expand it further to:

$$
\begin{align}
Q(\theta|\theta^{(t)}) &= E_{Z|X,\theta^{(t)}}[\log L(\theta;X,Z)] \\
&= E_{Z|X,\theta^{(t)}}[l(\theta;X,Z)] \\
&= E_{Z|X,\theta^{(t)}}[\log f(X,Z|\theta)] \\
&= \int_Z f(z|x,\theta^{(t)}) \log f(x,z|\theta) dz
\end{align}
$$

*Mixture model*

The mixture model has a pdf of the form:

$$
f(x; \Psi) = \sum_\limits{g=1}^G \pi_g f_g(x;\theta_g)
$$

where $$\Psi = (\pi_1,\pi_2,...,\pi_{G-1}, \theta_1,...,\theta_G)'$$

In our case $$G=2$$ and

$$\pi_1=\pi$$

$$\pi_2=1-\pi$$

$$\theta_g = \lambda$$ (we have only one parameter) 

**Solution**

*Observed (incomplete) data*

$$x_1,...,x_n$$ 

Represent the original data distributed according to zero-inflated Poisson.

*Unobserved data*

$$z_1,...,z_n$$ 

Tells us if the variable belongs to the flat or the Poisson distribution.

*Complete data*

$$\vec{y_i} = (x_i, \vec{z_i}^T)^T$$

*Incomplete data likelihood*

$$L(\pi,\lambda; \vec{x}) = \prod_\limits{i=1}^n  \pi I_{[x_i=0]} + (1-\pi) f_p(x_i;\lambda)$$

*Complete data likelihood*

$$
\begin{align}
L_c(\pi,\lambda;\vec{y}) &= \prod_{i=1}^n\prod_{g=1}^G [\pi_g f_g(x_i|\theta_g)]^{z_{ig}} \\
&=  \prod_\limits{i=1}^n (\pi I_{[x_i=0]})^{z_i}  ((1-\pi)f_p(x_i|\lambda))^{(1-z_i)} \\
\end{align}
$$

$$
\begin{align}
l_c(\pi,\lambda; \vec{y})  &=   \sum_\limits{i=1}^n \log[ (\pi I_{[x_i=0]})^{z_i} ((1-\pi)f_p(x_i|\lambda)^{(1-z_i)}] \\
&= \sum_\limits{i=1}^n z_i \log (\pi I_{[x_i=0]}) + (1-z_i) \log ((1-\pi)f_p(x_i|\lambda)) \\
&= \sum_\limits{i=1}^n z_i \log (\pi I_{[x_i=0]}) + (1-z_i) \log (1-\pi) +(1-z_i) \log f_p(x_i|\lambda) \\
\end{align}
$$

We shall leave it at this form as it allows us to nicely divide the sum into components dependent only on $$\pi$$ and only $$\lambda$$ and find their MLEs separately.

*E-step* 

$$
\begin{align}
Q(\pi,\lambda|\pi^{(0)},\lambda^{(0)}) &= E_{Z|X,\pi^{(0)},\lambda^{(0)}}\bigg[\sum_\limits{i=1}^n z_i \log (\pi I_{[x_i=0]}) + (1-z_i) \log (1-\pi) +(1-z_i) \log f_p(x_i|\lambda) \bigg] \\
&=  \sum_\limits{i=1}^n E_{z_i|x_i,\pi^{(0)},\lambda^{(0)}}[z_i] \log (\pi I_{[x_i=0]}) + (1-E_{z_i|x_i,\pi^{(0)},\lambda^{(0)}}[z_i]) \log (1-\pi) +(1-E_{z_i|x_i,\pi^{(0)},\lambda^{(0)}}[z_i]) \log f_p(x_i|\lambda) \\
&= \sum_\limits{i=1}^n z_i^{(1)} \log (\pi I_{[x_i=0]}) + (1-z_i^{(1)}) \log (1-\pi) +(1-z_i^{(1)}) \log f_p(x_i|\lambda) \\
\end{align}
$$

$$
\begin{align}
Q_1(\pi|\pi^{(0)}) &= \sum_\limits{i=1}^n z_i^{(1)} \log (\pi I_{[x_i=0]}) + (1-z_i^{(1)}) \log (1-\pi) \\
Q_2(\lambda|\lambda^{(0)}) &= \sum_\limits{i=1}^n (1-z_i^{(1)}) \log f_p(x_i|\lambda)
\end{align}
$$

As we don't know $$\pi$$ and $$\lambda$$ the only way we can further expand this equation is by computing 

$$
\begin{align}
z_i^{(1)} &= E_{z_i|x_i,\pi^{(0)},\lambda^{(0)}}[z_i] = \sum_\limits{z_i} z_i f(z_i|x_i,\lambda^{(0)},\pi^{(0)}) \\
\end{align}
$$

From Bayes rule have:

$$
\begin{align}
f(z_i|x_i,\lambda^{(0)},\pi^{(0)}) &=\frac{f(z_i,x_i|\lambda^{(0)},\pi^{(0)})}{f(x_i|\lambda^{(0)},\pi^{(0)})} \\
&= \frac{[\pi^{(0)} I_{[x_i=0]}]^{z_i} ((1-\pi^{(0)})f_p(x_i|\lambda^{(0)}))^{(1-z_i)}}{[\pi^{(0)} I_{[x_i=0]}] +((1-\pi^{(0)})f_p(x_i|\lambda^{(0)}))}
\end{align}
$$

Note that the denominator is simply what we agreed is $$f(x|\pi,\lambda)$$ at the beginning of the section.

Then because $$z_i \sim Binomial(p)$$ we have:

$$
\begin{align}
E_{z_i|x_i,\pi^{(0)},\lambda^{(0)}}[z_i] &= 0 \cdot \frac{[\pi^{(0)} I_{[x_i=0]}]^{0} ((1-\pi^{(0)})f_p(x_i|\lambda^{(0)}))^{(1-0)}}{[\pi^{(0)} I_{[x_i=0]}] +((1-\pi^{(0)})f_p(x_i|\lambda^{(0)}))}+1\cdot\frac{[\pi^{(0)} I_{[x_i=0]}]^{1} ((1-\pi^{(0)})f_p(x_i|\lambda^{(0)}))^{(1-1)}}{[\pi^{(0)} I_{[x_i=0]}] +((1-\pi^{(0)})f_p(x_i|\lambda^{(0)}))} \\
&= \frac{[\pi^{(0)} I_{[x_i=0]}]}{[\pi^{(0)} I_{[x_i=0]}] +((1-\pi^{(0)})f_p(x_i|\lambda^{(0)}))}
\end{align}
$$

*M-step* 

Now we are looking for parameters $$\pi$$ and $$\lambda$$ which maximize our expectation:

$$
\begin{align}
\frac{d Q_1}{d \pi} &= \frac{d}{d\pi} \sum_\limits{i=1}^n z_i^{(1)} \log (\pi I_{[x_i=0]}) + (1-z_i^{(1)}) \log (1-\pi) \\
&= \sum_\limits{i=1}^n \frac{z_i^{(1)}}{\pi} - \frac{(1-z_i^{(1)})}{1-\pi} \\
\end{align}
$$

$$
\begin{align}
\sum_\limits{i=1}^n \frac{z_i^{(1)}}{\pi} - \frac{(1-z_i^{(1)})}{1-\pi} &=0\\
\sum_\limits{i=1}^n \frac{(1-\pi)z_i^{(1)}-\pi(1-z_i^{(1)})}{\pi(1-\pi)}  &=0\\
\sum_\limits{i=1}^n (1-\pi)z_i^{(1)}-\pi(1-z_i^{(1)})  &=0\\
\sum_\limits{i=1}^n  z_i^{(1)}-\pi z_i^{(1)} -\pi + \pi z_i^{(1)} &=0\\
\pi^{(1)} = \hat{\pi} &= \frac{\sum_\limits{i=1}^n  z_i^{(1)}}{n} = \overline{Z}^{(1)}
\end{align}
$$


$$
\begin{align}
\frac{dQ_2}{d\lambda} &= \frac{d}{d\lambda} \sum_\limits{i=1}^n (1-z_i^{(1)}) \log \frac{e^{-\lambda}\lambda^{x_i}}{x_i!} \\
&= \frac{d}{d\lambda} \sum_\limits{i=1}^n (1-z_i^{(1)}) (-\lambda + x_i\log\lambda-\log x_i!) \\
&= \sum_\limits{i=1}^n (1-z_i^{(1)}) (\frac{x_i}{\lambda}-1)
\end{align}
$$

$$
\begin{align}
\sum_\limits{i=1}^n (1-z_i^{(1)}) (\frac{x_i}{\lambda}-1) &=0 \\
\sum_\limits{i=1}^n \frac{(1-z_i^{(1)})x_i}{\lambda} -(1-z_i^{(1)})&=0 \\
\sum_\limits{i=1}^n (1-z_i^{(1)})x_i &= \lambda \sum_\limits{i=1}^n (1-z_i^{(1)}) \\
\lambda^{(1)} = \hat{\lambda} &= \frac{\sum_\limits{i=1}^n (1-z_i^{(1)}) x_i}{\sum_\limits{i=1}^n (1-z_i^{(1)})}
\end{align}
$$

*k-th step*

A the k-th E-step we just need to substitute $$z_i^{(1)}$$ with $$z_i^{(k)}$$ given by:

$$
z_i^{(k+1)} = E_{z_i|x_i,\pi^{(0)},\lambda^{(0)}}[z_i] =  \frac{[\pi^{(k)} I_{[x_i=0]}]}{[\pi^{(k)} I_{[x_i=0]}] +((1-\pi^{(k)})f_p(x_i|\lambda^{(k)}))}
$$

and then maximize $$Q$$ by plugging in the updated values of the unobserved data:

$$
\begin{align}
\pi^{(k+1)} &= \overline{Z}^{(k+1)} \\
\lambda^{(k+1)} = \hat{\lambda} &= \frac{\sum_\limits{i=1}^n (1-z_i^{(k+1)}) x_i}{\sum_\limits{i=1}^n (1-z_i^{(k+1)})}
\end{align}
$$

Repeat E and M steps until convergence, e.g. by checking if

$$
|l(\pi^{(k+1)},\lambda^{(k+1)}; \vec{x})-l(\pi^{(k)},\lambda^{(k)}; \vec{x})| < \rho
$$

where $$\rho$$ is a specified tolerance, and also making sure that the incomplete data likelihood is not decreasing to guearantee the proper convergence:

$$
l(\pi^{(k+1)},\lambda^{(k+1)}; \vec{x}) \geq l(\pi^{(k)},\lambda^{(k)}; \vec{x})
$$


## (4.)

Refer to question 2, use R to implement your EM algorithm, and evaluate your method by performing a simulation study. Suppose that you generate 10 samples with each of size n = 100 count responses from the ZIP distribution with π = 0.3 and λ = 3. Find the MLEs of π and λ using the EM algorithm based on each of your simulated data set. Report the empirical mean, the variance, and the mean square error of your MLEs of π and λ.

**NOTE:** Change $$\lambda$$ to 1.5.

## (5.)

Three independent binomial experiments are conducted with $$n_1$$, $$n_2$$ and $$n_3$$ trails with $$x_1$$, $$x_2$$, and $$x_3$$ are the respective numbers of successes observed.

**Definitions**

Binomial distribution's pmf:

$$
f(k,n,p) = P(X=k) = {n \choose k}p^k(1-p)^{(n-k)}
$$

#### (a)

Suppose that the probability of success, $$P$$, is the same in each trial. Find the MLE of $$p$$, based on all three trials.

**Solution**

The joint distributions of the independent binomial experiments:

$$
f(\vec{x},\vec{n},p) = \prod_{i=1}^3 {n_i \choose x_i}p^{x_i}(1-p)^{(n_i-x_i)}
$$

We can convert it as follows:

$$
\begin{align}
f(\vec{x},\vec{n},p) &= \exp\{\sum_{i=1}^3 \log [ {n_i \choose x_i}p^{x_i}(1-p)^{(n_i-x_i)}]\} \\
&= \exp\{\sum_{i=1}^3 \log [ \frac{n!}{x_i!(n-x_i)!} p^{x_i}(1-p)^{(n_i-x_i)}]\} \\
&= \exp\{\sum_{i=1}^3 \log  \frac{n!}{x_i!(n-x_i)!} + x_i \log p + (n_i-x_i) \log (1-p)\} \\
&= \exp\{\sum_{i=1}^3 \log  \frac{n!}{x_i!(n-x_i)!} + x_i \log p + (n_i-x_i) \log (1-p)\} \\
&= \exp\{\sum_{i=1}^3 x_i \log \frac{p}{1-p} + \log  \frac{n!}{x_i!(n-x_i)!} + n_i \log (1-p)\} \\
\end{align}
$$

Which gives us the exponential family form where:

$$
\begin{align}
a(p) &=  \log \frac{p}{1-p} \\
b(x) &= \sum_{i=1}^3 x_i \\
c(x) &= \sum_{i=1}^3 \log  \frac{n!}{x_i!(n-x_i)!} \\
d(p) &= \log (1-p) \sum_{i=1}^3 n_i
\end{align}
$$

Now we can find MLE using the regular approach:

$$
\begin{align}
E(T) &= t \\
E(\sum_{i=1}^3 X_i) &= \sum_{i=1}^3 x_i \\
\sum_{i=1}^3 E(X_i) &= \sum_{i=1}^3 x_i \\
\end{align}
$$

$$X_i \sim \text{Binomial}(n_i,p)$$, so $$E(X_i)=n_i p$$. Therefore:

$$
\begin{align}
p \sum_{i=1}^3 n_i &= \sum_{i=1}^3 x_i \\
\hat{p} &= \frac{\sum_{i=1}^3 x_i}{\sum_{i=1}^3 n_i}
\end{align}
$$

#### (b)

Suppose the probability of success varies between trials, and is $$p$$, $$p + a$$, and $$p$$, respectively. Find MLEs of $$p$$ and $$a$$, based on all three trials.

**Solution**

In such case we'll have:

$$
\begin{align}
f(\vec{x},\vec{n},p,a) &= p^{x_1}(1-p)^{(n_1-x_1)} (p+a)^{x_2}(1-(p+a))^{(n_2-x_2)} p^{x_3}(1-p)^{(n_3-x_3)} \prod_{i=1}^3 {n_i \choose x_i} \\
l(p,a;\vec{x},\vec{n}) &= (x_1+x_3)\log p + (n_1+n_3-x_1-x_3) \log (1-p) +x_2 \log (p+a) +(n_2-x_2)\log(1-p-a) +\sum_{i=1}^3 \log  \frac{n!}{x_i!(n-x_i)!}  \\
\end{align}
$$

$$
\begin{align}
\frac{\part l}{\part p} &=  \frac{(x_1+x_3)}{p} - \frac{n_1+n_3-x_1-x_3}{1-p} + \frac{x_2}{p+a} -\frac{n_2-x_2}{1-p-a} 
\end{align}
$$


$$
\begin{align}
\frac{\part l}{\part a} = \frac{x_2}{p+a}-\frac{(n_2-x_2)}{1-p-a }
\end{align}

$$
To find MLE we need to solve the following system of equations:
$$

\begin{align}
&\begin{cases}
\frac{(x_1+x_3)}{p} - \frac{n_1+n_3-x_1-x_3}{1-p} + \frac{x_2}{p+a} -\frac{n_2-x_2}{1-p-a}  &= 0 \\
\frac{x_2}{p+a}-\frac{(n_2-x_2)}{1-p-a } &= 0 \\
\end{cases} \\\\
&\begin{cases}
\frac{(x_1+x_3)}{p} - \frac{n_1+n_3-x_1-x_3}{1-p}  &= 0 \\
\frac{x_2}{p+a}-\frac{(n_2-x_2)}{1-p-a } &= 0 \\
\end{cases} \\\\
&\begin{cases}
\frac{(x_1+x_3)}{p}   &= \frac{n_1+n_3-x_1-x_3}{1-p} \\
\frac{x_2}{p+a}-\frac{(n_2-x_2)}{1-p-a } &= 0 \\
\end{cases} \\\\
&\begin{cases}
(1-p)(x_1+x_3)   &= p (n_1+n_3-x_1-x_3) \\
\frac{x_2}{p+a}-\frac{(n_2-x_2)}{1-p-a } &= 0 \\
\end{cases} \\\\
&\begin{cases}
x_1+x_3   &= p (n_1+n_3) \\
\frac{x_2}{p+a}-\frac{(n_2-x_2)}{1-p-a } &= 0 \\
\end{cases} \\\\
&\begin{cases}
\hat{p}   &= \frac{x_1+x_3}{n_1+n_3} \\
x_2(1-p-a) &= (n_2-x_2)(p+a) \\
\end{cases} \\\\
&\begin{cases}
\hat{p}   &= \frac{x_1+x_3}{n_1+n_3} \\
x_2-x_2(p+a) &= n_2(p+a) -x_2(p+a) \\
\end{cases} \\\\
&\begin{cases}
\hat{p}   &= \frac{x_1+x_3}{n_1+n_3} \\
a &=\frac{x_2}{n_2}-p\\
\end{cases} \\\\
&\begin{cases}
\hat{p}   &= \frac{x_1+x_3}{n_1+n_3} \\
\hat{a} &=\frac{x_2}{n_2}-\frac{x_1+x_3}{n_1+n_3}\\
\end{cases} \\\\
\end{align}
$$




## (6.)

Let $$(X_1 , X_2 , . . . , X_n )$$ be a random sample from geometric distribution with pmf $$f (x) = p(1 − p)^x$$ , $$x = 0, 1, 2, ...$$ where $$0 < p < 1$$.

**Definitions**

*Cramer-Rao Lower Bound*

$$var(\hat{\tau})\geq \frac{[\tau'(\theta)]^2}{nJ(\theta)}$$

Or in other form:

$$var_{\theta}(T(X)) \geq \frac{(\frac{d}{d\theta} E_\theta T(X))^2}{n E_\theta (( \frac{\part}{\part \theta} \log f(X|\theta))^2)}$$

where $$E_\theta=E_{X|\theta}$$ and $$n$$ is for iid case

*Fisher Information*

$$J(\theta)=-E[\frac{\part^2 \log f}{\part \theta^2}]$$

*CRLB attainment condition*

If $$T(X)=T(X_1,...,X_n)$$ is any unbiased estimator of $$\tau(\theta)$$ then $$T(X)$$ attains the CRLB iff:

$$
a(\theta) [T(X)-\tau(\theta)] = \frac{\part}{\part \theta} \log L(\theta|X)
$$

for some function $$a(\theta)$$. 

#### (a)

Find the Cramer-Rao lower bound for the variance of unbiased estimate of the population mean, i.e., $$g(p) = E(X)$$.

So in our case $$\tau(p) = g(p) = E(X)$$.

$$
\begin{align}
E(X) &= \sum_{x=0}^\infty xp(1-p)^x\\
\end{align}
$$

We know that the sum of a series $$\sum_{x=0}^\infty xa^x = \frac{a}{(a-1)^2}$$. Therefore:

$$
E(X) = p\frac{1-p}{(1-p-1)^2}=p\frac{1-p}{p^2}=\frac{1-p}{p}\\
$$

$$
\\
\frac{d g(p)}{dp} = \frac{d}{dp} \frac{1-p}{p} = -\frac{1}{p^2} (1-p) - \frac{1}{p} = -\frac{1}{p^2}\\
\\
$$

$$
\begin{align}
\frac{\part^2 \log f}{\part p^2} &= \frac{\part^2}{\part p^2}\log p+ x\log(1-p) \\
&= \frac{d}{dp} \frac{1}{p} -\frac{x}{1-p} \\
&= -\frac{1}{p^2}-\frac{x}{(1-p)^2} \\
\end{align}
$$

$$
\begin{align}
E_X\bigg[\frac{\part^2 \log f}{\part p^2}\bigg] &= -E \bigg[-\frac{1}{p^2}-\frac{x}{(1-p)^2} \bigg] \\
&=\frac{1}{p^2}+\frac{E[x]}{(1-p)^2} \\
&= \frac{1}{p^2}+\frac{1}{(1-p)^2}\frac{1-p}{p} \\
&= \frac{1}{p^2}+\frac{1}{p(1-p)} \\
&= \frac{1-p+p}{p^2(1-p)} \\
&= \frac{1}{p^2(1-p)}
\end{align}
$$

$$
\begin{align}
var(\hat{\tau}) &\geq \frac{p^2(1-p)}{np^4} \\
var(\hat{\tau}) &\geq \frac{(1-p)}{np^2} \\
\end{align}
$$

#### (b)

Is there a function of $$p$$ for which there exists an unbiased estimator the variance of which coincides with the CRLB? If so, find it and compare your result with the result you got from a).

**Definitions**

The equality on CRLB holds only if:

$$
S(\theta) = \sum_\limits{i=1}^n \frac{\part}{\part \theta} \log f_\theta(x_i) = K(\theta,n) [\tau(T(X_1,...,X_n)) - \tau(\theta)]
$$

**Solution**

From the previous section our estimator $$\tau(p)=\frac{1-p}{p}$$. Then:

$$
\begin{align}
\sum_\limits{i=1}^n \frac{\part}{\part p} \log f_p(x_i) &= \sum_\limits{i=1}^n \frac{\part}{\part p} \log p+ x_i\log(1-p) \\
&= \sum_\limits{i=1}^n \frac{1}{p} -\frac{x_i}{1-p} \\
&= -\frac{1}{1-p} \bigg( \sum_\limits{i=1}^n x_i -\frac{1-p}{p} \bigg) \\
&=\frac{n}{p-1} \bigg( \frac{\sum_\limits{i=1}^n x_i}{n} -\frac{1-p}{p} \bigg) \\
&=\frac{n}{p-1} \bigg(\overline{X} -\frac{1-p}{p} \bigg)
\end{align}
$$

$$\overline{X}$$ is then the unbiased estimator which variance attains the CRLB for our $$\tau(p)$$.

