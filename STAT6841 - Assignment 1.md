---
layout: page
title: STAT6841 Assignment 1  
---

## STAT6841 Assignment 1

*Michal Lisicki*


## (1.) 

If $$P(A) = \frac{1}{3}$$ and $$P(\bar{B}) = \frac{1}{4}$$ can $$A$$ and $$B$$ be disjoint? Explain.

**Definitions**

If $$A_1,A_2,...,A_n \in \mathcal{F}$$ are pairwise disjoint, i.e. $$A_i \cap A_j = \empty \quad  \forall i \neq j  $$   then $$P(\bigcup_\limits{i=1}^{\infty} A_i) = \sum_\limits{i=1}^\infty P(A_i)$$.

$$P(A \cap \bar{B}) = P(A)-P(A \cap B)$$

$$A \subseteq B \Rightarrow P(A) \leq P(B) $$

**Solution**

Proof by contradiction. Assuming $$A$$ and $$B$$ belong to the same event field $$\mathcal{F}$$ if $$A$$ and $$B$$ are disjoint then $$A \subseteq \bar{B}$$ and therefore $$P(A)\leq P(\bar{B})$$, but $$\frac{1}{3} \nleq \frac{1}{4}$$.

$$\blacksquare$$ 



## (2.)

Suppose that a sample space S has n elements. Find the number of subsets that can be formed from the elements of S. Justify your answer.

**Solution**

The number of all subsets equals to the number of all combinations without repetitions for all the numbers of elements we might want to draw, including 0 for empty set:
$$
\begin{align}
\sum_{k=0}^n {n\choose k} 
\end{align}
$$
which can be calculated using binomial theorem:
$$
(x+y)^n = \sum_{k=0}^n {n \choose k}x^{n-k}y^k.
$$
In our case $x=1$ and $y=1$ which gives us:
$$
\sum_{k=0}^n {n\choose k} = 2^n. \\
$$
$$\blacksquare$$



## (3.) 

Let $$X$$ be a continuous random variable with pdf $$f(x)$$ and cdf $$F(x)$$. For a fixed number $$x_0$$, define the function
$$
g(x) =
\begin{cases}
\frac{f(x)}{1-F(x_0)}, & x\geq x_0 \\
0, & x<x_0
\end{cases}
$$
Prove that $$g(x)$$ is a pdf. (Assume $$F(x_0) < 1)$$.)

**Definitions**

pdf $$f(x)$$ is a function that satisfies: $$F_X(x) = \int_{-\infty}^x f_X(t) dt, \ \forall x$$

$$\frac{dF_X(x)}{dx} = f_X(x)$$, if $$f_X(x)$$ is continuous at $$x$$

**Solution**

In other words we need to prove that:
$$
\int_{-\infty}^{\infty} g(t) dt = 1 \Rightarrow \\
\int_{x_0}^{\infty} g(t) dt = 1 \\
$$
Notice that the denominator $$1-F(x_0)$$ is a complementary cdf, so it equals $$\int_{x_0}^{\infty} f(t) dt$$ . Because this term doesn't depend on x we only need to integrate over the numerator.
$$
\int_{x_0}^{\infty} g(t) dt = \frac{\int_{x_0}^\infty f(t)dt}{\int_{x_0}^\infty f(t)dt}=1
$$
$$\blacksquare$$





## (4.)

A certain river floods every year. Suppose that the low-water mark is set at 1 and the high-water mark, denoted by $$Y$$ , has distribution function
$$
F_Y(y) = P(Y\leq y)=1-\frac{1}{y^2}, \ 1\leq y < \infty
$$

### (a) 

Verify that $$F_Y (y)$$ is a cdf.

**Definitions**

$$\lim_\limits{x\rightarrow -\infty} F_X(x) = 0$$

$$\lim_\limits{x\rightarrow \infty} F_X(x) = 1$$

$$F_X(x)$$ is non-decreasing; if $$x_1<x_2, F_X(x_1)\leq F_X(x_2)$$  

$$F_X(x)$$ is right continuous, that is $$\lim_\limits{x\rightarrow x_0^+} F_X(x) = F_X(x_0)$$ 

**Solution**

$$P(y<1)=0 \Rightarrow \lim_\limits{y\rightarrow -\infty} F_Y(y) = 0$$

$$\lim_\limits{y\rightarrow \infty} 1-\frac{1}{y^2} = 1 \Rightarrow \lim_\limits{x\rightarrow \infty} F_X(x) = 1$$

Derivative $$\frac{d F_Y (y)}{d y} = \frac{2}{y^3}$$ is positive on the range $$1\leq y \leq \infty$$ and $$0$$ when it approaches $$\infty$$, so the function is non-decreasing. 

### (b) 

Find the pdf, $$f_Y (y)$$, of $$Y$$.

**Solution**

$$f_Y(y) = \frac{d F_Y (y)}{d y} = \frac{2}{y^3}$$

### (c) 

If the low-water mark is reset at 0, and we use a unit of measurement that is $$\frac{1}{10}$$ of that given previously, the high-water mark becomes $$Z = 10(Y-1)$$. Find the $$F_Z(z)$$.

**Definitions**

Change of variables

$$Z=g(Y) \Rightarrow Y=g^{-1}(Z)=h(Z)$$

**Solution**

As we checked before, the function should be monotonically increasing on its entire domain, and therefore $$g(Y)$$ is one to one transformation. Then:
$$
h(z) = \frac{z}{10}+1 \\
F_Z(z) = F_Y(h(z)) = 1-\frac{1}{(\frac{z}{10}+1)^2} = 1 - \frac{100}{(z+10)^2}
$$



## (5.)

Let $$X_1$$ and $$X_2$$ be independent exponential ($$\lambda$$) distributed random variables.

**Definitions**

Exponential distribution is given as $$f_X(x) = \lambda e^{-\lambda x}$$,	$$0\leq x < \infty$$,	$$\lambda>0$$

Then from the independence assumption we have $$f_X(x_1,x_2) = \lambda^2 e^{-\lambda (x_1+x_2)}$$.

### (a) 

Find the joint density of $$Z_1 = \frac{X_1}{X_1+X_2}$$, $$Z_2=X_1+X_2$$.

**Solution**

We need to use the change of variables and compute the distribution $$f_{Z_1,Z_2}(z_1,z_2)$$

First we need to find out if our function is globally invertible (bijective / defines a one to one relation) or only locally (many to one) and thus has to be partitioned into regions in order to be analyzed. In a single variable case this could be achieved by checking if the function is monotone, i.e. checking if its derivative is always positive or always negative on a specified domain.

Let's calculate the Jacobian and see what happens
$$
\begin{align}
J_z(x_1,x_2) &= {\begin{bmatrix}
\frac{\partial z_1}{\partial x_1} & \frac{\partial z_1}{\partial x_2}\\
\frac{\partial z_2}{\partial x_1} & \frac{\partial z_2}{\partial x_2}\\
\end{bmatrix}} \\
&= {\begin{bmatrix}
\frac{x_2}{(x_1+x_2)^2} & -\frac{x_1}{(x_1+x_2)^2} \\
1 & 1 \\
\end{bmatrix}} \\
\end{align}
$$
Then the determinant is:
$$
|J_z|= {\begin{vmatrix}
\frac{x_2}{(x_1+x_2)^2} & -\frac{x_1}{(x_1+x_2)^2} \\
1 & 1 \\
\end{vmatrix}}=\frac{x_2}{(x_1+x_2)^2}+\frac{x_1}{(x_1+x_2)^2}=\frac{x_1+x_2}{(x_1+x_2)^2} = \frac{1}{x_1+x_2}
$$
which can be 0 when $$x_1=-x_2$$, but the function is undefined in this region, so we can ignore it. Otherwise the function is monotonically decreasing or increasing, and therefore bijective. We can distinguish two regions:
$$
A_1=\{x_1,x_2: x_1<-x_2\} \\
A_2=\{x_1,x_2: x_1>-x_2\} \\
$$
However, because both $x_1$ and $x_2$ are always positive, the first situation will never occur, and therefore the function is one to one on our region of interest.

The inverse functions are defined as follows:
$$
h_1(z_1,z_2) = z_1 z_2 \\
h_2(z_1,z_2) = z_2-z_1z_2
$$

$$
\begin{align}
J_x(z_1,z_2) &= {\begin{bmatrix}
\frac{\partial x_1}{\partial z_1} & \frac{\partial x_1}{\partial z_2}\\
\frac{\partial x_2}{\partial z_1} & \frac{\partial x_2}{\partial z_2}\\
\end{bmatrix}} \\
&= {\begin{bmatrix}
z_2 & z_1 \\
-z_2 & 1 - z_1 \\
\end{bmatrix}} \\
\end{align}
$$

Therefore the determinant equals to:
$$
|J_x|=z_2(1-z_1)+z_1z_2 = z_2 -z_1z_2+z_1z_2=z_2
$$
Then the searched distribution is:
$$
\begin{align}
f_Z(z_1,z_2) &= f_X(h_1(z_1,z_2),h_2(z_1,z_2))|J_x| \\
&= \lambda^2 z_2 e^{-\lambda(z_1z_2+z_2-z_1z_2)} \\
&= \lambda^2 z_2 e^{-\lambda z_2}
\end{align}
$$
$$\blacksquare$$

### (b) 

Show that $$Z_1$$ and $$Z_2$$ are independent.

**Definitions**

$$f_{Z_1,Z_2}(z_1,z_2) = f_{Z_1}(z_1)f_{Z_2}(z_2)$$ for independent variables

$$f_{Z_1}(z_1) = \int_\limits{-\infty}^\infty f_{Z_1,Z_2}(z_1,z_2) dz_2$$

$$\int f'(x) g(x) dx = f(x) g(x) - \int f(x) g'(x) dx$$

**Solution**
$$
\begin{align}
f_{Z_1}(z_1) &= \int_\limits{-\infty}^\infty f_{Z_1,Z_2}(z_1,z_2) dz_2 \\
&= \int_\limits{0}^\infty \lambda^2 z_2 e^{-\lambda z_2} dz_2 \\
\end{align}
$$
From integration by parts we have:
$$
\begin{align}
\int -\lambda x e^{-\lambda x} dx &= e^{-\lambda x} x - \int e^{-\lambda x} dx \\
\int \lambda x e^{-\lambda x} dx &= - e^{-\lambda x} (x+\frac{1}{\lambda}) \\
\int \lambda^2 x e^{-\lambda x} dx &= - e^{-\lambda x} (\lambda x+1)
\end{align}
$$
Coming back tour integral:
$$
\begin{align}
\int_\limits{0}^\infty \lambda^2 z_2 e^{-\lambda z_2} dz_2 &= -e^{-\lambda z_2} (\lambda z_2+1) \bigg|_0^\infty = 1 
\end{align}
$$
And the second pdf:
$$
\begin{align}
f_{Z_2}(z_2) &= \int_\limits{-\infty}^\infty f_{Z_1,Z_2}(z_1,z_2) dz_1 \\
&= \int_\limits{0}^1 \lambda^2 z_2 e^{-\lambda z_2} dz_1 \\
&= \lambda^2 z_2 e^{-\lambda z_2} z_1 \bigg|_0^1 \\
&= \lambda^2 z_2 e^{-\lambda z_2}
\end{align}
$$
Putting it all together:
$$
f_{Z_1}(z_1)f_{Z_2}(z_2) = \lambda^2 z_2 e^{-\lambda z_2} = f_{Z_1,Z_2}(z_1,z_2)
$$
$$\blacksquare$$

### (c) 

Use the following result
$$
f_{X+Y}(t) = \int_\limits{-\infty}^t f_{X,Y}(x,t-x)dx
$$
to find the pdf of $$Z_2$$.
$$
\begin{align}
f_{Z_2}(z_2) &= f_{X_1+X_2}(z_2) \\
&= \int_\limits{0}^{z_2} f_{X_1,X_2}(x,z_2-x)dx \\
&= \int_\limits{0}^{z_2} \lambda^2 e^{-\lambda (x+z_2-x)}dx \\
&= \lambda^2 e^{-\lambda z_2}x \bigg|_{0}^{z_2} \\
&= \lambda^2 z_2 e^{-\lambda z_2} \\
\end{align}
$$
$$\blacksquare$$



## (6.)

Suppose that $$X$$ is a random variable with pdf $$f_\theta(x) = \frac{\theta}{x^{\theta+1}},\ \theta>0$$ for $$x\geq 1$$. For what values of $$\theta$$ such that $$E(X)$$ and $$var(X)$$ exist.

**Solution**

*Expected value*
$$
\begin{align}
E(X) &= \int_\limits{1}^{\infty} x f_\theta(x) dx \\
&= \int_\limits{1}^{\infty} x \frac{\theta}{x^{\theta+1}} dx \\
&= \theta \int_\limits{1}^{\infty} \frac{1}{x^\theta} dx \\
&= -\frac{\theta}{(\theta-1)x^{\theta-1}}\bigg|_1^\infty \\
&= \frac{\theta}{(\theta-1)}
\end{align}
$$
This function is well-defined and $$< \infty$$ for $$\theta \neq 1$$  and $$-\infty<\theta<\infty$$.

*Variance*
$$
\begin{align}
Var(X) &= E[X^2]-(E[X])^2 \\
&=\int_\limits{1}^\infty x^2 f_\theta(x) dx - \frac{\theta^2}{(\theta-1)^2} \\
&=\int_\limits{1}^\infty x^2 \frac{\theta}{x^{\theta+1}} dx - \frac{\theta^2}{(\theta-1)^2} \\
&=\int_\limits{1}^\infty \frac{\theta}{x^{\theta-1}} dx - \frac{\theta^2}{(\theta-1)^2} \\
&= -\frac{\theta}{(\theta-2)x^{\theta-2}}\bigg|_1^\infty - \frac{\theta^2}{(\theta-1)^2} \\
&= -\frac{\theta}{(\theta-2)} - \frac{\theta^2}{(\theta-1)^2} \\
&= -\frac{\theta}{(\theta-2)(\theta-1)^2} \\
\end{align}
$$
which is well defined and  $$< \infty$$ for $$\theta \neq 1$$ and $$\theta \neq 2$$ and $$-\infty<\theta<\infty$$.

So in overall $$\theta \in \mathcal{R}\backslash \{1,2\}$$





## (7.)

 Suppose $$X_1,...,X_n$$ is a random sample from $$X$$ with pdf $$f(x)$$ and cdf $$F(x)$$. Let $$X_{(1)} \leq X_{(2)} \leq...\leq X_{(n)} \leq$$ be the order statistics. 

a) Find the pdf of $$X_{(k)}-X_{(j)}$$, where $$k \gt j$$.

**Definitions**

Joint probability of ordered statistics:

$$f_{X_{(j)}, X_{(k)}} = \frac{n!}{(j-1)!(k-j-1)!(n-k)!} [F_X(x)]^{j-1} f_X(x) [F_X(y)-F_X(x)]^{k-j-1} f_X(y) [1-F_X(y)]^{n-k}$$ 



**Solution**

We can change the variables so that:
$$
Z_1 = X_{(k)}-X_{(j)} \\
Z_2 = X_{(j)}
$$


Then the inverse functions are:
$$
X_{(j)} = Z_2 \\
X_{(k)} = Z_1+Z_2
$$

$$
f_{X_{(j)},X_{(k)}}(x,y) = f_{Z_1}(z_1) = \int_{-\infty}^\infty f_{ {Z_1},{Z_2} }(z_1,z_2) dz_2 = \int_{0}^\infty f_{X_{(k)},X_{(j)}}(z_1+z_2,z_2) dz_2 = \\
= \int_{0}^\infty \frac{n!}{(j-1)!(k-j-1)!(n-k)!} [F_X(z_2)]^{j-1} f_X(z_2) [F_X(z_1+z_2)-F_X(z_2)]^{k-j-1} f_X(z_1+z_2) [1-F_X(z_1+z_2)]^{n-k} dz_2
$$



b) Suppose $$X \sim \text{Unif}(0,1)$$, and the pdf of $$X_{(4)}-X_{(2)}$$ when the sample size $$n = 5$$.

*pdf*
$$
f(x) =
\begin{cases}
\frac{1}{1-0} = 1,& x\in[0,1] \\
0, & \text{otherwise}
\end{cases}
$$
*cdf*
$$
F(x) =
\begin{cases}
0 & \text{for } x<0 \\
\frac{x-0}{1-0}=x & \text{for } 0 \leq x \leq 1 \\
1 & \text{for } x>1
\end{cases}
$$

$$
\int_{0}^1 \frac{n!}{(j-1)!(k-j-1)!(n-k)!} [F_X(z_2)]^{j-1} f_X(z_2) [F_X(z_1+z_2)-F_X(z_2)]^{k-j-1} f_X(z_1+z_2) [1-F_X(z_1+z_2)]^{n-k} dz_2 = \\
\int_{0}^1 \frac{5!}{1!1!1!} (z_2)^1 1 (z_1+z_2-z_2)^11(1-z_1-z_2)^1 dz_2 = \\
\int_{0}^1 120 z_1 z_2 (1-z_1-z_2)dz_2 = \\
120z_1 \bigg[z_2 - \frac{z_1z_2^2}{2} - \frac{z_2^3}{3} \bigg]_0^1 = \\
120z_1 (1-\frac{1}{2}-\frac{1}{3}) = \\
20z_1 = \\
20(X_{(4)}-X_{2})
$$
