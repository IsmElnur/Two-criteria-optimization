# Two-criteria-optimization

The basic principles of multicriteria optimization in combination with methods of random and direct passive search are considered in relation to the problem of filtering a discrete signal using the weighted moving average method.

## Formulation of the problem

A signal is specified on the interval, where the discrete sequence of samples is: 
```math
x_k=x_{min}+\frac{k(x_{max}-x_{min})}{K},
k=0, ..., K;
```
$K$ is number of counts.

Discrete uniform noise $σ = (σ_0,...,σ_K)$ with zero average and uniformly distributed amplitude over the interval $[-a,a]: \tilde{f}_k=f_k+σ_k; σ_k=rnd(-a,a)$ is superimposed on the signal. It's necessary to filter the signal by the weighted moving average method – geometric mean:
```math
\overline{f}_k(α)=\displaystyle\prod_{j=k-M}^{k+M} \tilde{f}_j^{α_{j+M+1-k}}
```

In the formula, $\overline{f}_k$ is the values of the filtered signal; $r=2M+1$ is averaging window size; $α=(α_1,...,α_r)$ is normalized weight coefficients, such that
```math
\displaystyle\sum_{j=1}^{r} α_j=1; α_j \geq 0
```

The set of weights α should provide optimization of the filtered signal.

When calculating the noise and proximity criteria, the Manhattan geometry is used.

Noise criterion, $\omega$
```math
\omega=\displaystyle\sum_{k=0}^{K} |\overline{f}_k-\overline{f}_{k-1}|
```

Proximity criterion, $\delta$
```math
\delta=\displaystyle\sum_{k=0}^{K} |\overline{f}_k-\overline{f}_{k}|
```

Since both criteria $\omega$ and $\delta$ are mutually contradictory, to solve the problem of selecting weights $\alpha$ it is necessary to use multi-objective optimization methods. In this work, linear convolution of criteria should be used:
$$J=\lambda\omega+(1-\lambda)\delta\to \underset{\alpha}{min}; \lambda\in[0,1]$$

It's necessary for different values of weights $\lambda_l=\frac{l}{L}$ $(l=0,...,L)$, using the random search method $J(\alpha)$, to search for a minimum $J(\alpha)$ with a given probability of falling into the vicinity of the extremum $P$ with an acceptable length of the uncertainty interval $\epsilon$. The number of tests $N$ is estimated using the formula:

```math
N=\frac{ln(1-P)}{ln(1-\frac{\epsilon}{x_{max}-x_{min}})}
```
