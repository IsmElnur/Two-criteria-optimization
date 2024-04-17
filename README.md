# Two-criteria-optimization

The basic principles of multicriteria optimization in combination with methods of random and direct passive search are considered in relation to the problem of filtering a discrete signal using the weighted moving average method.

## Formulation of the problem

A signal is specified on the interval, where the discrete sequence of samples is: 

```math
x_k=x_{min}+\frac{k(x_{max}-x_{min})}{K},
k=0, ..., K;
```
K  is  number  of  counts.

Discrete uniform noise $σ = (σ_0,...,σ_K)$ with zero average and uniformly distributed amplitude over the interval $[-a,a]: \tilde{f}_k=f_k+σ_k; σ_k=rnd(-a,a)$ is superimposed on the signal. It's necessary to filter the signal by the weighted moving average method – geometric mean:
```math
\overline{f}_k(α)=\displaystyle\prod_{j=k-M}^{k+M} \tilde{f}_j^{α_{j+M+1-k}}
```

