# Two-objective-optimization

The basic principles of multi-objective optimization in combination with methods of random and direct passive search are considered in relation to the problem of filtering a discrete signal using the weighted moving average method.

## Formulation of the problem

A signal is specified on the interval, where the discrete sequence of samples is: 
```math
x_k=x_{min}+\frac{k(x_{max}-x_{min})}{K},
k=0, ..., K;
```
$K$ is number of counts.

Discrete uniform noise $\alpha = (\alpha_0,...,\alpha_K)$ with zero average and uniformly distributed amplitude over the interval $[-a;a]: \tilde{f}_k=f_k+\alpha_k; \alpha_k=rnd(-a;a)$ is superimposed on the signal. It's necessary to filter the signal by the weighted moving average method – geometric mean:
```math
\overline{f}_k(α)=\displaystyle\prod_{j=k-M}^{k+M} \tilde{f}_j^{\alpha_{j+M+1-k}}
```

In the formula, $\overline{f}_k$ is the values of the filtered signal; $r=2M+1$ is averaging window size; $\alpha=(\alpha_1,...,\alpha_r)$ is normalized weight coefficients, such that
```math
\displaystyle\sum_{j=1}^{r} \alpha_j=1; \alpha_j \geq 0
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
Random values of weights are taken symmetrical with respect to the central weight $\alpha_{M+1}$ and calculated taking into account the normalization condition, id est, uniformly distributed over a sequence of residual intervals.
$$\alpha_{M+1}=rnd(0,1)$$
$$\alpha_{M}=\alpha_{M+2}=0.5rnd(0,1-\alpha_{M+1})$$
$$...$$
$$\alpha_m=\alpha_{r-m+1}=0.5rnd(0,1-\displaystyle\sum_{s=m+1}^{r-m} \alpha_s)$$
$$...$$
$$\alpha_1=\alpha_r=0.5(1-\displaystyle\sum_{s=2}^{r-1} \alpha_s)$$

The final goal of the work is to find the optimal weight $\lambda^{·}$ (by direct passive search on a grid $\lambda_l$), which minimizes the distance from the approximately found optimal value of the integral criterion $J^{·}(\omega^{·},\delta^{·})$ to the ideal point $\hat{J}(\hat{\omega},\hat{\delta})=\hat{J}(0;0)=0$:
$$dist(J^{·},\hat{J})\to\underset{\lambda}{min}$$

Formula for calculating distance:
$$dist(J^{·},\hat{J})=|\omega|+|\delta|$$

## Model problem

The original signal is:
$$f_k=\sin(x_k)+0.5$$
$$x_k=x_{min}+\frac{k(x_{max}-x_{min})}{K}; k=0,...,K; K=100$$
$$x_{min}=0; x_{max}=\pi$$

The amplitude of the uniform noise is $2a=0.5$

Sampling the convolution weight is $\lambda_l=\frac{l}{L}$ $(l=0,...,L);L=10$

The probability of falling into the vicinity of the extreme is $P=0.95$

The uncertainty interval is $\epsilon=0.01$

Sliding window size is $r=3$; $r=5$

## Progress of work

For an averaging window of size $r=3$, we calculate for each value $\lambda_l$ the minimum: $J$, $dist(J^{·},\hat{J})$, values of the weight $\alpha$, the value of the noise criterion $\omega$, and the proximity criterion $\delta$. We'll enter the measurement results into a table.

| $\lambda_l$ | $J$   | $dist(J^{·},\hat{J})$ | $\alpha$ | $\omega$ | $\delta$ |
| :---:       |:---:  |         :---:         |  :---:   | :---:    | :---:    |
| 0           |13.276 |         19.298        | 0.3 0.4 0.3 |6.022  | 13.276   |
| 0.1         |12.547 |         19.22         | 0.304 0.392 0.304 |5.939  | 13.281|
|0.2	        |11.8	  |         18.978        |	0.319 0.362 0.319	|5.638	|13.34|
|0.3	        |11.021	|18.889|	0.326 0.347 0.326	|5.504	|13.385|
|0.4	        |10.229	|18.871|	0.330 0.341 0.330	|5.466|	13.405|
|0.5	        |9.435	|18.870	|0.331 0.337 0.331|	5.452|	13.418|
|0.6	        |8.637	|18.874	|0.334 0.332 0.334|	5.435|	13.439|
|0.7	        |7.836	|18.875	|0.334 0.331 0.334|	5.434|	13.441|
|0.8	        |7.035	|18.880	|0.335 0.329 0.335|	5.431|	13.449|
|0.9	        |6.231	|18.918	|0.341 0.319 0.341|	5.424|	13.494|
|1	          |5.425	|18.914	|0.34 0.32 0.34	|5.425|	13.489|

Using the direct passive search method, we can conclude that the optimal value of the weight is $\lambda_5=0.5$. Minimum distance is $dist(J^{·},\hat{J})=18.87$

We'll build graphs of signals, where the original signal $f_k$, the filtered signal $\overline{f}_k$ and noise $\tilde{f}_k$.

![image](https://github.com/IsmElnur/Two-criteria-optimization/assets/37519575/e0fb0fe0-511b-4dc1-9acd-a51918943349)

We'll construct a graphical display of the found approximations to the optimal criteria in the coordinate system $(\omega,\delta)$ depending on the weights $\lambda_l$.

![image](https://github.com/IsmElnur/Two-criteria-optimization/assets/37519575/5fddf818-c01b-4918-af60-7d062623f7bc)

For an averaging window of size $r=5$, we calculate for each value $\lambda_l$ the minimum: $J$, $dist(J^{·},\hat{J})$, values of the weight $\alpha$, the value of the noise criterion $\omega$, and the proximity criterion $\delta$. We'll enter the measurement results into a table.

| $\lambda_l$ | $J$   | $dist(J^{·},\hat{J})$ | $\alpha$ | $\omega$ | $\delta$ |
| :---:       |:---:  |         :---:         |  :---:   | :---:    | :---:    |
|0|	10.052|	21,167|	0.004 0.472 0.047 0.472 0.004|	11.116	|10.052|
|0.1|	10.102|	19,801|	0.011 0.434 0.11 0.434 0.011	|9.648	|10.152|
|0.2|	9.974	|18,612	|0.035 0.389 0.152 0.389 0.035	|8.194	|10.418|
|0.3|	9.606	|17,346	|0.069 0.32 0.222 0.32 0.069	|6.341	|11.006|
|0.4|	9.037	|16,829	|0.113 0.079 0.216 0.079 0.113|	5.3	|11.529|
|0.5|	8.393	|16,768	|0.123 0.263 0.228 0.263 0.123|	5.046|	11.740|
|0.6|	7.646	|16,903	|0.165 0.226 0.219 0.226 0.165|	4.426|	12.477|
|0.7|	6.828	|17,041	|0.18 0.213 0.214 0.213 0.18	|4.288	|12.753|
|0.8|	5.957	|17,249	|0.194 0.196 0.221 0.196 0.194|	4.179	|13.07|
|0.9|	5.055	|17,688	|0.221 0.178 0.202 0.178 0.221|	4.108	|13.58|
|1|	4.111	|17,556	|0.214 0.183 0.207 0.178 0.221|	4.111	|13.444|

Using the direct passive search method, we can conclude that the optimal value of the weight is $\lambda_5=0.5$. Minimum distance is $dist(J^{·},\hat{J})=16.768$

We'll build graphs of signals, where the original signal $f_k$, the filtered signal $\overline{f}_k$ and noise $\tilde{f}_k$.

![image](https://github.com/IsmElnur/Two-criteria-optimization/assets/37519575/6fe69e81-cc81-4238-8c13-173468c3e605)

We'll construct a graphical display of the found approximations to the optimal criteria in the coordinate system $(\omega,\delta)$ depending on the weights $\lambda_l$.

![image](https://github.com/IsmElnur/Two-criteria-optimization/assets/37519575/be1f68ae-81c5-4e38-939e-7f92db362602)

The algorithm for solving two-objective optimization is implemented in the C++ programming language.
