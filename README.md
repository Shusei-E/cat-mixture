
### Finite Mixture Models with Categorical Data

See [cat-mixture.pdf](cat-mixture.pdf) for the writeup for the model and stan code.  Source is [cat-mixture.tex](cat-mixture.tex)

* It pulls the stan code from [finite-mixture_stan-05-03.stan](finite-mixture_stan-05-03.stan).
* Simulation data is created in [01-sim-data.R](01-sim-data.R)
* The stan model is fit in [03_fit-stan_on-sim.R](03_fit-stan_on-sim.R)
* The EM model is fit in [04_fit-EM.R](04_fit-EM.R) and diagnosed in [05_fit-EM.R](05_fit-EM.R).
* Simulation data and stan output is not tracked on git but some plots are kept in [figures](figures/).