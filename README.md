
# `nimblewomble`: An R package for Bayesian wombling with `nimble`

<!-- badges: start -->
<!-- badges: end -->

The goal of `nimblewomble` is to perform Bayesian Wombling (boundary analysis) using `nimble`. For more details refer to: _Bayesian Wombling_: Sudipto Banerjee and Alan E. Gelfand<https://doi.org/10.1198/016214506000000041> and _Bayesian Modeling with Curvature Processes_: Aritra Halder, Sudipto Banerjee and Dipak K. Dey <https://doi.org/10.1080/01621459.2023.2177166>

## Installation

You can install the development version of `nimblewomble` like so:

``` r
devtools::install_github("arh926/nimblewomble")
```

## Example
![Rplot02](https://github.com/user-attachments/assets/9143789f-8e12-4373-bb46-5fa8eef622c3)

This is a basic example which shows you the workflow on a simulated data:

### Fitting a Gaussian Process

``` r
require(nimble)
require(nimblewomble)

set.seed(1)

N = 1e2
tau = 1
coords = matrix(runif(2 * N, -10, 10), ncol = 2); colnames(coords) = c("x", "y")
y = rnorm(N, 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), tau)

mc_sp = gp_fit(coords = coords, y = y, kernel = "matern1")

zbeta = zbeta_samples(y = y, coords = coords,
                      model = mc_sp$mcmc,
                      kernel = "matern1")
```
### Predicitve Inference for Rates of Change
``` r
xsplit = ysplit = seq(-10, 10, by = 1)[-c(1, 21)]
grid = as.matrix(expand.grid(xsplit, ysplit), ncol = 2)
colnames(grid) = c("x", "y")
gradients = sprates(grid = grid,
                    coords = coords,
                    model = zbeta,
                    kernel = "matern1")
                    require(ggplot2)
require(ggplot2)
require(cowplot)
require(MBA)
require(metR)

p1 = sp_ggplot(data_frame = data.frame(grid,
                                       z = gradients$estimate.sx[,"50%"],
                                       sig = gradients$estimate.sx$sig))
p1
```
![Rplot](https://github.com/user-attachments/assets/7cac91d6-ce28-4106-84d7-cba45031d1b8)


### Wombling Measure: Predicitve inference on Line Inetgrals
``` r
curve = # Pick a curve from the surface that is interesting to you
wm = spwombling(coords = coords,
                curve = curve,
                model = zbeta,
                kernel = "matern1")

col.pts = sapply(wm$estimate.wm$sig, function(x){
  if(x == 1) return("green")
  else if(x == -1) return("cyan")
  else return(NA)
})

p2 = sp_ggplot(obs, legend.key.height = 0.7, legend.key.width = 0.4, text.size = 10)

p2 + geom_path(curve, mapping = aes(x, y), linewidth = 2) + 
  geom_path(curve, mapping = aes(x, y), colour = c(col.pts, NA), linewidth = 1, na.rm = TRUE)

```
![Rplot01](https://github.com/user-attachments/assets/13153dcb-7263-4dfa-9371-7da0bd1ed587)

