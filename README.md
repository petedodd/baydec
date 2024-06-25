# baydec
Bayesian deconvolution


# installation

```r
install.packages("devtools")
devtools::install_github("petedodd/baydec")
```

# example

Load relevant packages:

```r
library(baydec)      #this package
library(data.table)  #for handling data
library(ggplot2)     #for plotting
```


Some example data (the `darmanis`` data from BEDwARS):

```r
data(darmanis)

```


Run the MCMC deconvolution:

```r
nrep <- 1e3
smps <- deconvolve(darmanis$X,darmanis$Y,chains=4,cores=4,iter=nrep,seed=1234)

```

Summarize the outputs:
```r
smy <- postprocess.baydec(smps)

```

Plot against truth:

```r

ggplot(smy,aes(cells,value))+
  geom_point(shape=1)+
  geom_point(data=darmanis$truth,shape=1,col=2)+
  facet_wrap(~expts)

CF <- merge(darmanis$truth[,.(expts,cells,true=value)],
            smy[,.(expts,cells,estimate=value)],
            by=c('expts','cells'))
CF$cells <- as.factor(CF$cells)


ggplot(CF,aes(estimate,true,shape=cells,col=cells))+
  facet_wrap(~expts)+
  geom_abline(slope=1,col=2,intercept=0,lty=2)+
  geom_point(size=2)+
  theme(legend.position='top')

```



## License

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
