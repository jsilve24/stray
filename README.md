<img align="left" width="160" height="160" src="inst/stray.png" />

# stray
Multinomial Logistic-Normal Linear Regression Models

## Citation ##
The paper is currently being written. Cite how you like till then :)

## License ##
All source code freely availale under [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Installation ##
While not strictly necessary, we recommend you edit the default Makevars (on Unix systems) located at `~/.R/Makevars` to include the following line
```
CXXFLAGS = -O3 -march=native
```
this will enable full compiler optimization and vectorization. 


Whether or not this step is done, mongrel can be installed with the following:

``` r
devtools::install_github("jsilve24/stray", build_vignettes=TRUE)
```
Or to download the development version

``` r
devtools::install_github("jsilve24/stray", ref="develop", build_vignettes=TRUE)
```


## Bugs/Feature requests ##
I appreciate bug reports and feature requests. Please post to the github issue tracker [here](https://github.com/jsilve24/stray/issues). 


