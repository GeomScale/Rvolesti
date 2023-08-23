# Volume computation and sampling

## About
The `volesti` package provides [R](https://www.r-project.org/) with functions for volume estimation and sampling. In particular, it provides an R interface for the C++ library [**volesti**](https://github.com/GeomScale/volesti).

`volesti` computes approximations of volume of polytopes given as a set of points or linear inequalities or as a Minkowski sum of segments (zonotopes). There are algorithms for volume approximation as well as algorithms for sampling, rounding and rotating polytopes. Last but not least, `volesti` provides implementations of geometric algorithms to compute the score of a portfolio given asset returns and to detect financial crises in stock markets.

##  Download and install

* The latest stable version is available from CRAN.
* The latest development version is available on Github `https://github.com/GeomScale/Rvolesti`

To use the development version of `volesti` that includes the C++ code, you need to clone the repository and fetch the submodule. Follow these steps:

* Clone the main `Rvolesti` repository.
* Fetch the submodule from the [volesti](https://github.com/GeomScale/volesti) repository:
```
git submodule update --recursive --init --remote
```
* Build package by running:
```
Rcpp::compileAttributes()
devtools::build()
```
* Install Rvolesti by running:
```
install.packages("volesti")
```

The package-dependencies are: `Rcpp`, `RcppEigen`, `BH`.

## Documentation

* [Using the R Interface](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/r_interface.md)
* [Wikipage with Tutorials and Demos](https://github.com/GeomScale/volesti/wiki)
* [Tutorial given to PyData meetup](https://vissarion.github.io/tutorials/volesti_tutorial_pydata.html)

## How to update the volesti R package?

The C++ source code is retrieved from [volesti](https://github.com/GeomScale/volesti) package and placed in `src/include` of the current repository. To update the current C++ code we have to follow two steps:

- Update the `cran_include` branch in [volesti](https://github.com/GeomScale/volesti)
    - Clone the main `volesti` repository
    - Checkout to the `cran_include` branch
    - Update your code and open a PR similar to https://github.com/GeomScale/volesti/pull/277
- Retrieve the new `include` directory using submodule
```
git submodule update --recursive --init --remote
```

*Note:* it is possible the this update will brake the R interface, thus this operation should be processed with care. 

## Credits

* [Contributors and Package History](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/credits.md)
* [List of Publications](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/publications.md)

Copyright (c) 2012-2023 Vissarion Fisikopoulos

Copyright (c) 2018-2023 Apostolos Chalkis

You may redistribute or modify the software under the GNU Lesser General Public License as published by Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.