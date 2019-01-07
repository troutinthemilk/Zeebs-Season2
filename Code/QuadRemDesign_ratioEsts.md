Ratio estimates of density
================
2018-08-27

#### Quad sampling: assuming perfect detection

For each transect *i* ($i = 1, \\dotsc, T$), we observe mussel counts in *n*<sub>*i*</sub> quadrats. The total number of mussels for each transect is *y*<sub>*i*</sub> (summed over all quads) and the total area sampled in a transect is *a*<sub>*i*</sub> = *n*<sub>*i*</sub> × 0.5<sup>2</sup> square-meters. The total area and number of quads sampled varies across quads, so we will use a ratio estimate of density:
$$\\hat{D} = \\dfrac{\\sum\_{i=1}^T y\_i}{\\sum\_{i=1}^T a\_i}$$

We can use the `survey` package to obtain this ratio estimate and SE for each lake. The design assumes a SRS of transects around the lake with mussel counts and area recorded fro each transect.

#### Transect sampling (no distance): assuming perfect detection

For each transect *i* ($i = 1, \\dotsc, T$), let *l*<sub>*i*</sub> be the length of the transect sampled. The total number of mussels for each transect observed by **both** divers is *y*<sub>*i*</sub>. The total area sampled in a transect is *a*<sub>*i*</sub> = *l*<sub>*i*</sub> × *w* × 2 square-meters where *w* is the (half-width) distance from the transect that divers looked for mussels. The total area and number of quads sampled varies across quads, so we will use a ratio estimate of density:
$$\\hat{D} = \\dfrac{\\sum\_{i=1}^T y\_i}{\\sum\_{i=1}^T a\_i}$$

Results for transects 1-8, assuming *w* = 0.5m:

### Try using `unmarked` for double observer removal data

This model uses multinomial distribution for observed removal counts and a Poisson model for transect population mussels counts, denoted *N*<sub>*i*</sub>. Let *y*<sub>*i**j*</sub> be the number of mussels (not clusters) removed at time *j* at transect *i*, and **y**<sub>**i**</sub> = (*y*<sub>*i*1</sub>, *y*<sub>*i*2</sub>, *y*<sub>*i*0</sub>) where *y*<sub>*i*0</sub> = *N*<sub>*i*</sub> − ∑<sub>*j*</sub>*y*<sub>*i**j*</sub> is the number of unobserved mussels. Then the model is
$$N\_i \\sim Pois(\\lambda\_i) \\ \\ \\\\
\\pmb{y\_i} \\mid N\_i \\sim Multinom(N\_i, \\pmb{\\pi})$$
 Since transects vary in area, we allow mean abundance *λ*<sub>*i*</sub> to depend on the area of the site using the model
log(*λ*<sub>*i*</sub>)=*β*<sub>0</sub> + *β*<sub>1</sub>*a*<sub>*i*</sub>
 This model is fit for Burgan and Little Birch since transect lengths varied. Area is not used in the model for Florida since all transects were 30m in length.

Estimated density is formed from transect level estimates of abundance $\\hat{N}\_i$:
$$\\hat{D} = \\dfrac{\\sum\_{i}\\hat{N}\_i}{A}$$

To look at:

-   *β*<sub>1</sub> estimates are negative?? When excluding intercept it is positive and the $\\hat{D}$ estimates are slightly lower for Burgan (LBL is about the same)
-   use negative binomial model for *N*

GOF assessment:

#### Using neg binom abundance model

GOF assessment:

#### from ddf removal sampling: assumes *w* = 0.5??

Since area is equal to length, assuming that half width is 0.5m. Use the `mrds` package to estimate density as
$$\\hat{D} = \\dfrac{n \\bar{s}}{\\hat{P}\_dA}$$
 where $\\bar{s}$ is the mean mussel count per detection, $\\hat{P}\_d$ is the estimated probability of detection on a transect (proportion of actual clusters that were detected), and *A* is the total area surveyed.

#### Design based removal estimates

Using the fisheries package `FSA`, we can get design based (model independent) estimates of abundance for each transect, independently of one another. Then total transect level estimates and divide by area surveyed to get a density estimate and SE.

#### Results

Note: detection prob rate for `multiPois` is at the mussel level while it is at the cluster level for `MRDS` models.

``` r
> kable(arrange(ests.df, Lake), digits = 3)
```

| Lake              |    Dhat|     SE| method                            |  detProb|
|:------------------|-------:|------:|:----------------------------------|--------:|
| Lake Burgan       |   0.559|  0.210| quads, ratio                      |    1.000|
| Lake Burgan       |   0.213|  0.070| Double no distance, ratio         |    1.000|
| Lake Burgan       |   0.275|  0.054| Double no distance, multiPois     |    0.952|
| Lake Burgan       |   0.221|     NA| Double no distance, MRDS          |    0.963|
| Lake Burgan       |   0.308|  0.132| Double no distance, multiNegBinom |    0.948|
| Lake Burgan       |   0.220|  0.019| Double no distance, design based  |    0.875|
| Lake Florida      |   0.071|  0.047| quads, ratio                      |    1.000|
| Lake Florida      |   0.011|  0.007| Double no distance, ratio         |    1.000|
| Lake Florida      |   0.012|  0.006| Double no distance, multiPois     |    0.750|
| Lake Florida      |   0.012|     NA| Double no distance, MRDS          |    0.960|
| Lake Florida      |   0.012|  0.008| Double no distance, multiNegBinom |    0.750|
| Lake Florida      |   0.011|  0.002| Double no distance, design based  |    0.917|
| Little Birch Lake |  24.465|  9.420| quads, ratio                      |    1.000|
| Little Birch Lake |  10.085|  2.899| Double no distance, ratio         |    1.000|
| Little Birch Lake |  11.721|  0.984| Double no distance, multiPois     |    0.643|
| Little Birch Lake |  10.851|     NA| Double no distance, MRDS          |    0.929|
| Little Birch Lake |      NA|     NA| Double no distance, multiNegBinom |       NA|
| Little Birch Lake |  10.296|  0.083| Double no distance, design based  |    0.829|
