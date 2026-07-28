# Citation

If Fstic contributes to work you publish, please cite the software release.

## Software

> Ruiz-Rodriguez, P. and Coscolla, M. *Fstic: allele-frequency-based genetic
> distance calculator*. Zenodo. [doi:10.5281/zenodo.16813662](https://doi.org/10.5281/zenodo.16813662)

```bibtex
@software{fstic,
  author    = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  title     = {Fstic: allele-frequency-based genetic distance calculator},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.16813662},
  url       = {https://github.com/PathoGenOmics-Lab/fstic}
}
```

The DOI above always resolves to the latest release. Zenodo also mints a
version-specific DOI for each release, which is the better one to cite if the
exact version matters for reproducibility.

## The estimators

Fstic implements published estimators. Cite the original alongside the software
for whichever you used.

| `--formula` | Reference |
| --- | --- |
| `gst`, `fst` | Nei, M. (1973) Analysis of gene diversity in subdivided populations. *PNAS* 70(12):3321-3323. |
| `jost_d` | Jost, L. (2008) GST and its relatives do not measure differentiation. *Molecular Ecology* 17(18):4015-4026. |
| `reynolds` | Reynolds, J., Weir, B.S. and Cockerham, C.C. (1983) Estimation of the coancestry coefficient. *Genetics* 105(3):767-779. |
| `nei` | Nei, M. (1972) Genetic distance between populations. *The American Naturalist* 106(949):283-292. |
| `chord` | Cavalli-Sforza, L.L. and Edwards, A.W.F. (1967) Phylogenetic analysis: models and estimation procedures. *Evolution* 21(3):550-570. |
| `rogers` | Rogers, J.S. (1972) Measures of genetic similarity and genetic distance. *Studies in Genetics VII*, University of Texas Publication 7213:145-153. |
| `bray-curtis` | Bray, J.R. and Curtis, J.T. (1957) An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs* 27(4):325-349. |

!!! note "Say which FST"

    `--formula fst` is the Nei (1973) formulation summed per locus, not the
    Weir & Cockerham (1984) theta estimator. Worth stating explicitly in a
    methods section, since readers assume the latter.

## Authors

Paula Ruiz-Rodriguez and Mireia Coscolla, at
[I2SysBio](https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/),
University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public
Health, Valencia, Spain.

Part of [PathoGenOmics Lab](https://github.com/PathoGenOmics-Lab).
