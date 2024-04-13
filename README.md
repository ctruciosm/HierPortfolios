# HierPortfolios

This second release is of this R package is already available on CRAN.

Four hierarchical portfolio allocation strategies are implemented, namely:

- Hierarchical Risk Parity (De Prado, 2016)
- Hierarchical Clustering-Based Asset Allocation (Raffinot, 2017)
- Hierarchical Equal Risk Controbution (Raffinot, 2018) 
- A Constrained Hierarchical Risk Parity Algorithm with Cluster-based Capital Allocation (Pfitzingera and Katzke, 2019)

Each strategy was implemented in an easy-to-use function: `HRP_Portfolio`, `HACC_Portfolio`, `HERC_Portfolio` and `DHRP_Portfolio`.

# References

- De Prado, M. L. (2016). Building diversified portfolios that outperform out of sample. _The Journal of Portfolio Management_, 42(4), 59-69.
- Raffinot, T. (2017). Hierarchical clustering-based asset allocation. _The Journal of Portfolio Management_, 44(2), 89-99.
- Raffinot, T. (2018). The hierarchical equal risk contribution portfolio. _Available at SSRN_ 3237540.
- Pfitzingera, J. and Katzke, N. (2019). A Constrained Hierarchical Risk Parity Algorithm with Cluster-based Capital Allocation. _Available at_  \url{https://www.ekon.sun.ac.za/wpapers/2019/wp142019/wp142019.pdf}


# Installation

To install the latest version of this package use the following commands:

install.packages("devtools")
devtools::install_github("ctruciosm/HierPorfolios")

