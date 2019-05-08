# network-bi-partitioning

links.shp contains the geographical information of the Melbourne network in question. The entire simulation model of Melbourne is open sourced and can be accessed via https://www.cityxlab.com/dta.html.

densities.csv contains the simulated link density data for the 6-10 morning peak period with a 15-min interval. Exclusive public transport links have zero densities throughout the simulation and should be removed from the data set.

symnmf_newton.m applies the symmetric nonnegative matrix factorization given a similarity matrix as input. This code is developed by Kuang et al. in line with their paper: Kuang, D., Ding, C., Park, H., 2012. Symmetric nonnegative matrix factorization for graph clustering, Proceedings of the 2012 SIAM International Conference on Data Mining, pp. 106-117.

pareto_select.m is the main code that implements the bi-partitioning approach. As a function, it requires several input arguments. example_data.dat provides such an exmaple input data set.
