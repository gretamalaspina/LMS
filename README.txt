LMS - Splitted Levenberg-Marquardt Method for Large-Scale Sparse Problems

--------------------------------------------------------------------------------------------------------------------------------

LMS is a modification of Levenberg-Marquardt method for the solution of large scale least squares adjustment problems, which exploits the sparsity structure of the problems to partition the variaables into almost-independent subsets and approximate the linear system that arises at each iteration of LM method with a set of independent systems. To improve the quality of the approximation, the method employes a correction strategy that involves a modification of the right hand side of the resulting systems.


This python implementation of the method relies on METIS [1] for the partitioning of the variables, and on PyPardiso [2] for the solution of the lienar systems.

The Jupyter Notebook Demo_LMS.ipynb contains a comprehensive example of how to run all the phases of the code. 

The considered testcases are available at https://cloud.pmf.uns.ac.rs/s/GaSNnns9fdJeXqD.

--------------------------------------------------------------------------------------------------------------------------------

References

[1] G. Karypis and V. Kumar, Graph Partitioning and Sparse Matrix Ordering System, University of Minnesota, 2009
[2] A. Haas, https://github.com/haasad/PyPardisoProject

--------------------------------------------------------------------------------------------------------------------------------
This work is supported by the European Union’s Horizon 2020 programme under the Marie Sklodowska-Curie Grant Agreement no. 8129