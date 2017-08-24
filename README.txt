This package contains codes implementing the algorithms introduced in 'Adjoint Map Representation for Shape Analysis and Matching' by Ruqi Huang and Maks Ovsjanikov, SGP 2017. 
We put all dependencies but CVX (which you should download and install separately) in this package, all the data sets (Triangular meshes) involved are included as well. 
This package is currently only pre-compiled for Mac OSX operating system, we will try to make it avaiable for other platforms in the future. 

To run this code, just simply run the scripts in the main folder. The '_F[num]', say, '_F1', in the end of each script name indicates which figure in the paper is re-produced by this script. 

Some small Remarks: 
1) The data we directly load to generate Figure 2 and Figure 3 is produced by the script 'Adjoint_regularization_F1.m'.  

2) Due to some mis-displacement in the early version of the codes, we unfortunately are not able to re-produce exactly Figure 8 in the paper, as we lost the randomly selected pairs of deformed spheres producing the result. However, in the reproduced figure, the relative performance of different methods is similar to that in Figure 8 in the paper. 

3) Due to the discrepancy in how area is computed, the x-axis of Figure 8 and 9 in the paper is not consistent with that of the other figures (Figure 1, 10~12), but the relative performance of different methods remains the same. In the code, we rescale the errors in Figure 8 and 9 by a factor of 1/sqrt(6), so that the reprodcued figures are the same as the ones presented in paper. 

Please let us know (rqhuang88@gmail.com) if you have any issue implementing this package. 