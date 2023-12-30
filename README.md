# D-RBF-PU-method for solving an elliptic PDE

    By: Davoud Mirzaei, Uppsala University, Sweden 

------------

This file contains the Matlab code for the D-RBF-PU method of
 
"D. Mirzaei, The Direct Radial Basis Function Partition of Unity (D-RBF-PU) Method for Solving PDEs, SIAM J. Sci. Comput. (2021)"

------------

Execute 'StartRun.m' to see the results.

------------
Some points: 

- Initially, this code provides the results for a 2D example but all functions except PolyMat and ScatPoints2D work in all dimensions. 

- The similar functions ScatPoints1D and ScatPoints3D can be simply developed by the user.

- To generalize PolyMat to other dimensions you just need to switch over different cases for the MultiIndex vector. Or, you may use an integer partitioning algorithm to produce this vector in arbitrary dimensions for a given polynomial order. 

- The RBF used in the code is a polyharmonic spline kernel, leveraging the scalability property inherent to this kernel type. Nevertheless, users have the flexibility to incorporate alternative kernels into the Frbf.m function.

- The code addresses the elliptic PDE \( -\Delta u + u = f \) with both Dirichlet and Neumann boundary conditions. Nonetheless, users can adapt it for other types of PDEs with minimal modifications to the StartRun.m function.  

- Compared to the results of the paper, the present code assigns patch radii to maintain a constant number of points within each patch.

------------
