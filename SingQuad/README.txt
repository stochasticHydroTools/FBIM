main driver: gausselim_dr.f
main computational code: stokesGE.f

Gaussian elimination based
code for solving Stokes equation with velocity 
boundary conditions in simply/multiply, interior/exterior
problems.

The input file name is now geom_input.
The first line is still k0,k
k0 =0 is for interior problems and k0 = 1
k is the number of obstacles (excluding the outer boundary for the interior case)
k is the number of obstacle components which are assumed to be within
the outer boundary for the multiply connedcted case. For the exterior
case they are assumed to be everywhere.

After the first line, the code would expect k-k0+1 lines of input with 5 inputs each:

number of points to discretize the boundary, length of semi major axis, length of semi minor axis, x coordinate of center of ellipse, y coordinate of center of ellipse

If it is an interior problem, the outer boundary should be specified first.

The generates the velocity on the boundary due to a known
bounded solution in the interior/exterior and solves the velocity
boundary value problem.

To run the code, the make file is makestokesn

Details about stokesGE.f:

v here is represented as a vector in the following form:
-vy1, vx1, -vy2, vx2, ...
Same for single-layer, it is
-fy1, fx1, -fy2, fx2, ...

Delete and ignore the last row and column in these and the matrix, it corresponds to uniform motion.

In these codes ``single layer'' means
S_{\Gamma} \sigma - 1/(8 \pi) \int_{\Gamma} \sigma dS
where S_{\Gamma} \sigma is the true Stokes single layer potential.

Let swgt be the weight in assemblesl, then the single layer corresponds to

Single layer  = swgt * assemblesl + 0.25*swgt * assembleintw

All the parameters mentioned above and the Sherman Lauricella representation is for solving interior/exterior velocity boundary value problem for Stokes. The parameters above are the weights of different components in the representation. If you set dwgt=cwgt=pwgt = 0, the matrix formed will be just corresponding to single layer potential, which is then overwritten to invert it in the code. However you can just call the subroutine assemblesl separately to generate the matrix corresponding to the single layer.  

The subroutine assemblesl does the following. Given density \sigma = (\sigma_1 (x), \sigma_2(x)), the routine computes the matrix corresponding to the single layer for Stokes with density \sigma. The densities are ordered as (-\sigma_2(x_i), \sigma_1(x_i))  and the velocity on output is ordered as (-u_2(x_i), u_1(x_i)). 

I would recommend order = 16. The general entry of the matrix is K(i,j) w_j  where K(i,j) is the Green's function and w_j is the weight. If |i-j| >= stencil_size then w_j = trapezoidal weight and the stencil_size is directly proportional to order. The error in computing the single layer potential would be h^(order) where h is the trapezoidal spacing. 
