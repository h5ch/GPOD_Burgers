# GPOD_Burgers
A one-dimensional burgers equation case for Gappy POD. 

In this case, Gappy POD was applied for multi-fidelity flow reconstruction. The process of reconstruction is summarized as follows: 

### 1. Data generation:

The viscous Burgers equation is written as:

  $$ 
  \frac {\partial u} {\partial t} + \frac {\partial } {\partial t} (\frac {1}{2} u^2) = \frac {1}{Re} \frac{\partial^2 u}{\partial x^2}
  $$

An analytical solution [^1] to standard Burgers equation was used to generated high-fidelity samples. The solution represents a nonlinear wave propagating to the right with decreasing amplitude due to the viscosity effect.

  $$
  u(x,t) = \frac {x/(t+1)} {1 + \sqrt{(t+1)/A_0} exp{Re[x^2/(4t+4)]}}
  $$

Then the viscous term in the Burgers equation was omitted to generate low-fidelity data, and the invisid Burgers euation was solved using fifth-order TENO scheme [^2]. 
   
   ![image](https://github.com/h5ch/GPOD_Burgers/blob/main/Results/Burgers_solution.png)

### 2. Data reconstruction:

Use SVD to find the POD basis of snapshot which was used to reconstruct the solution of invisic Burgers equation with improved accuracy. 

   ![image](https://github.com/h5ch/GPOD_Burgers/blob/main/Results/Reconstruction.png)

*Note: This case also referenced [^3][^4]. 

## References
[^1]: Maleewong, M., & Sirisup, S. (2011). On-line and Off-line POD Assisted Projective Integral for Non-linear Problems: A Case Study with Burgers-Equation. International Journal of Mathematical and Computational Sciences, 5(7), 984-992.

[^2]: Fu, L., Hu, X. Y., & Adams, N. A. (2016). A family of high-order targeted ENO schemes for compressible-fluid simulations. Journal of Computational Physics, 305, 333-359.

[^3]: Pawar, S., & San, O. (2019). CFD Julia: A learning module structuring an introductory course on computational fluid dynamics. Fluids, 4(3), 159.

[^4]: Li, T., Buzzicotti, M., Biferale, L., Bonaccorso, F., Chen, S., & Wan, M. (2023). Multi-scale reconstruction of turbulent rotating flows with proper orthogonal decomposition and generative adversarial networks. Journal of Fluid Mechanics, 971, A3.
