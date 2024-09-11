# GPOD_Burgers
A one-dimensional burgers equation case for Gappy POD. 

In this case, Gappy POD was applied for multi-fidelity flow reconstruction. The process of reconstruction is summarized as follows: 

1. Data generation:

   The viscous Burgers equation is written as:

  $$ 
  \frac {\partial u} {\partial t} + \frac {\partial } {\partial t} (\frac {1}{2} u^2) = \frac {1}{Re} \frac{\partial^2 u}{\partial x^2}
  $$

   An analytical solution [4] to standard Burgers equation was used to generated high-fidelity samples. 

  $$
  u(x,t) = \frac {x/(t+1)} {1 + \sqrt{(t+1)/A_0} exp{Re[x^2/(4t+4)]}}
  $$

  Then the viscous term in the Burgers equation was omitted to generate low-fidelity data, and the invisid Burgers euation was solved using fifth-order TENO scheme [2]. 

2. Data reconstruction:

   Use SVD to find the POD basis of snapshot which was used to reconstruct the solution of invisic Burgers equation with improved accuracy. 

## References
[1] Pawar, S., & San, O. (2019). CFD Julia: A learning module structuring an introductory course on computational fluid dynamics. Fluids, 4(3), 159.

[2] Fu, L., Hu, X. Y., & Adams, N. A. (2016). A family of high-order targeted ENO schemes for compressible-fluid simulations. Journal of Computational Physics, 305, 333-359.

[3] Maleewong, M., & Sirisup, S. (2011). On-line and Off-line POD Assisted Projective Integral for Non-linear Problems: A Case Study with Burgers-Equation. International Journal of Mathematical and Computational Sciences, 5(7), 984-992.

[4] Li, T., Buzzicotti, M., Biferale, L., Bonaccorso, F., Chen, S., & Wan, M. (2023). Multi-scale reconstruction of turbulent rotating flows with proper orthogonal decomposition and generative adversarial networks. Journal of Fluid Mechanics, 971, A3.
