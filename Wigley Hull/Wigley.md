## Reconstruction of flow past a ship hull based on multi-fidelty data

Quick and accurate prediction of wave resistance is important during structural and geometrical design and loads evaluation. 
Methods such as potential-flow panel method and computational fluid dynamics method have been applied to predict ship waves and resistance. 
Potential-flow panel method solves the boundary integral equation based on the use of a Green function such as Havelock source which satisfies the boundary condition, providing efficient tool for routine applications. 
However, this kind of method leads to a loss of accuracy due to the assumption that the fluid is inviscid. 
CFD method solves the Navier-Stokes equations that include most of the relevant flow physics. 
But CFD method is computationally expensive, which do not meet the requirement of efficiency. 

In this work, multi-fidelity Gappy proper orthogonal decomposition (POD) based method was applied to predict pressure distribution of a Wigley ship hull. 
Results were post-proceeded to obtain wave resistance and validated by the experimental values. 
Gappy POD method gave a more accurate prediction than panel method and reduced computational effort compared to CFD simulation.
The optimum number of modes for reconstruction and the applicability for both interpolation and extrapolation was also evaluated in this work. 


### 1. Data collection:
   Implement panel method (based on invisid equaiton) to generate low-fidelity data and acquire samples with higher fidelity using STAR-CCM+ (based on N-S equation). 
  <div align="center">
	<img src="https://github.com/h5ch/GPOD_Burgers/blob/main/Wigley%20Hull/imgs/wave_pressure.png" alt="Editor" width="800">
  </div>

    fig.1 Comparison of low-fidelity and high-fidelity data at Fr=0.267: (a) wave profile; (b) pressure coefficient distribution 
    obtained from panel method; (c) pressure coefficient distribution obtained from STAR-CCM+.



### 2. Data reconstruction:
   Apply Gappy POD to reconstruct pressure coefficient distribution. 
  <div align="center">
	<img src="https://github.com/h5ch/GPOD_Burgers/blob/main/Wigley%20Hull/imgs/reconstruction.png" alt="Editor" width="800">
  </div>
  
    fig. 2 Reconstruction of Gappy POD: (a) reconstructed pressure coefficient distribution; (b) reconstruction error. 


### 3. Number of snapshot:
  Use three different sets of training data which includes different numbers of snapshots. 
  <div align="center">
    
  | Set  | Training |
  | ------------- | ------------- |
  | 9 snapshots   | 0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4  |
  | 5 snapshots   | 0.2,0.25,0.3,0.35,0.4  |
  | 3 snapshots   | 0.2,0.3,0.4  |
  
   </div>

  <div align="center">
	<img src="https://github.com/h5ch/GPOD_Burgers/blob/main/Wigley%20Hull/imgs/pressure_coefficient267.png" alt="Editor" width="800">
  </div>
  
    fig. 3 Pressure profile predicted by Gappy POD at $Fr=0.267$ with (a-c): 9 snapshots; (d-f): 5 snapshots; (g-i): 3 snapshots
    at the bow, middle, and stern of the Wigley hull respectively. 

  <div align="center">
	<img src="https://github.com/h5ch/GPOD_Burgers/blob/main/Wigley%20Hull/imgs/eigen.png" alt="Editor" width="800">
  </div>
  
    fig. 4 Reconstruction error as a function of the number of POD modes and corresponding minimum eigenvalues. (a)(d): 3 snapshots; 
    (b)(e): 5 snapshots; (c)(f): 9 snapshots. 

### 4. Comparison of interpolation test case and extrapolation test case: 
  Give prediction of interpolation test case ($Fr \in [0.2,0.4]$) and extrapolation test case ($Fr \notin [0.2,0.4]$).

  <div align="center">
	<img src="https://github.com/h5ch/GPOD_Burgers/blob/main/Wigley%20Hull/imgs/pdfe.png" alt="Editor" width="800">
  </div>
  
    fig.5 The p.d.f. of normalized error obtained from Gappy POD with different numbers of modes for (a): interpolation; 
    (b): extrapolation.
