## Reconstruction of flow past a ship hull based on multi-fidelty data

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
