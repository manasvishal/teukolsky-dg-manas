# Teukolsky_DG----> A discontinuous Galerkin solver written in MatLab for solving 1+1D reduced Teukolsky equation.


# How to run the code:
* First add the 00th folder and its subfolders to path in matlab. This contains all the scripts needed for the DG grid as well as some custom written scripts for different system. 
* Then change the directory to whatever system you want to solve. 
* Run the driver code associated with the system.


# Folders :
* [00 folder](/00_essential_codes/) : Contains all the essential scripts to run the solver. Must be added in the path!
* [01 folder](/01_advection_code_without_delta/) : First attemmpt to implement the 1D advection equation without a delta term.
* [02 folder](/02_wave_code_without_delta/) : A wave code solver with given initial conditions.
* [03 folder](/03_advection_code_with_delta/) : Introduce a delta source term in the 1D Advection equation.
* [04 folder](/04_wave_code_with_delta/) : Introduce a delta source term in the 1D wave equation.
* [05 folder](/05_scalar_wave_in_Schw/) : Solve perturbed scalar field equation in Schwarzschild.
* [06 folder](/06_scalar_wave_in_Schw_modified/) : Removed Global 1D dependence in the code -- much faster from now on.
* [07 folder](/07_scalar_wave_in_Kerr/) : First attempt at solving scalar wave in Kerr. Did not work out because of the coupling in RHS.
* [08 folder](/08_wave_code_with_rhs_utt/) : Perturbed the 1D wave equation by adding u_tt term in RHS and the solution was unstable.
* [09 folder](/09_advection_code_with_diff_wave_speed/) : Applied FD to the RHS source term in the previous case and the solution was unstable.
* [10 folder](/10_FD_test_for_coupled_equation/) : Not important. Can be deleted.
* [11 folder](/11_scalar_wave_in_kerr_lmax4/) : Matrix implementation of the coupled scalar wave equation with LF flux in Kerr. Modes include \ell=0,2 lmax=4.
* [12 folder](/12_scalar_wave_in_kerr_lmax6/) : Matrix implementation of the coupled scalar wave equation in Kerr. Modes include \ell=0,2,4 lmax=6.
* [13 folder](/13_scalar_wave_in_kerr_lmax8/) : Matrix implementation of the coupled scalar wave equation in Kerr. Modes include \ell=0,2,4,6 lmax=8.
* [14 folder](/14_m1_kerr_lmax6/) : Solves scalar wave eqn in Kerr with modes (l,m)=[(1,1),(3,1),(5,1)]
* [15 folder](/15_m2_kerr_lmax6/) : Solves scalar wave eqn in Kerr with modes (l,m)=[(2,2),(4,2),(6,2)]
* [16 folder](/16_m2_kerr_lmax8/) : Solves scalar wave eqn in Kerr with modes (l,m)=[(2,2),(4,2),(6,2),(8,2)]
* [17 folder](/17_hyp_layers_wave/) : First attempt at introducing hyperboloidal layers. Not important. Can be deleted.
* [18 folder](/18_hyp_layers_wave2/) : Working hyperboloidal layers in the ordinary 1D wave equation with a potential.
* [19 folder](/19_hyp_layers_sch/) : Introduced hyperboloidal layers in the scalar wave equation in the Schwarzschild geometry.
* [20 folder](/20_hyp_kerr_lmax4/) : Hyperboloidal layers in Kerr with modes l=0,2
* [21 folder](/21_hyp_kerr_lmax8/) : Hyperboloidal layers in Kerr with modes l=0,2,4,6
* [22 folder](/22_m2_hyp_kerr_lmax8/) : Hyperboloidal layers in Kerr with modes (l,m)=[(2,2),(4,2),(6,2),(8,2)]
* [23 folder](/23_m1_hyp_kerr_lmax9/) : Hyperboloidal layers in Kerr with modes (l,m)=[(2,1),(4,1),(6,1),(8,1)]
* [24 folder](/24_m2_hyp_symgrid_lmax8/) : Attempted symmetric hyperboloidal layers in the grid. Not working properly as of July 19th.
* [25 folder](/25_source_tern_ell4em2/) : Added the source term in Kerr with modes (l,m)=[(2,2),(4,2)] and generalized the matrix. Uses an adaptive grid allowing the user to place a different number of subdomains to the left/right of the particle. Note : the potentials and matrices are wrong!!
* [26 folder](/26_verify_hyp_layers/) : Several tests to study how the solution behaves depending on the hyperboloidal layers' parameters
* [27 folder](/27_ori_khanna_test/) : Script to study the experiment suggested by Prof. Ori and Prof. Khanna
* [28 folder](/28_scalar_fluxes_one_mode) : Takes the 18th folder as base where we use the symmetric hyperbolic approach to solve the wave equation. Here instead we use the Schwarzschild potential and a source term.
* [29 folder](/29_scalar_fluxes_sch_two_modes) : Takes the 25th folder as base and change the matrices for the new symmetric hyperbolic approach.


#


* Books and reviews:
  * [Nodal Discontinuous Galerkin Methods](https://www.springer.com/gp/book/9780387720654) by Jan Hesthaven, and Tim Warburton.

* Notes and derivations
  * [notes](https://drive.google.com/drive/folders/1BHSqf32-VC3zhgFmNFhq1b9BqqWLfw-t?usp=sharing)
  
* For further queries and/or suggestions, email vishalmanas28@gmail.com or contact on twitter [@manas_vishal](https://twitter.com/manas_vishal) 
