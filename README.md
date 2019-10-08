# Code supplement for "A P2P-dominant Distribution System Architecture"

This repository contains supplementary data and Julia code for "A P2P-dominant Distribution System Architecture" by Jip Kim and Yury Dvorkin.
The following packages must be installed:

  - Gurobi
  - Ipopt
  - JuMP
  - JLD    
  - MatpowerCases

This code includes two execution files:

SC.jl: Solving the proposed P2P model with the System-centric configuration
PC.jl: Solving the proposed P2P model with the Peer-centric configuration

To run the code, execute SC.jl/PC.jl, or include() it from a Julia prompt.

The input data is imported from MatpowerCases Package.

"AP15busDN" testsystem is given as follows in data folder:
  - Node.csv: Distribution network data
  - Line.csv: Distribution line data
  - Generator.csv: Distributed generation data
  - We also specify additional input data in the beginning of SC.jl/PC.jl

The results of the solve are saved as result.jld file using JLD(Julia Data format) package.
