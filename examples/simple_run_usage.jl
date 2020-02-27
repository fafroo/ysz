################## Typical basic usage ######################
#
# use this funtion to call simple simulation of the experiments
#   ysz_fitting.simple_run()
#
# parameters:
#   TC, pO2, (bias=0.0) >>> the experimental setting
#         note: experimental data are for
#                 TC \in    {700, 750, 800, 850}
#                 pO2 \in   {0, 20, 40, 60, 80, 100}
#                 bias \in  {-1.0, -0.5, 0.0, 0.5, 1.0} but for 700 is missing {-1.0, 1.0}
#         note: You can write an array instead of 1 value for this parameters (as for prms_values described below)
#   simulations     >>> =["EIS"] or =["CV"] or =["EIS", "CV"] for perform both methods for all experimental settings
#   pyplot          >>> 0 = plot nothing, 1 = finall plot of the experiment, 2 = plot details
#   prms_names      >>> in array or tuple of strings specify, which parametr of the model should be changed
#   prms_values     >>> spedify, to which value it will be changed
#         note: You can change any of the parameters in YSZParameters class by typing its name, e.g. x_frac = 0.50
#         note: You can do small parametrical study by typing an array of values instead of one value
#              For example:  prms_names = ["DGA", "DGR", "DGO"]
#                            prms_values = [0.4, 0.5, 0.3]                       # 1 set of parameters
#                            prms_values = [0.4, [0.4, 0.5, 0.6], 0.3]           # 3 sets of parameters
#                            prms_values = [0.4, collect(0.4 : 0.1 : 0.6), 0.3]  # 3 sets of parameters
#   show_experiment >>> if the experimental setting is in the database, it shows fitting error and (if pyplot > 1) plot experimental data
#   fig_size        >>> standard PyPlot tuple (x,y)
#
# side-notes:
#   - this is using LoMA model
#   - expA \in {0, 1} determines, if the kinetics of the reaction will be switched to pure exponential form (from LoMA)
#   - "L" inductance
#   - TC = (700, 750, 800, 850)°C  => use => DD = ( 1.277, 2.92, 5.35, 9.05)e-13
#       - I know it should influence also the other parameters like \nu, but this is the first approximation to match overall resistance of the system
##############################################################

includet("../src/ysz_fitting.jl")

# example 1
ysz_fitting.simple_run(TC=850, pO2=[20, 40], simulations=["EIS", "CV"], pyplot=1, 
prms_names=["betaA", "betaR", "betaO", "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[0.5, 0.5, 0.5,       0, 0, 0,       21.0, 20.0, 20.5,        0.7, 0.0, 0.0,     9.05e-13])

# example 2 (this is equivalent to example 1)
ysz_fitting.simple_run(TC=850, pO2=[20, 40], simulations=["EIS", "CV"], pyplot=1, 
prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[      21.0, 20.0, 20.5,        0.7, 0.0, 0.0,     9.05e-13])

# example 3
ysz_fitting.simple_run(TC=850, pO2=20, simulations=["EIS"], pyplot=1, 
prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], 
prms_values=[      collect(20.5 : 0.5 : 21.5), collect(19.5 : 0.5 : 20.5), 20.5,        0.7, 0.0, 0.0,     9.05e-13])
