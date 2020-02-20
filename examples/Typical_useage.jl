################## Typical basic usage ######################
#
# use these two funtions to call simple simulation of the experiments
#   ysz_fitting.EIS_simple_run()
#   ysz_fitting.CV_simple_run()
#
# parameters:
#   pyplot        >>> 0 = plot nothing, 1 = finall plot of the experiment, 2 = plot details
#   prms_names    >>> in array or tuple of strings specify, which parametr of the model should be changed
#   prms_values   >>> spedify, to which value it will be changed
#       note: You can change any of the parameters in YSZParameters class by typing its name, e.g. x_frac = 0.50
#       note: You can do small parametrical study by typing an array of values instead of one value
#             For example:  prms_names = ["DGA", "DGR", "DGO"]
#                           prms_values = [0.4, 0.5, 0.3]                       # 1 set of parameters
#                           prms_values = [0.4, [0.4, 0.5, 0.6], 0.3]           # 3 sets of parameters
#                           prms_values = [0.4, collect(0.4 : 0.1 : 0.6), 0.3]  # 3 sets of parameters
#   TC, pO2, (EIS_bias=0.0) >>> the experimental setting
#       note: You can also write an array instead of 1 value for this parameters
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

includet("ysz_fitting.jl")


ysz_fitting.EIS_simple_run(pyplot=1, prms_names=["betaA", "betaR", "betaO", "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], prms_values=[0.5, 0.5, 0.5,       0, 0, 0,       13.8, 16.9, 17.95,        0.65, -0.9, collect(-0.4 : -0.05 : -0.6),     9.05e-13], TC=850, pO2=20)

ysz_fitting.EIS_simple_run(pyplot=1, prms_names=["betaA", "betaR", "betaO", "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], prms_values=[0.5, 0.5, 0.5,       0, 0, 0,       13.8, 16.9, 17.95,        0.65, -0.9, -0.54,     9.05e-13], TC=850, pO2=[20, 40])


ysz_fitting.CV_simple_run(pyplot=1, prms_names=["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"], prms_values=[13.8, 16.9, 17.95,        0.65, -0.9, -0.54,       9.05e-13], TC=850, pO2=20)

ysz_fitting.CV_simple_run(pyplot=1, prms_names=["betaA", "betaR", "betaO", "SA", "SR", "SO", "expA", "expR", "expO", "A0", "R0", "K0",  "DGA", "DGR", "DGO", "DD"], prms_values=[0.5, 0.5, 0.5,       0, 0, 0,      0, 0, 0,       15.0, 18., 18.0,              0.4, -0.4, [-0.3, -0],          9.05e-13], TC=850, pO2=40, show_experiment=true)
