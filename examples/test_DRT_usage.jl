################## Typical basic usage ######################
#
# three modes:  "EEC" - equivalent electric circuit of RC & RCPE in series
#                       R1, C1        = paralel RC circuit
#                       R2, C2, alpha = paralel R + CPE circuit
#               "exp" - experimetnal data
#               "sim" - simulation
# lambda - regularization parameter ... minimization of F = || A*x - b || + || lambda*I*x ||
#
##############################################################
includet("../src/ysz_fitting.jl")

#
ysz_fitting.test_DRT(mode="EEC", R1=0, C1=0.001, R2=1.0, C2=0.01, alpha=1.0, lambda=0.0)

#
ysz_fitting.test_DRT(mode="EEC", R1=0, C1=0.001, R2=1.0, C2=0.01, alpha=0.9, lambda=0.01)

#
ysz_fitting.test_DRT(mode="exp", TC=800, pO2=40, bias=0.0, lambda=0.001)

#
ysz_fitting.test_DRT(mode="sim", TC=800, pO2=100, bias=0.0,  
              prms_names=["expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu"     ], 
              prms_values=[1, 1, 1,       19.0, 19.8, 18.2,        0.0, -0.8, -0.0,     [90]*1.0e-14, [0.85]     ],
              lambda=0.001);

# new data_set called "POLY" for polycrystaline. The monocrystaline is "MONO"
ysz_fitting.test_DRT(mode="exp", TC=800, pO2=20, bias=collect(0.0 : 0.05 : 0.2), lambda=0.0, 
       backward_check=false, plot_option="DRT Nyq Bode Rtau", data_set="POLY", 
       tau_min_fac=1.1, tau_max_fac=1.1, tau_range_fac=20.0, f_range=(0.6, 9000, 1.2))
       
# richer options
ysz_fitting.test_DRT(mode="exp", TC=800, pO2=20, bias=0.0, lambda=0.0, 
       backward_check=false, plot_option="DRT Nyq Bode Rtau", data_set="POLY_OCV_test", 
       tau_min_fac=1.1, tau_max_fac=1.1, tau_range_fac=20.0, f_range=(0.6, 9000, 1.2))
       
