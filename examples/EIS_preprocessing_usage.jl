################## Typical basic usage ######################
#
# EIS data preprocessing tool
#
##############################################################

include("../src/ysz_fitting.jl")




ysz_fitting.EIS_view_experimental_data(
                                        TC = collect(800 : 100 : 800),
                                        pO2 = [100],
                                        bias = collect(-0.2 : 0.05 : 0.2), 
                                        data_set = ["MONO_110"], 
                                        use_DRT = false
                                        ,fig_num=7
                                        ,EIS_preprocessing_control = ysz_fitting.EIS_preprocessing_control(
                                                f_interval="auto", 
                                                #add_inductance=0,
                                                #scale_factor=0.1,
                                                #trim_inductance=true, 
                                                leave_only_inductance_points=2,     # must be at least >= 2 
                                                outlayers_threshold=2.0,                                    
                                                use_DRT=false, DRT_control=ysz_fitting.DRT_control_struct()
                                         )
                                      );

# Notes to EIS_preprocessing_control
# 
# - f_interval 
#     -> "auto"
#     -> (10, 10000)  
# - add_inductance=1 will add L = 1 H to the Nyquist
# - scale_factor just scales the EIS data for both RE and IM part
# - trim_inductance leaves only one point under the real axis for high frequencies
# - "leave_only_inductance_points=x" leaves only (x) points above the real axis for high frequencies
#   -> x >= 2 
# - outlayers_threshold filters out the single "outlayers"
#     -> the average AVG difference between neighbouring EIS points is computed. 
#     -> then when a single point is more distant from its neigbours than "AVG * outlayers_threshold" number, it will be removed from EIS data
# - replace original data with DRT projection (to obtain more smooth data)
#     -> DRT_control can be specified
