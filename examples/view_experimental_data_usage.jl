################## Typical basic usage ######################
#
# 
#
##############################################################

includet("../src/ysz_fitting.jl")

# CV ... parameters are (TC, pO2)
ysz_fitting.CV_view_experimental_data([800, 850], [20, 40])

# EIS ... 
ysz_fitting.EIS_view_experimental_data(TC=700, pO2=[0], bias=collect(-1 : 0.1 : -0.8), 
                   use_checknodes=false, data_set=["MONO_NEW"], plot_legend=true)
