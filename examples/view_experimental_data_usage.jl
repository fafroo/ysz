################## Typical basic usage ######################
#
# 
#
##############################################################

includet("../src/ysz_fitting.jl")

# CV ... parameters are (TC, pO2)
ysz_fitting.CV_view_experimental_data([800, 850], [20, 40])

# EIS ... parameters are (TC, pO2, bias)
ysz_fitting.EIS_view_experimental_data([800, 850], 0, 0.0)
