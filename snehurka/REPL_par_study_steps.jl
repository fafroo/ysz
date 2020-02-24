includet("../examples/ysz_fitting.jl")

par_study_name=par_study_name

EIS_hypercube, prms_lists, scripted_tuple, prms_names, TC, pO2, EIS_bias = ysz_fitting.import_par_study_from_metafile(name=par_study_name, EIS_bool=true);
EIS_err = ysz_fitting.EIS_get_error_projection_to_prms(EIS_hypercube, prms_lists, TC, pO2, EIS_bias);
ysz_fitting.display_the_best(EIS_err, EIS_hypercube, prms_lists, TC, pO2, EIS_bias, 10);
ysz_fitting.display_err_projection(EIS_err, prms_lists, prms_names, 200);


#ysz_fitting.EIS_simple_run(pyplot=true, prms_names=["A0", "R0", "K0", "SA", "SR", "SO", "DGA", "DGR", "DGO", "betaA", "betaR", "betaO", "DD"], prms_values=[19.7, 19.7, 18.6,    1, 1, 1,    0.7, -0.8, -0.3,      0.5, 0.5, 0.5,    5.35e-13])
