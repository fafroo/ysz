includet("../examples/ysz_fitting.jl")
name=ARGS[1]
EIS_hypercube, prms_lists, scripted_tuple, prms_names = ysz_fitting.import_par_study_from_metafile(name=name, EIS_bool=true);
EIS_err = ysz_fitting.EIS_get_error_projection_to_prms(EIS_hypercube, prms_lists, );
ysz_fitting.display_the_best(EIS_err, EIS_hypercube, prms_lists, 10);
ysz_fitting.display_err_projection(EIS_err, prms_lists, prms_names, 200);


# ysz_fitting.EIS_simple_run(pyplot=true, prms_names=["A0", "R0", "K0","SA","SR","SO", "DGA", "DGR", "DGO", "betaA", "betaR", "betaO", ], prms_values=[19.8, 19, 16.828,    1, 1, 1,    -0.5, 0.2, 0.2,      0.5, 0.5, 0.5,  ], dx_exp=-8)
