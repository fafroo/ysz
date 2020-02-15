includet("../examples/ysz_fitting.jl")
EIS_hypercube, prms_lists, scripted_tuple, prms_names = ysz_fitting.import_par_study_from_metafile(name="prvni_vrh_exp_ads_2_level_ONLY_EIS", EIS_bool=true);
EIS_err = ysz_fitting.EIS_get_error_projection_to_prms(EIS_hypercube, prms_lists);
ysz_fitting.display_the_best(EIS_err, EIS_hypercube, prms_lists, 10);
ysz_fitting.display_err_projection(EIS_err, prms_lists, prms_names, 100);
