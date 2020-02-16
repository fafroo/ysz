includet("../examples/ysz_fitting.jl")
EIS_hypercube, prms_lists, scripted_tuple, prms_names = ysz_fitting.import_par_study_from_metafile(name="prvni_gas_ONLY_EIS", EIS_bool=true);
EIS_err = ysz_fitting.EIS_get_error_projection_to_prms(EIS_hypercube, prms_lists, );
ysz_fitting.display_the_best(EIS_err, EIS_hypercube, prms_lists, 10);
ysz_fitting.display_err_projection(EIS_err, prms_lists, prms_names, 100);

println("EIS_hypercube, prms_lists, scripted_tuple, prms_names = ysz_fitting.import_par_study_from_metafile(name=\"prvni_gas_ONLY_EIS\", EIS_bool=true);")
println("EIS_err = ysz_fitting.EIS_get_error_projection_to_prms(EIS_hypercube, prms_lists, );")
println("ysz_fitting.display_the_best(EIS_err, EIS_hypercube, prms_lists, 10);")
println("ysz_fitting.display_err_projection(EIS_err, prms_lists, prms_names, 100);")
println("optional >>> ysz_fitting.EIS_test_checknodes_range(omega_range=ysz_fitting.EIS_get_shared_omega_range())")
