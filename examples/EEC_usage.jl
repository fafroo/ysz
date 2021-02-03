includet("../src/ysz_fitting.jl")

# view_EEC function for trying several EECs ...
ysz_fitting.view_EEC(
           EEC_structure="RL-R-L-RCPE-RCPE",
           prms_values=[0.7, 0.0011,      1.7, 0,        0. , 0.001, 1.0,    1.0, 0.01, [0.8, 1.]], 
           f_range=(0.01, 10000, 1.2), # (lowest frequency, highest frequency, growing_factor)
           print_bool=true,
           plot_legend=true);
           
               
# fixed_prms_values example              
computed_data = ysz_fitting.run_EEC_fitting(TC=[700,], pO2=[100], bias=collect(-0.0 : -0.2 : -0.4), data_set=["MONO_110"], 
           #   extra_tokens=Dict("rep"=>[1,2], "shaky_table" => ["on", "off"]),
              f_interval="auto", succes_fit_threshold=0.007, trim_inductance=false,
              #init_values = [0.69402504, 1.6663523, 0.033978099, 0.05, 0.8, 0.012615565, 0.05, 0.9],
              fixed_prms_names=["L2"], fixed_prms_values=[1.75],
              alpha_low=0.2, alpha_upp=1,            
              which_initial_guess="both",
              save_file_bool=false, file_name="test_alg.txt", save_R1_file=false,
              plot_bool=true, plot_legend=true, plot_best_initial_guess=false, plot_fit=true,
              show_all_initial_guesses=true,
              use_DRT=false);

### Notes:
### - extra_tokens -> you can specify additional parameters to import specific file. For example. 
###         (i) repetitions of experiment ... "rep" => [1,2]
###         (ii) if there was some shaky_table ... "shaky_table" => ["on", "off"]  
###               (**) you can specify string of characters, not only numbers
###   * syntax for this parameter in terminal is then
###       -> extra_tokens=Dict("rep"=>[1,2],     "shaky_table" => ["on", "off"])
###   * and you can use this additional parameters in defining the path in "import_experimental_data.jl" file. For example
###       -> fNAME = string("(........)/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(extra_tokens["rep"]).z")
###     or with shorter form
###       -> fNAME = string("(........)/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(ET["rep"]).z")
###     In this case, all predefined possibilities for "rep" will be investigated, i.e. [1,2]
###   * if you do not use some extra token (i.e. "shaky_table") in the importing path, both variants (i.e. ["on", "off"]) will be   
###     anyway investigated (fitted) and written in the output file. But the investigated data will be exactly the same for both cases.
###     Which means... meaningless waste of computation resoures :)
### 
### - run_EEC_fitting works at the time for structrure "R-L-{(RCPE-)_Ntimes}{-C}" EEC
### - parameters are ["R1", "L2", "R3", "C3", "alpha3", "R4", "C4", "alpha4", ..., C(2 + 3*N + 1)]
### - if initial values are not specified, algorithm chooses quite good initial values by itself for
###   - supported structure is "R-L-{(RCPE-)_Ntimes}{-C}" EEC
###     - "R-L" is obligatory
###     - {(RCPE-)_ntimes} means for example "RCPE-RCPE-RCPE-". ---> can be even 0-times
###     - {-C} at the end is optional parametr for capacitance
### - f_interval crop frequencies used from experimental measurement. 
###       e.g. f_exp = [1, 2, 3, 4, 5, 6, 7, 8] and f_interval = [3, 6], then only [3, 4, 5, 6] frequencies are used 
###       NOTE: f_interval = "auto" ------> (not bad) f_interval is found automatically
###                                 ------> this is the default option 
### - trim_inductance -> trims all "positive Im(Z) points" in high frequencies except of one to estabilish the intersection of Z with real axis.
### - succes_fit_threshold decides which fitting error are too big and a warning should be displayed
### - plot_bool: if Nyquist should be plotted
### - plot_legend: ... clear :)
### - which_initial_guess ... initial guesses are
###                                   - algorithmic initial guess
###                                   - previous fitted values (unless it is the first fitted Nyquist :) )
###                       = "alg"       ---> only use the algorithmic one
###                       = "previous"  ---> use the previous near bias if it is possible
###                       = "both"      ---> use both and take the better one
### - show_all_initial_guesses ... plot all initial guesses together with its fitted values 
###                            ... print to terminal values of these
### - plot_best_initial_guess ... clear as well [:
### - plot_best_fit ... plots the best fit :)
### - save_file_bool decides, if the output is printed into terminal (=false) or is saved to a file (=true)
### - default values are << save_to_folder="../data/EEC/" >> and  << file_name="default_output.txt" >>
### - save_R1_file -> if true, it saves file with only one column with R1 values. The name is "$(file_name)_R1" (before the extension)
### - fixed_prms_names: defines which parameters should NOT be fitted. e.g. ["R1", "L2"]
### - fixed_prms_values: defines the values of the fixed_prms_names. e.g. [0.4, 1.4e-6]    
### - alpha_low=0.2, alpha_upp=1.0 (default values) defines lower and upper bound for fitting alphas         

# saving and loading
ysz_fitting.save_EEC_data_holder(computed_data, folder="../data/EEC/", file_name="default.txt")
computed_data2 = ysz_fitting.load_EEC_data_holder(folder="../data/EEC/", file_name="default.txt")


                      
# view of the results
R1_data = ysz_fitting.plot_EEC_data_general(computed_data, x_name="bias", y_name="R1",
                                            TC_interval = [700, 850],
                                            pO2_interval = [0, 100],
                                            bias_interval = [-Inf, Inf],
                                            data_set = Nothing,                              
                                            #
                                            fig_num=102,
                                            plot_legend=true, plot_all_prms=true);
            
### Notes:
### - computed_data are precomputed results
### - x_name defines what should be on x-axis
### - y_name defines what should be on y-axis
### - *_interval serve as a crop tool for displayed data. 
###       NOTE: It is an ""inteval"" for cropping. It is NOT a list of values to be plotted.
### - if data_set = Nothing, than all data_sets in computed_data are ploted
### - plot_all_prms: write all plotted parameters in the title of the figure
###
### !!! the function returns the last ploted set of y_values !!!
###

ysz_fitting.display_fit_vs_exp(computed_data, 
                                   TC = [700],
                                   pO2 = [60],
                                   bias = [-0.05],
                                   data_set = ["MONO_110"],
                                   f_interval="auto", use_DRT=true)
                                   
### Notes:
### f_interval="auto" behaves the same as in ysz_fitting.run_EEC_fitting( ... )

############################################################
############## Temperature fitting stuff ###################
TC_data = [700, 750, 800, 850]
fit_info = findfit_from_for_R(TC_data, R1_data)
fit_minimizer = fit_info.fit_minimizer



# concatening data via bias!
ysz_fitting.get_joint_EEC_data_via_bias(EEC_data_1, EEC_data_2)

