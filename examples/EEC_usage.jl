includet("../src/ysz_fitting.jl")

# view_EEC function for trying several EECs ...
ysz_fitting.view_EEC(
           EEC_structure="RL-R-L-RCPE-RCPE",
           prms_values=[0.7, 0.0011,      1.7, 0,        0. , 0.001, 1.0,    1.0, 0.01, [0.8, 1.]], 
           f_range=(0.01, 10000, 1.2), # (lowest frequency, highest frequency, growing_factor)
           print_bool=true,
           plot_legend=true);
           

# run_EEC_fitting               
computed_data = ysz_fitting.run_EEC_fitting(TC=[700], pO2=[0], bias=collect(-1. : 0.1 : 1), data_set="MONO_110", 
               f_interval=[5, 20000], succes_fit_threshold=0.004,
               #init_values = [0.69402504, 1.6663523, 0.033978099, 0.05, 0.8, 0.012615565, 0.05, 0.9],
               save_file_bool=false, save_to_folder="../data/EEC/", file_name="monocrystaline.txt", with_errors=false,
               plot_bool=true, plot_legend=true, plot_initial_guess=false, plot_fit=true, use_DRT=false)        
               
               
# fixed_prms_values example              
computed_data = ysz_fitting.run_EEC_fitting(TC=[850], pO2=[100], bias=collect(-1.0 : 0.05 : -0.8), data_set="MONO_110", 
                      f_interval="auto", succes_fit_threshold=0.007,
                      #init_values = [0.69402504, 1.6663523, 0.033978099, 0.05, 0.8, 0.012615565, 0.05, 0.9],
                      fixed_prms_names=["L2"], fixed_prms_values=[1.75],
                      save_file_bool=false, save_to_folder="../data/EEC/", file_name="monocrystaline.txt",
                      plot_bool=true, plot_legend=false, plot_initial_guess=true, plot_fit=true, use_DRT=false);

### Notes:
### - run_EEC_fitting works at the time only for "R-L-RCPE-RCPE" EEC
### - parameters are ["R1", "L2", "R3", "C3", "alpha3", "R4", "C4", "alpha4"]
### - if initial values are not specified, algorithm chooses quite good initial values by itself
###
### - with_errors: this switch decides the way of fitting: 
###                     true  = LsqFit package for curve fitting with relative errors (standard deviations)
###                     false = Optim package -> general non-linear fit - no relative errors .... but works much better
### - f_interval crop frequencies used from experimental measurement. 
###       e.g. f_exp = [1, 2, 3, 4, 5, 6, 7, 8] and f_interval = [3, 6], then only [3, 4, 5, 6] frequencies are used 
###       NOTE: f_interval = "auto" ------> (not bad) f_interval is found automatically
###                                 ------> this is the default option 
### - succes_fit_threshold decides which fitting error are too big and a warning should be displayed
### - plot_bool: if Nyquist should be plotted
### - plot_legend: ... clear :)
### - plot_initial_guess ... clear as well [:
### - save_file_bool decides, if the output is printed into terminal (=false) or is saved to a file (=true)
### - default values are << save_to_folder="../data/EEC/" >> and  << file_name="default_output.txt" >>
### - fixed_prms_names: defines which parameters should NOT be fitted. e.g. ["R1", "L2"]
### - fixed_prms_values: defines the values of the fixed_prms_names. e.g. [0.4, 1.4e-6]                      
                      
# view of the results
ysz_fitting.plot_EEC_data_general(computed_data, x_name="bias", y_name="R1",
                                            TC_interval = [700, 850],
                                            pO2_interval = [0, 100],
                                            bias_interval = [-Inf, Inf],
                                            data_set = Nothing,                              
                                            #
                                            fig_num=102, 
                                            plot_legend=true);
            
### Notes:
### - computed_data are precomputed results
### - x_name defines what should be on x-axis
### - y_name defines what should be on y-axis
### - *_interval serve as a crop tool for displayed data. 
###       NOTE: It is an ""inteval"" for cropping. It is NOT a list of values to be plotted.
### - if data_set = Nothing, than all data_sets in computed_data are ploted
