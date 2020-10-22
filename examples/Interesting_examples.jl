


# invariance test --- ysz_model_GAS_LoMA
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [0.001], 0.0, plot_option="RC"), pyplot=1,  
                     prms_names=[ "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu"     ], 
                     prms_values=[1,1,1,  collect(19 : 1.0 : 21), collect(19 : 1.0 : 21), collect(19 : 1.0 : 21),         0.0, -0.0, 0.0,     [90]*1.0e-14, [0.85]     ], 
                     use_experiment=false);

####################### recall, semicircles, Rf ########################                     
ysz_fitting.test_DRT(mode="EEC", 
           R_ohm=0.0, R1=1.0, C1=0.1, R2=1.0, C2=0.01, alpha=1.0, lambda=0.0, 
           backward_check=false, draw_semicircles=false, plot_option="Nyq Bode RC DRT", 
           tau_min_fac=50, tau_range_fac=5.0, f_range=(0.1, 90000, 1.10))                     
                     
                     
###################### CPE vs exp ##########################
# exp
ysz_fitting.test_DRT(mode="exp", 
           TC=850, pO2=[100], bias=[0.0], data_set="MONO_110", 
           lambda=0.0, 
           backward_check=true, draw_semicircles=true, plot_option="DRT RC Nyq Bode", 
           tau_min_fac=1, tau_range_fac=2.0, f_range=(0.6, 9000, 1.1));


# EEC .. zmenit tau_range_fac na 10
ysz_fitting.test_DRT(mode="EEC", 
           R_ohm=0.0, R1=0.0, C1=0.1, R2=1.0, C2=0.01, alpha=0.98, lambda=0.0, 
           backward_check=false, draw_semicircles=true, plot_option="Nyq Bode RC DRT", 
           tau_min_fac=50, tau_range_fac=2.0, f_range=(0.1, 90000, 1.10))
                     
###################### trendy v experimentalnich datech ################                     
ysz_fitting.test_DRT(mode="exp", 
                  TC=800, pO2=[20, 40, 60], bias=[0.0], data_set="MONO_110", 
                  lambda=0.0,
                  backward_check=false, draw_semicircles=false, plot_option="Rf RC Nyq Bode", 
                  tau_min_fac=1, tau_range_fac=2.0, f_range=(0.6, 9000, 1.1));                     
                     
                     
####################### peaky ########################                     
## 3 peaky pro ysz_shn.jl model  -> change "separate_vacancy": "true" on "false"
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, plot_option="DRT RC"), pyplot=1,  
                     prms_names=["kappaA", "kappaR", "kappaO", 
                                 "rA", "rR", "rO",     "rB", "rC",     
                                 "DGA", "DGR", "DGO",     
                                 "DD", "nu", "separate_vacancy"   ], 
                     prms_values=(1.0, 0.0, 0.0,
                                  collect(21.65 : 0.05 : 21.65), collect(20.0 : -0.1 : 20.0), 19.7,     10, 10,       
                                  collect(-0.5 : 0.05 : -0.5), 0.0, 0.05,    
                                  [90]*1.0e-14, [0.85], true    ), 
                     use_experiment=false);

## 4 peaky !!!!!!!  -> change "separate_vacancy": "true" on "false"
novy = 1; ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, plot_option="Nyq Rf RC"), pyplot=1,  
                     prms_names=["kappaA", "kappaR", "kappaO", 
                                 "rA", "rR", "rO",     "rB", "rC",     
                                 "DGA", "DGR", "DGO",     
                                 "DD", "nu", "sites_Om0", "sites_Om1", "separate_vacancy"    ], 
                     prms_values=(1.0, 0.0, 0.0,
                                  collect(21.65 : 0.05 : 21.65), collect(20.0 : -0.1 : 20.0), 19.7,     collect(18 : 0.05 : 18), 18,       
                                  collect(-0.5 : 0.05 : -0.5), 0.0, 0.05,    
                                  [90]*1.0e-14, [0.85], 2.0, collect(0.07 : 0.01 : 0.07), true    ), 
                     use_experiment=false);
                     
## 4 peaks -> clear process counterpart to one peak:
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, DRT_draw_semicircles=false, plot_option="Nyq Rf RC"), pyplot=1,  
                     prms_names=["kappaA", "kappaR", "kappaO", 
                                 "rA", "rR", "rO",     "rB", "rC",     
                                 "DGA", "DGR", "DGO",     
                                 "DD", "nu", "sites_Om0", "sites_Om1"     ], 
                     prms_values=(1.0, 0.0, 0.0,
                                  collect(21.85 : 0.3 : 21.85), 20.0, collect(19.0 : 0.1 : 20.),     collect(18 : 0.05 : 18), 18,       
                                  collect(-0.75 : -0.05 : -0.75), 0.0, 0.05,    
                                  [90]*1.0e-14, [0.85], 2.0, collect(0.07 : 0.01 : 0.07)     ), 
                     use_experiment=true);

## 4 peaks -> intersting exchange of roles with rB                      
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, plot_legend=false, DRT_draw_semicircles=false, plot_option="Nyq Rf RC"), pyplot=1,  
                            prms_names=["kappaA", "kappaR", "kappaO", 
                                        "rA", "rR", "rO",     "rB", "rC",     
                                        "DGA", "DGR", "DGO",     
                                        "DD", "nu", "sites_Om0", "sites_Om1"     ], 
                            prms_values=(1.0, 0.0, 0.0,
                                         collect(21.85 : 0.3 : 21.85), 20.0, collect(19.0 : 0.1 : 19.),  collect(17.0 : 0.1 : 20), 17,       
                                         collect(-0.75 : -0.05 : -0.75), 0.0, 0.05,    
                                         [12.60]*1.0e-11, [0.85], 1.0, collect(0.07 : 0.01 : 0.07)     ), 
                            use_experiment=true);
                     
################ ocasek #######################                  
## zatoceni low freq !!!!!!!!!!!!!!!!!!!
novy = 1; ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, 80, 0.0, DRT_draw_semicircles=false, plot_option="Nyq Bode Rf RC"), pyplot=1,  
                     prms_names=["kappaA", "kappaR", "kappaO", 
                                 "rA", "rR", "rO",     "rB", "rC",     
                                 "DGA", "DGR", "DGO",     
                                 "DD", "nu", "sites_Om0", "sites_Om1"     ], 
                     prms_values=(1.0, 0.0, 0.0,
                                  collect(22.25 : 0.1 : 22.25), collect(19.7 : 0.05 : 19.7), collect(20.38 : 0.01 : 20.38),    18,  collect(24.4 : 0.05 : 24.75),       
                                  collect(-1.05 : -0.05 : -1.05), 0.0, 0.05,    
                                  [90]*1.0e-14, 0.85, 2.0, collect(0.08 : 0.002 : 0.08)     ), 
                     use_experiment=false);
                     
## experiment -> ukazka trendu ocasku
ysz_fitting.test_DRT(mode="exp", TC=850, pO2=[0, 20, 40, 60, 80, 100], bias=[0.0], data_set="MONO_110", lambda=0.0, 
       backward_check=false, draw_semicircles=false, plot_option="Nyq Bode", tau_min_fac=1, tau_range_fac=20.0, f_range=(0.6, 9000, 1.1))
       
## experiment -> rozliseni skoro shodnych krivek
ysz_fitting.test_DRT(mode="exp", TC=800, pO2=[0, 20, 40, 60, 80, 100], bias=[0.0], data_set="MONO_110", lambda=0.0, 
       backward_check=false, draw_semicircles=false, plot_option="Nyq Bode", tau_min_fac=1, tau_range_fac=20.0, f_range=(0.6, 9000, 1.1))
                     
############### Solartron #####################
ysz_fitting.test_DRT(mode="exp", TC=750, pO2=[100], bias=[0.0], data_set="Solartron", lambda=0.0, 
       backward_check=true, draw_semicircles=true, plot_option="Rf RC Nyq Bode", tau_min_fac=1, tau_range_fac=10.0, f_range=(0.2, 64000, 1.1))

       
#####################   NEW NEW NEW  ###################
########################################################

# HebbWagner => nejaky trend tam je
ysz_fitting.test_DRT(mode="exp", TC=collect(600 : 20 : 720), pO2=[100], bias=[0.3], data_set="HebbWagner", lambda=0.0, 
              backward_check=false, draw_semicircles=false, plot_option="DRT RC Nyq Bode", tau_max_fac=2., tau_range_fac=10.0, f_range=(20., 65000, 1.1))
    
# Solartron ... peak_merge_tol
ysz_fitting.test_DRT(mode="exp", TC=collect(600), pO2=[100], bias=[0.3], data_set="TEST DRT - Zahner - dummy cell.z", lambda=0.0, 
              backward_check=false, draw_semicircles=true, plot_option="DRT RC Nyq Bode", 
              tau_max_fac=2., tau_range_fac=10.0, f_range=(1.2, 65000, 1.1), peak_merge_tol=0.5)
              
              
######## 2020-05-21 NEW  ###################
# zajimava anomalie pro jednu teplotu (ale je konzistentni pro okolni biasy)
ysz_fitting.run_EEC_fitting(TC=[700, 750, 800, 850], pO2=[100], bias=0.0, data_set="MONO_110", 
               f_interval=[40, 20000], succes_fit_threshold=0.004,
               #init_values = [0.69402504, 1.6663523, 0.033978099, 0.05, 0.8, 0.012615565, 0.05, 0.9],
               save_file_bool=false, file_name="monocrystaline.txt", with_errors=false,
               plot_bool=true, plot_legend=true, plot_initial_guess=false, plot_fit=true)

# dooost zvlastni Nyquisty !!!!!! :DDD
ysz_fitting.run_EEC_fitting(TC=[700, 750, 800, 850], pO2=[100], bias=1.0, data_set="MONO_110", 
               f_interval=[40, 20000], succes_fit_threshold=0.004,
               #init_values = [0.69402504, 1.6663523, 0.033978099, 0.05, 0.8, 0.012615565, 0.05, 0.9],
               save_file_bool=false, file_name="monocrystaline.txt", with_errors=false,
               plot_bool=true, plot_legend=true, plot_initial_guess=false, plot_fit=true)


#################### e_fac testovani ... spolu s B C legracemi ################
ysz_fitting.EIS_view_experimental_data( 
                                   TC = [750],
                                   pO2 = [0],
                                   bias = collect(-0.5 : 0.1 : 0.3),
                                   data_set = "MONO_110",
                                   )
                                   
# vliv e_fac
ysz_fitting.simple_run(ysz_fitting.EIS_simulation([850], 100,  collect(0.1 : 0.1 : 0.1), data_set="POLY"), pyplot=1,  
                            prms_names=["ARO_mode", 
                                        "kappaA", "kappaR", "kappaO", 
                                        "rA", "rR", "rO",         "rB", "rC",     
                                        "DGA", "DGR", "DGO",     
                                        "DD", "nu", 
                                        "separate_vacancy", "e_fac",       
                                        "sites_Om0", "sites_Om1"  
                                        ], 
                            prms_values=(false,
                                         0.0, 0.0, 0.0,
                                         21.5, 20.9, 20.7,        1, 20, #collect(21 : -0.5 : 21),       
                                         -0.0, -0.0, collect(-0.0 : 0.01 : 0.0),     
                                         [90].*1.0e-14, collect(0.85 : 0.05 : 0.85), 
                                         true, collect(0.0 : 0.01 : 0.05),       
                                         1/4, 1/2    
                                         ), 
                            use_experiment=false);
               
# podivne preklopeni pod realnou osu ... zmizeni standardniho RC oblouku
ysz_fitting.simple_run(ysz_fitting.EIS_simulation([850], 100,  collect(0.1 : 0.1 : 0.1), data_set="POLY"), pyplot=1,  
                            prms_names=["ARO_mode", 
                                        "kappaA", "kappaR", "kappaO", 
                                        "rA", "rR", "rO",         "rB", "rC",     
                                        "DGA", "DGR", "DGO",     
                                        "DD", "nu", 
                                        "separate_vacancy", "e_fac",       
                                        "sites_Om0", "sites_Om1"  
                                        ], 
                            prms_values=(false,
                                         0.0, 0.0, 0.0,
                                         21.5, 20.9, 20.7,        1, 1, #collect(21 : -0.5 : 21),       
                                         -0.0, -0.0, collect(-0.0 : 0.01 : 0.0),     
                                         [90].*1.0e-14, collect(0.85 : 0.05 : 0.85), 
                                         true, collect(0.02 : 0.0005 : 0.025),       
                                         1/4, 1/2    
                                         ), 
                            use_experiment=false);
                            
                            
                            
                            
########################################################################################################################
########################################################################################################################
########################################################################################################################
# Fitting ... 

# separate_vacancy ... best fit for {800 , 100, 0.0, OLD_MONO_100}
ysz_fitting.simple_run(
              ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false, DRT_draw_semicircles=false, plot_option="Nyq Bode Rf RC"),
              #ysz_fitting.CV_simulation(850, 80, 0.0),
              pyplot=1,  
                                   prms_names=["kappaA", "kappaR", "kappaO", 
                                               "rA", "rR", "rO",     "rB", "rC",     
                                               "DGA", "DGR", "DGO",     
                                               "DD", "nu", "separate_vacancy", "sites_Om0", "sites_Om1"     ], 
                                   prms_values=(0.0, 0.0, 0.0,
                                                21.79, 27.89, 19.74,   0, 0,
                                                0.097, 0.008, -0.218, #collect(-0.3 : 0.01 : 0.1),
                                                [9.3]*1.0e-11, [0.85], true, 1.2, 0.5  ), 
                                   use_experiment=true);

# separate_vacancy ... {800, [40, 60, 80, 100], 0.0, OLD_MONO_100}  ... BUT no CV fit
ysz_fitting.simple_run(
       ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false, DRT_draw_semicircles=false, plot_option="Nyq Bode Rf RC"),
       #ysz_fitting.CV_simulation(850, 80, 0.0),
       pyplot=1,  
                            prms_names=["kappaA", "kappaR", "kappaO", 
                                        "rA", "rR", "rO",     "rB", "rC",     
                                        "DGA", "DGR", "DGO",     
                                        "DD", "nu", "separate_vacancy", "sites_Om0", "sites_Om1"     ], 
                            prms_values=(0.0, 0.0, 0.0,
                                         21.29, 27.89, 18.54,   0, 0,
                                         0.097, 0.011, -0.218, #collect(-0.3 : 0.01 : 0.1),    
                                         [9.0]*1.0e-11, [0.85], true, 0.1, 0.5  ), 
                            use_experiment=true);

                                
# separate_vacancy ... {800, [40, 80], 0.0, OLD_MONO_100} .... with CVs !!!!!!!!!!!!!!!!!!!
ysz_fitting.simple_run(
       ysz_fitting.EIS_simulation(800, [40, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false, DRT_draw_semicircles=false, plot_option="Nyq Bode Rf RC"),
       #ysz_fitting.CV_simulation(800, [40, 80], 0.0),
       pyplot=1,  
                            prms_names=["kappaA", "kappaR", "kappaO", 
                                        "rA", "rR", "rO",     "rB", "rC",     
                                        "DGA", "DGR", "DGO",     
                                        "DD", "nu", "separate_vacancy", "sites_Om0", "sites_Om1"     ], 
                            prms_values= (0.0, 0.0, 0.0, 21.7474, 26.7855, 20.6466, 0.0, 0.0, 0.027508, 0.106477, -0.175718, 9.3e-11, 0.85, true, 2.60081, 0.5)
                            
                            ,use_experiment=true);
                            
# separate_vacancy ... {800, [40, 60, 80, 100], 0.0, OLD_MONO_100} ... with CVs !!!! 40 not quite good  ----> <<<"e_fac">>>
ysz_fitting.simple_run(
                     ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false, DRT_draw_semicircles=false, plot_option="Nyq Bode Rf RC"),
                     #ysz_fitting.CV_simulation(850, 80, 0.0),
                     pyplot=1,  
                                          prms_names=["kappaA", "kappaR", "kappaO", 
                                                      "rA", "rR", "rO",     "rB", "rC",     
                                                      "DGA", "DGR", "DGO",     
                                                      "DD", "nu", "separate_vacancy", "sites_Om0", "sites_Om1",  "e_fac"     ], 
                                          prms_values=(0.0, 0.0, 0.0, 21.411013641457583, 27.481801343839734, 19.774755928324968, 0.0, 0.0, -0.11464160143584941, -0.0684415679565699, 0.17278366094279843, 9.3e-11, 0.85, true, 1.014446180003061, 0.5, 0.16057108498298814)
                                          
                                          ,use_experiment=true);

                                          
                                          
                                          
                                          
# Hebb-Wagner 2 ... extended search... biases .... hm... 
Hebb_new = ysz_fitting.run_EEC_fitting(TC=[800], pO2=[2], bias=collect(0.3 : 0.05 : 0.6), data_set="HebbWagner_110", 
                      f_interval="auto", succes_fit_threshold=0.004,
                      #init_values = [6.689345, 0.0, 20.0, 0.009169754, 0.99551207, 76.7828, 0.00026448714, 0.71030164],
                      #fixed_prms_names=["L2"], fixed_prms_values=[0.0],
                      save_file_bool=false, save_to_folder="../data/EEC/", file_name="testov.txt",
                      plot_bool=true, plot_legend=false, plot_best_initial_guess=false, plot_fit=true, use_DRT=false);
                                          
                                          
                                          
                                          
                                          
                                          
######################################################################################################                                          
################ CORRECTED ######################## ysz_model_GAS_LoMA ###############################
#  ! ! ! ! ! ! use     ysz_model_GAS_LoMA !!!!   ->>> all good... 4 EIS_and CV  -->>> capturing the non-monotony !!!
##### EXP
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                            prms_names=["expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC"     ], 
                            prms_values=(1.0, 1.0, 1.0, 21.8106, 27.898, 20.4851, -0.697806, -0.434529, 0.229967, 9.3e-11, 0.85, 3.04504)
                            ,use_experiment=true);

# 4 CV ... convex fail ... but the best with EXP      
#  # # #  ##  # # #                                             ##  # # # ##  ##  ##  !!! !!! WEIRD discontinuous shift in EIS >>>
#
ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC"     ], 
                                   prms_values=(1.0, 1.0, 1.0, 27.345, 27.5299, 20.3249, -0.590819, 0.0340098, 0.758359, 9.3e-11, 0.85, 0.0374857)
                                   ,use_experiment=true);          
                                   
###### LoMA
########## # # # # # #
# 4 EIS ... not capturing the non-monotony ... even with extended search area ##### # # ## ## ## # # # 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC"     ], 
                                   prms_values=(0.0, 0.0, 0.0, 23.536793805295844, 26.44428245339161, 21.308492340939917, 0.5500221642340698, 0.029112401256050482, -0.020889037864824146, 9.3e-11, 0.85, 2.0593360762011694)
                                   ,use_experiment=true);
                                    
# 4 CV ... not capturing the shape ... convex fail
ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC"     ], 
                                   prms_values=(0.0, 0.0, 0.0, 24.115354988076593, 25.42144932757581, 22.139558838796805, -0.7999773681257123, -0.5118625772334413, -0.020915586937639292, 9.3e-11, 0.85, 0.5824816856107942)
                                   ,use_experiment=true);

################# e_fac --------- PAR STUDY -------------- ###################                                   
# e_fac ... R_ohm static ... example
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["e_fac", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC"     ], 
                                   prms_values=([0.0, 0.1, 0.2, 0.3, 0.4], 1.0, 1.0, 1.0, 20.8106, 27.898, 20.4851, -0.097806, -0.434529, 0.229967, 9.3e-11, 0.85, 3.04504)
                                   ,use_experiment=false);                                   
                                   
# e_fac ... R_ohm varying !!!! 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["e_fac", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC"     ], 
                                   prms_values=([0.0, 0.1, 0.2, 0.3, 0.4], 1.0, 1.0, 1.0, 23.8106, 27.898, 20.4851, -0.097806, -0.434529, 0.229967, 9.3e-11, 0.85, 3.04504)
                                   ,use_experiment=false);
# e_fac ... HF right turn
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["e_fac", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC"     ], 
                                   prms_values=([0.0, 0.1, 0.2, 0.3, 0.4], 1.0, 1.0, 1.0, 21.6106, 27.898, 20.4851, -0.497806, -0.234529, 0.229967, 9.3e-11, 0.85, 3.04504)
                                   ,use_experiment=false);
                                   

##################### FINALL CORRECT MODEL #######################
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                          prms_names=["expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC", "e_fac"     ], 
                                          prms_values=(0.0, 1.0, 0.0, 24.4631, 21.7136, 21.5347, -0.0113085, -0.0926479, 0.179689, 9.3e-11, 0.85, 3.70135, 0.0)
                                          ,use_experiment=true);
                                          
# EEL ... but pretty much the same as the above LEL                                          
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                          prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "e_fac"], 
                                          prms_values=(1.0, 1.0, 0.0, 24.1783, 21.7137, 21.5419, -0.0107187, -0.0862911, 0.185543, 9.3e-11, 0.85, 3.74939, 0.0)
                                          ,use_experiment=true);                                      
                    
                                   
                                   
# # # # # # 3 peaks with ARO separate_vacancy
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false, plot_option="Nyq Bode Rtau RC"), pyplot=1,  
                  prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(true,   
                                1.0, 1.0, 1.0,   
                                21.9783, 22.0137, 21.5419,    
                                -0.107187, -0.02911, -0.085543,  
                                9.3e-11, 0.85,      10., 15.,  0.0)
                  ,use_experiment=false);
           
# # # # EXP -- A.DG variations makes bigger semiarc
new = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
           prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "nus", "e_fac"],
            prms_values=(true,   
                         1.0, 1.0, 1.0,   
                         21.9783, 21.8137, 20.9419,    
                         [0.0, 0.2, 0.4, 0.6, 0.7], 0.72911, 0.285543,  
                         9.5e-11, 0.85, [8.0], [4.], 0.85, 0.0)
           ,use_experiment=false);
           
###### overvoltage doing looping :D
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                  prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(1.0, 1.0, 1.0, 21.9679, 21.5724, 20.9371, 0.293345, -0.241793, -0.05, 9.3e-11, 0.85, 6.1177, 27.7853, 0.3)
                  ,use_experiment=true);
                  
                  
                  
###############  R.r test ... HF ocasek
new_GAS_LoMA_shared = ysz_fitting.simple_run(
                                                ysz_fitting.EIS_simulation(800, [80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                #ysz_fitting.CV_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                
                                                    pyplot=1,  
                                                     prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(0.0, 0.0, 0.0, 22.5, collect(22.076839213753235 : 0.1 : 23.8), 27.69994335225577, 0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, 0.3033398196955781)
                                                     ,use_experiment=true);
                  
# EEE EIS CV ... almost fit (manually translated from EEL fit)
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                         prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(1.0, 1.0, 1.0, 21.8528, 21.6247, 20.8442, 0.322754, -0.120249, -0.0687233, 9.3e-11, 0.85, 1.75728, 10.9117, 0.0)
                         ,use_experiment=true);
                         
# EEE EIS CV ... relatively good fit for both CV and EIS 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                  prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(1.0, 1.0, 1.0, 21.94051183287039, 21.553329968593776, 20.965273571151613, 0.09841737344413018, -0.091875316601428, 0.04652480321433385, 9.3e-11, 0.85, 6.478278331551995, 7.546431173856936, 0.0)
                  ,use_experiment=true);
                  
# EEE EIS CV .. fitnes factos {1, 10} ... good fit for both ... only 2 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                  prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(1.0, 1.0, 1.0, 24.383816208538917, 21.94866414517787, 20.68161992890677, 0.1380699169577035, -0.22933373067502555, -0.23743440380456354, 9.3e-11, 0.85, 2.7481343742171562, 2.930552174417333, 0.0)
                  ,use_experiment=true);

# EEE EIS CV ... fitness factors {1, 1} ... very good EIS fit ... 3 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                  prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(1.0, 1.0, 1.0, 21.96912506564684, 21.569512049276458, 20.93977987957153, 0.1241112223034025, -0.10660629726868907, -0.006580383278506099, 9.3e-11, 0.85, 5.848806563958082, 8.27404376722011, 0.0)
                  ,use_experiment=true);

# EEE EIS .... CV fitted automatically :) ... very good fit for .. with e_fac !!                  
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                         prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(1.0, 1.0, 1.0, 21.9445, 21.3697, 20.8594, 0.101682, -0.178855, 0.115187, 9.3e-11, 0.85, 6.22922, 5.75616, 0.138391)
                         ,use_experiment=true);

# EEL EIS .... CV hooked ... sensitive on pO2 ... 3 peaks... with 3_fac
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                         prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(1.0, 1.0, 0.0, 21.9258, 21.2852, 21.7927, 0.0277078, -0.280821, 0.362658, 9.3e-11, 0.85, 7.54141, 1.65368, 0.296044)
                         ,use_experiment=true);
                         
# LEE EIS .... CV hooked
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(0.0, 1.0, 1.0, 24.253505704719995, 27.895308798638613, 20.689058014097917, 0.5435406241193925, 0.2124972805439068, -0.025572361570387683, 9.3e-11, 0.85, 2.1001141086263937, 0.8016033006905938, 0.004387336514558263)
                                ,use_experiment=true);

# ELE EIS .... CV hooked ... the worse sensitivity to pO2 ... and the worse shape of Nyquisty
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1.0, 0.0, 1.0, 22.004601154444064, 21.969771021412477, 20.64499552886603, -0.5553757496916619, -0.04999999999933978, 0.23448570623706932, 9.3e-11, 0.85, 5.925253670629901, 9.734714051365351, 0.7504653628307916)
                                ,use_experiment=true);

#### LLL CV .... husty!!! ... trefa
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(0.0, 0.0, 0.0, 24.2177, 21.6481, 22.6917, 0.698797, 0.252055, -0.122604, 9.3e-11, 0.85, 1.46538, 4.51906, 0.285749)
                                ,use_experiment=true);
                                
#### LLL EIS CV .... husty!!! ... obe nejsou spatni
# # new_GAS_LoMA_shared = ysz_fitting.simple_run([
# #                                   ysz_fitting.EIS_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
# #                                   ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
# #                                   ],
# #                                       pyplot=1,  
# #                                        prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
# #                                         prms_values=(0.0, 0.0, 0.0, 22.5013, 21.6143, 26.7452, 0.0289421, -0.05, 0.121954, 9.3e-11, 0.85, 1.17611, 8.74441, 0.324685)
# #                                        ,use_experiment=true);
#                                        
new_GAS_LoMA_shared = ysz_fitting.simple_run(
                                         #ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                         ysz_fitting.CV_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                         
                                             pyplot=1,  
                                              prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(0.0, 0.0, 0.0, 22.576839213753235, 22.011902293420093, 27.69994335225577, 0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, 0.3033398196955781)
                                              ,use_experiment=true);

#### LLL EIS CV .... zkus BEZ e_fac                                       
 # --->> nic moc :)
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(true, 0.0, 0.0, 0.0, 24.069480130382026, 27.84830853588366, 22.089454852605634, 0.7999985765425789, -0.5103390865926676, -0.022909285494503932, 9.3e-11, 0.85, 0.9567488020824775, 2.851374822442069, 0.0)
                           ,use_experiment=true);

#### EEE EIS .... shared vac 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                          prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC",  "ms_par","e_fac"     ], 
                                          prms_values=(false, 1.0, 1.0, 1.0, 21.8485, 20.5134, 20.9458, -0.05, -0.147443, -0.21487, 9.3e-11, 0.85, 19.4005, 9.96125, 1.31778)
                                          ,use_experiment=true);
                                          
#### LLL EIS .... shared vac
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                          prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC",  "ms_par","e_fac"     ], 
                                          prms_values=(false, 0.0, 0.0, 0.0, 23.6526, 20.9531, 22.7513, -0.163142, -0.300821, -0.40848, 9.3e-11, 0.85, 26.642, 9.97322, 0.919853)
                                          ,use_experiment=true);
                                          

### EEE EIS CV.... shared vac ... not bad, not the best ... 2 peaks
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                          prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC",  "ms_par","e_fac"     ], 
                                          prms_values=(false, 1.0, 1.0, 1.0, 21.6797, 21.0495, 20.7467, 0.0212226, -0.423638, 0.798779, 9.3e-11, 0.85, 6.75588, 9.95628, 0.305956)
                                          ,use_experiment=true);
                                          
### LLL EIS CV ... shared vac ...                                           
#garbage . . .
                           
### HAND FIT ### :) LLL EIS ... quite good EIS... bad CV
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 
                                                     22.277, 22.312, 22.3, 
                                                     0.048, -0.271, 0.566, 
                                                     1.05e-11, 0.35, 24., 2.544, 0.303)
                                       ,use_experiment=true);
                                       
#                                                    #################################################
### corrected HAND FIT ... LLL EIS ... WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOW WOWOWOWWOOW WO !!!!!! BRUTALISK !!!!! EIS cooool !!!!
#                                                    #################################################
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.325513444696398, 22.100276559559344, 22.382691084143286, 0.04243114428504254, -0.22279524777844634, 0.5430896356175956, 1.05e-11, 0.35, 34.397339021979015, 4.503994521917953, 0.338333333916094)
                                       ,use_experiment=true);
                                       
######### corrected HAND FIT for more pO2 ... LLL EIS 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true);
#### .... continuation from previos... to pO2 60
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.519825220050315, 22.101114333751948, 22.537917125481112, 0.046779285103424426, -0.036495661417036834, 0.6645252416469626, 1.06e-11, 0.35, 66.81244689470016, 2.899812097038787, 0.7357178655306754)
                                       ,use_experiment=true);
                           
##### ....... and ... then pO2 40, 60, 80,  again                           
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.49042860963425, 22.18758070122782, 22.535374543362902, 0.04626751671867485, -0.02797383209333109, 0.6910672027823139, 1.06e-11, 0.35, 51.10551528037439, 3.2209400174059213, 0.6775179360491119)
                                       ,use_experiment=true);                           
                           
#### ...... shared !!! 2 EIS CV ...not bad :)   >>>>>>> almost THE SAME result as for the separate vac
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(false, 0.0, 0.0, 0.0, 22.5615, 22.113, 27.0515, 0.0456777, -0.33161, 0.684794, 9.3e-11, 0.85, 75.0419, 9.5712, 0.311755)
                                              ,use_experiment=true);
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           

                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
###########################################################                           
#####################  PAR_STUDY  #########################                           
# separate VS shared vacancy ... 3 vs 2 peaks
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(true, 1.0, 1.0, 1.0, 
                                       21.6106, 21.898, 21.4851, 
                                       0.497806, collect(-0.3 : -0.02 : -0.3), 0.229967, 
                                       9.3e-11, 0.85, collect(1.0 : 1 : 14), 3.5, 0.0)
                           ,use_experiment=false);

                           
# e_fac descussion ... without inductance
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(0.0, true, 0.0, 0.0, 0.0, 
                                       23.6106, 22.598, 23.4851, 
                                       -0.497806, collect(-0.3 : -0.02 : -0.3), 0.229967, 
                                       9.3e-11, 0.85, 1.0, 1.0, collect(0.0 : 0.1 : 0.5))
                           ,use_experiment=false);
                           
                           
                           
                           
                           
                           
# OC study ... 3 clear peaks  LLL   !!!!                       
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(0.0, [true, false], 0.0, 0.0, 0.0, 
                                       22.2106, 22.098, 21.5851, 
                                       -0.1497806, 0.0, 0.029967, 
                                       9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.1 : 0.0))
                           ,use_experiment=false); 
# separate_vacancy study ... 3 clear peaks ... LLL !!!
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, [true, false], 0.0, 0.0, 0.0, 
                                              22.2106, 21.798, 21.1851, 
                                              -0.1497806, -0.0, -0.19967, 
                                              9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.1 : 0.0))
                                  ,use_experiment=false);
                           
# OC study ... 3 clear peaks ... EEE !!!!                          
                           ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                              22.1106, 22.098, 21.5851, 
                                              0.1497806, -0.0, 0.029967, 
                                              9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                  ,use_experiment=false);
                           
                           
########################################################################################################                           
#  .... .. .. 4 clear peaks ... EEE !!! bulk peak!
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                       prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                       prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                                   21.2, 21., 21.0, 
                                                                   0., 0., 0., 
                                                                   [0.05]*1.0e-11 , [0.9], 30., 10., collect(0.0 : 0.1 : 0.0))
                                                       ,use_experiment=false);
########################################################################################################
                           
                           
                           
# e_fac study .... 3 clear peaks ... down loping !!!! :D
# the same looping is appearing using EXP kinetics. 
#### when A or O is fast, middle peak starts to be the boundary one.... and makes a "inductance tail" 
# --- dva nebo tri obrazky
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, true, 0.0, 0.0, 0.0, 
                                              22.2106, 22.098, 21.5851, 
                                              -0.1497806, 0.0, 0.029967, 
                                              9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.2 : 0.5))
                                  ,use_experiment=false);


                                  
                                  
###########################################################################################
##########################        USED IN ARTICLE      ####################################
###########################################################################################
###########################################################################################
###########################################################################################

# par_stud_e_fac_LF
new_GAS_LoMA_shared = ysz_fitting.simple_run(
                                                       ysz_fitting.EIS_simulation(800, [80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                       #ysz_fitting.CV_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                       
                                                           pyplot=1,  
                                                            prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                                             prms_values=(0.0, 0.0, 0.0, 22.5, collect(22.9 : 0.1 : 22.9), 27.69994335225577, 0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, collect(0.0 : 0.05 : 0.2))
                                                            ,use_experiment=false);
                                                            
# par_stud_e_fac_HF
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                          prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"], 
                                          prms_values=(1.0, 1.0, 1.0, 21.6106, 27.898, 20.4851, -0.497806, -0.234529, 0.229967, 9.3e-11, 0.85, 3.04504, 0.05, [0.0, 0.1, 0.2, 0.3, 0.4])
                                          ,use_experiment=false);

                                          
# par_study_sha_sep_2_peaks >>> LaTeX tab. parameters: id 0.1
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                         prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                         prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                     22.0, 22., 21.5, 
                                                     0., 0.0, 0., 
                                                     9.3e-11, 0.5, 300., 100., collect(0.0 : 0.1 : 0.0))
                                         ,use_experiment=false);

# par_study_sha_sep_3_peaks >>> LaTeX tab. parameters: id 0.2
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                         prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                         prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                     22.0, 22., 21.5, 
                                                     0.15, 0.0, -0.629967, 
                                                     9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                         ,use_experiment=false);

# VERY interesting !!!
# separare_vacancy interchanged the order of reactions <-> peaks bijections >>> LaTeX tab.: id 0.2
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                            22.0, [22., 29.], 21.5, 
                                                            0.15, 0.0, -0.629967, 
                                                            9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false);

# par_study_DGR >>> LaTeX tab. parameters: id 0.1                                         
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                            22.0, 22., 21.5, 
                                                            0., collect(-0.3 : 0.1 : 0.1), 0., 
                                                            9.3e-11, 0.5, 300., 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false);
                                                
# HOT candidate !!! corrected_hand_guess -> expanded to 3 pO2 ..... CV rubbish
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true);
