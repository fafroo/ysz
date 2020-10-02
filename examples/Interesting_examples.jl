


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
# 4 EIS ... not capturing the non-monotony ... even with extended search area
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC"     ], 
                                   prms_values=(0.0, 0.0, 0.0, 23.536793805295844, 26.44428245339161, 21.308492340939917, -0.5500221642340698, 0.029112401256050482, -0.020889037864824146, 9.3e-11, 0.85, 2.0593360762011694)
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
