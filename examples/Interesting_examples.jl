


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

                                          
################ CORRECTED ######################## ysz_model_GAS_LoMA ############################
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

# e_fac ... R_ohm static ... example
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["e_fac", "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC"     ], 
                                   prms_values=([0.0, 0.1, 0.2, 0.3, 0.4], 1.0, 1.0, 1.0, 20.8106, 27.898, 20.4851, -0.097806, -0.434529, 0.229967, 9.3e-11, 0.85, 3.04504)
                                   ,use_experiment=false);                                   
                                   
# e_fac ... R_ohm varying !!!! 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["e_fac", "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC"     ], 
                                   prms_values=([0.0, 0.1, 0.2, 0.3, 0.4], 1.0, 1.0, 1.0, 23.8106, 27.898, 20.4851, -0.097806, -0.434529, 0.229967, 9.3e-11, 0.85, 3.04504)
                                   ,use_experiment=false);
# e_fac ... HF right turn
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["e_fac", "expA", "expR", "expO", "A0", "R0", "K0", "DGA", "DGR", "DGO",      "DD", "nu", "OC"     ], 
                                   prms_values=([0.1, 0.2, 0.3, 0.4], 1.0, 1.0, 1.0, 21.8106, 27.898, 20.4851, -0.697806, -0.434529, 0.229967, 9.3e-11, 0.85, 3.04504)
                                   ,use_experiment=false);
