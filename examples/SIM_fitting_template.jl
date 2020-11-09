#
#   TODO: - Some instructions for this
#         - finish the Snehurka control
#
#

# # # # # # function get_snehurka_control()
# # # # # #   #######################################################  
# # # # # #   # preparing bash output ###############################
# # # # # #   
# # # # # #   ### if true, no script is called! Just direclty run_par_study_script_wrap()
# # # # # #   direct_bool = true
# # # # # #   
# # # # # #             SIM_fitting_mode = true    #!#!#!#!#!#!#!#!#!#!
# # # # # #   
# # # # # #   bash_command = "sbatch"
# # # # # #   #bash_command = "echo"
# # # # # #   #bash_command = "julia"
# # # # # #   
# # # # # #   #mode = "test_one_prms"
# # # # # #   #mode = "only_print"
# # # # # #   mode = "go"
# # # # # #   
# # # # # #   express3_bool = true
# # # # # # 
# # # # # # 
# # # # # # 
# # # # # #   ##### TODO!!!! tohle bych mohl udelat taky obecne, at skrz jeden skriptovaci soubor muze projit vse
# # # # # #   
# # # # # #   if SIM_fitting_mode
# # # # # #     if express3_bool
# # # # # #       run_file_name = "../snehurka/run_EX3_ysz_fitting_SIM_fitting.jl"
# # # # # #     else
# # # # # #       run_file_name = "../snehurka/run_ysz_fitting_SIM_fitting.jl"
# # # # # #     end
# # # # # #   else
# # # # # #     if express3_bool
# # # # # #       run_file_name = "../snehurka/run_EX3_ysz_fitting_par_study-prms-.jl"
# # # # # #     else
# # # # # #       run_file_name = "../snehurka/run_ysz_fitting_par_study-prms-.jl"
# # # # # #     end  
# # # # # #   end
# # # # # #   
# # # # # #   return snehurka_control
# # # # # # end


function assemble_meta_SIM_fitting_snehurka(;only_return_SIM_fitting=false)  
  #######################################################
  #######################################################
  #######################################################
  ########### SIM_fitting definition ######################
 
 
  TC = 800
  pO2 = [40, 60]
  bias = 0.0

  data_set = "OLD_MONO_100"
  simulations = ["EIS"]
  fitness_factors = [1.0]
  
  physical_model_name = "ysz_model_GAS_LoMA_shared"
  
  #####
  
  prms_names=["separate_vacancy",
              "A.exp", "R.exp", "O.exp", 
              "A.r", "R.r", "O.r",              
              "A.DG", "R.DG", "O.DG",     
              "conductivity", "nu",      "OC", "ms_par", "e_fac"  ]
    
                                           
   prms_lists = (
      true,
      0.0, 0.0, 0.0,
      # rX
      collect(22. : 10.0 : 23.3),  
      collect(22. : 10.0 : 23.3),  
      collect(22. : 10.0 : 23.3), 
      # DGX
      collect(-0.3 : 10.20 : 0.3), 
      collect(-0.3 : 10.15 : 0.3),
      collect(-0.3 : 10.15 : 0.3),

      # hint: TC = (700, 750, 800, 850)  => DD = ( ??, 2.97, 7.27, 12.3)e-11 for "MONO_110"
      # hint: TC = (700, 750, 800, 850)  => DD = ( ??, 2.97, 9.3, 12.3)e-11 for "OLD_MONO_100"
      # hint: conductivity OLD_MONO_100 -> TC = [700, 750, 800, 850] = [1.02, 2.07,  3.72, 5.85]
      # hint: conductivity OLD_MONO_100 -> TC = [700, 750, 800, 850] = [1.015, 2.06,  3.63, 5.93]   # new version :) 
      1.02,
      0.65,

      20.0,
      20.0,
      0.0
    )
  
  
  mask          =(0,
                  0, 0, 0,
                  1, 1, 1,
                  1, 1, 1,
                  0, 1,       1, 1, 0)
  lower_bounds=(0.0, 
                0.0, 0.0, 0.0,
                15.5, 15.9, 15.7,       
                -0.8, -0.8, -0.8,
                5.71, 0.7,     0.0, 0.1, 0.0)
  upper_bounds=(1.0,
                1.0, 1.0, 1.0,
                27.5, 27.9, 27.7,              
                0.8, 0.8, 0.8,
                5.99, 0.94,      100, 100, 2.00)
                

  scripted_tuple =(1,
                  1, 1, 1,
                  1, 1, 1,       
                  1, 1, 1,
                  1, 1,           1, 1, 1)
                  
  
  #######################################################
  
  name = "cosi_kdesi"
  
  #######################################################
  #############  SIM_fitting construction ###############
  #######################################################
  #######################################################
  
  SIM_fitting = ysz_fitting.build_SIM_fitting(
                                    TC=TC,
                                    pO2=pO2,
                                    bias=bias,
                                    data_set=data_set,
                                    simulations=simulations,
                                    fitness_factors=fitness_factors,
                                    physical_model_name=physical_model_name,
                                    #
                                    prms_names=prms_names,
                                    #x0=output_prms_lists,
                                    mask=mask,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds,
                                    #
                                    print_to_file=false,
                                    save_dir="../data/SIM_fitting/temp/"*name*"/", 
                                    file_name="SIM_fitting_default.txt",
                                    #
                                    bboptimize_bool=false, 
                                    iteration_count=1000,
                                    )
                                    
  pyplot = false
  plot_each_x_th = 20
  print_only_result = true
  
  # 
  if only_return_SIM_fitting
    return SIM_fitting
  else
    return ysz_fitting.meta_run_par_study(only_return_SIM_fitting=false,
                            prms_lists=prms_lists,
                            pyplot=pyplot,
                            plot_each_x_th=plot_each_x_th,
                            print_only_result=print_only_result,
                            SIM_fitting=SIM_fitting,
                            scripted_tuple=scripted_tuple,
                          
                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = false,
  
                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!
                  
                            #bash_command = "sbatch",
                            #bash_command = "echo",
                            bash_command = "julia",
                            
                            #mode = "test_one_prms",
                            #mode = "only_print",
                            mode = "go",
                            
                            express3_bool = true
                            ) 
  
  
    return SIM_fitting
  end
end




function assemble_meta_SIM_fitting_no_beta_S(;only_return_SIM_fitting=false)  
  #######################################################
  #######################################################
  #######################################################
  ########### SIM_fitting definition ######################
 
 
  TC = 800
  pO2 = [40, 60]
  bias = 0.0

  data_set = "OLD_MONO_100"
  simulations = ["EIS"]
  fitness_factors = [1.0]
  
  physical_model_name = "ysz_model_GAS_LoMA_shared"
  
  #####
  
  prms_names=["separate_vacancy",
              "A.exp", "R.exp", "O.exp", 
              "A.r", "R.r", "O.r",              
              "A.DG", "R.DG", "O.DG",     
              "conductivity", "nu",      "OC", "ms_par", "e_fac"  ]
 
  prms_lists=(true, 1.0, 1.0, 1.0, 21.96912506564684, 21.569512049276458, 20.93977987957153, 0.1241112223034025, -0.10660629726868907, -0.006580383278506099, 9.3e-11, 0.85, 5.848806563958082, 8.27404376722011, 0.0)
  
  prms_lists=(1, 0.0, 0.0, 0.0, 21.9464, 21.5159, 21.4928, 0.0032943, 0.0277383, -0.0838499, 1.02, 0.811907, 10.472, 6.00553, 0.0718344)
  # for 800 >> 40, 60 >> 0.0
  prms_lists=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 3.62, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)

  prms_lists=(true, 0.0, 0.0, 0.0, 
                                        22.86783400557656, 22.1237966317968, 22.40693968907697, 
                                        -0.20520337660564685, -0.173209183842112, 0.16770403513536281, 
                                        3.62, 0.9, 5.08895790653341, 3.269908550672186, 0.06182027814202923)
    
                                           
#     prms_lists = (
#       true,
#       0.0, 0.0, 0.0,
#       # rX
#       collect(21.5 : 1.0 : 21.5),
#       collect(21.5 : 1.0 : 21.5),
#       collect(21.5 : 1.0 : 21.5),
#       # DGX
#       collect(-0.25 : 0.05 : 0.25),
#       collect(-0.25 : 0.05 : 0.25),
#       collect(-0.25 : 0.05 : 0.25),
# 
#       # hint: TC = (700, 750, 800, 850)  => DD = ( ??, 2.97, 7.27, 12.3)e-11 for "MONO_110"
#       # hint: TC = (700, 750, 800, 850)  => DD = ( ??, 2.97, 9.3, 12.3)e-11 for "OLD_MONO_100"
#       # hint: conductivity OLD_MONO_100 -> TC = [700, 750, 800, 850] = [1.02, 2.07,  3.72, 5.85]
#       # hint: conductivity OLD_MONO_100 -> TC = [700, 750, 800, 850] = [1.015, 2.06,  3.63, 5.93]   # new version :) 
#       1.02,
#       0.65,
# 
#       20.0,
#       20.0,
#       0.0
#     )
  
  
  mask          =(0,                
                  0, 0, 0,
                  1, 1, 1,
                  1, 1, 1,
                  0, 1,       1, 1, 1)
                  
  lower_bounds=(0.0,                
                0.0, 0.0, 0.0,
                15.5, 15.9, 15.7,       
                -0.8, -0.8, -0.8,
                [1]*1.0e-13, 0.85,     0.01, 0.01, 0.0)
                
  upper_bounds=(1.0,
                1.0, 1.0, 1.0,
                27.5, 27.9, 27.7,              
                0.8, 0.8, 0.8,
                [1]*1.0e-8, 0.94,      Inf, Inf, 0.5)
                

  scripted_tuple =(1,
                  1, 1, 1,
                  1, 1, 1,       
                  1, 1, 1,
                  1, 1,           1, 1, 1)
                  
  
  #######################################################
  
  name = "cosi_kdesi"
  
  #######################################################
  #############  SIM_fitting construction ###############
  #######################################################
  #######################################################
  
  SIM_fitting = ysz_fitting.build_SIM_fitting(
                                    TC=TC,
                                    pO2=pO2,
                                    bias=bias,
                                    data_set=data_set,
                                    simulations=simulations,
                                    fitness_factors=fitness_factors,
                                    physical_model_name=physical_model_name,
                                    #
                                    prms_names=prms_names,
                                    #x0=output_prms_lists,
                                    mask=mask,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds,
                                    #
                                    print_to_file=false,
                                    save_dir="../data/SIM_fitting/temp/"*name*"/", 
                                    file_name="SIM_fitting_default.txt",
                                    #
                                    bboptimize_bool=false, 
                                    iteration_count=100,
                                    )
                                    
  pyplot = true
  plot_each_x_th = 20
  print_only_result = true
  
  # 
  if only_return_SIM_fitting
    return SIM_fitting
  else
    return ysz_fitting.meta_run_par_study(only_return_SIM_fitting=false,
                            prms_lists=prms_lists,
                            pyplot=pyplot,
                            plot_each_x_th=plot_each_x_th,
                            print_only_result=print_only_result,
                            SIM_fitting=SIM_fitting,
                            scripted_tuple=scripted_tuple,
                          
                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = true,
  
                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!
                  
                            bash_command = "sbatch",
                            #bash_command = "echo",
                            #bash_command = "julia",
                            
                            #mode = "test_one_prms",
                            #mode = "only_print",
                            mode = "go",
                            
                            express3_bool = true
                            ) 
  
  
    return SIM_fitting
  end
end






function assemble_meta_SIM_fitting_TEMPERATURE__OLD(;only_return_SIM_fitting=false)  
  #######################################################
  #######################################################
  #######################################################
  ########### SIM_fitting definition ######################
 
 
  TC = [700, 750, 800, 850]
  pO2 = [40, 60]
  bias = 0.0

  data_set = "OLD_MONO_100"
  simulations = ["EIS"]
  fitness_factors = [1.0]
  
  physical_model_name = "ysz_model_GAS_LoMA_Temperature"
  
  #####
  
  prms_names=["separate_vacancy", "e_fac", "A.exp", "R.exp", "O.exp",
                   "A.r_B", "A.r_C",        "A.DG_B", "A.DG_C",
                   "R.r_B", "R.r_C",        "R.DG_B", "R.DG_C",
                   "O.r_B", "O.r_C",        "O.DG_B", "O.DG_C",
                   "nu_B", "nu_C",     "CO_B", "CO_C",     "COmm_B", "COmm_C"]
 
 
  prms_lists=(1, 0.27, 0.0, 0.0, 0.0, 
                   collect(-0. : 10.4 : 0.2), 21.92,            collect(-0. : 10.1 : 0.1), 0.055,
                   collect(-0. : 10.4 : 0.2), 21.21,            collect(-0. : 10.1 : 0.1), 0.048,
                   collect(-0. : 10.4 : 0.2), 21.809,           collect(-0. : 10.1 : 0.1), -0.04,
                   -0.2, 0.76,       5.0, 37.97,      5.0, 5.56*0.15)
  
  
  mask          =(0, 0, 0, 0, 0,               
                  1, 0,           1, 0,
                  1, 0,           1, 0,
                  1, 0,           1, 0,                  
                  1, 0,     1, 0,     1, 0)
                  
  lower_bounds=(0.0, 0.0, 0.0, 0.0, 0.0,
                -3, 20.5,         -0.5, -0.8, 
                -3, 20.5,         -0.5, -0.8, 
                -3, 20.5,         -0.5, -0.8, 
                -0.3, 0.01,     -10, 0.05,    -10, 0.05)      

                  
  upper_bounds=(1.0, 0.5, 1.0, 1.0, 1.0,
                3, 26.5,         0.5, 0.8, 
                3, 26.5,         0.5, 0.8, 
                3, 26.5,         0.5, 0.8, 
                0.3, 0.95,     10, 100.0,    10, 100.0)
                

  scripted_tuple =(1, 1, 1, 1, 1,
                  1, 1,           1, 1,
                  1, 1,           1, 1,       
                  1, 1,           1, 1,
                  1, 1,       1, 1,     1, 1)
                  
  
  #######################################################
  
  name = "cosi_kdesi"
  
  #######################################################
  #############  SIM_fitting construction ###############
  #######################################################
  #######################################################
  
  SIM_fitting = ysz_fitting.build_SIM_fitting(
                                    TC=TC,
                                    pO2=pO2,
                                    bias=bias,
                                    data_set=data_set,
                                    simulations=simulations,
                                    fitness_factors=fitness_factors,
                                    physical_model_name=physical_model_name,
                                    #
                                    prms_names=prms_names,
                                    #x0=output_prms_lists,
                                    mask=mask,
                                    lower_bounds=lower_bounds,
                                    upper_bounds=upper_bounds,
                                    #
                                    print_to_file=false,
                                    save_dir="../data/SIM_fitting/temp/"*name*"/", 
                                    file_name="SIM_fitting_default.txt",
                                    #
                                    bboptimize_bool=false, 
                                    iteration_count=1000,
                                    )
                                    
  pyplot = true
  plot_each_x_th = 100
  print_only_result = true
  
  # 
  if only_return_SIM_fitting
    return SIM_fitting
  else
    return ysz_fitting.meta_run_par_study(only_return_SIM_fitting=false,
                            prms_lists=prms_lists,
                            pyplot=pyplot,
                            plot_each_x_th=plot_each_x_th,
                            print_only_result=print_only_result,
                            SIM_fitting=SIM_fitting,
                            scripted_tuple=scripted_tuple,
                          
                              ### if true, no script is called! Just direclty run_par_study_script_wrap()
                              direct_bool = true,
  
                            SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!
                  
                            #bash_command = "sbatch",
                            #bash_command = "echo",
                            bash_command = "julia",
                            
                            #mode = "test_one_prms",
                            #mode = "only_print",
                            mode = "go",
                            
                            express3_bool = true
                            ) 
  
  
    return SIM_fitting
  end
end
  
  



# function assemble_meta_SIM_fitting_with_beta_and_S(;only_return_SIM_fitting=false)  
#   #######################################################
#   #######################################################
#   #######################################################
#   ########### par_study definition ######################
#   
#   name = "cosi_kdesi"
#   
#   
#   
#   
#   prms_names=["separate_vacancy",
#               "A.beta", "R.beta", "O.beta",
#               "A.S", "R.S", "O.S",
#               "A.exp", "R.exp", "O.exp", 
#               "A.r", "R.r", "O.r",              
#               "A.DG", "R.DG", "O.DG",     
#               "conductivity", "nu",      "OC", "ms_par", "e_fac"  ]
#  
#   prms_lists=(true, 1.0, 1.0, 1.0, 21.96912506564684, 21.569512049276458, 20.93977987957153, 0.1241112223034025, -0.10660629726868907, -0.006580383278506099, 9.3e-11, 0.85, 5.848806563958082, 8.27404376722011, 0.0)
#   
#   prms_lists=(1, 
#               0.5, 0.5, 0.5,
#               0, 0, 0,
#               0.0, 0.0, 0.0, 21.9464, 21.5159, 21.4928, 0.0032943, 0.0277383, -0.0838499, 1.02, 0.811907, 10.472, 6.00553, 0.0718344)
#   
#   prms_lists=(1, 0.0, 0.0, 0.0, 22.4569, 22.1153, 23.7591, -0.0321725, -0.0408594, 0.600798, 3.58, 0.23498, 40.5525, 20.8275, 0.414237)
#  
# #     prms_lists = (
# #       true,
# #       0.0, 0.0, 0.0,
# #       # rX
# #       collect(21.5 : 1.0 : 21.5),
# #       collect(21.5 : 1.0 : 21.5),
# #       collect(21.5 : 1.0 : 21.5),
# #       # DGX
# #       collect(-0.25 : 0.05 : 0.25),
# #       collect(-0.25 : 0.05 : 0.25),
# #       collect(-0.25 : 0.05 : 0.25),
# # 
# #       # hint: TC = (700, 750, 800, 850)  => DD = ( ??, 2.97, 7.27, 12.3)e-11 for "MONO_110"
# #       # hint: TC = (700, 750, 800, 850)  => DD = ( ??, 2.97, 9.3, 12.3)e-11 for "OLD_MONO_100"
# #       # hint: conductivity OLD_MONO_100 -> TC = [700, 750, 800, 850] = [1.02, 2.07,  3.58, 5.85]
# # 
# #       1.02,
# #       0.65,
# # 
# #       20.0,
# #       20.0,
# #       0.0
# #     )
#   
#   
#   mask          =(0,
#                   1,1,1,
#                   1,1,1,
#   
#                   0, 0, 0,
#                   0, 0, 0,
#                   0, 0, 0,
#                   0, 0,       0, 0, 0)
#   lower_bounds=(0.0, 
#                 0.0, 0.0, 0.0,
#                 -3, -3, -3,
#                
#                 0.0, 0.0, 0.0,
#                 15.5, 15.9, 15.7,       
#                 -0.8, -0.8, -0.8,
#                 [1]*1.0e-13, 0.01,     0.0, 0.1, 0.0)
#   upper_bounds=(1.0,
#                 1, 1, 1,
#                 3, 3, 3,
#   
#                 1.0, 1.0, 1.0,
#                 27.5, 27.9, 27.7,              
#                 0.8, 0.8, 0.8,
#                 [1]*1.0e-8, 0.99,      Inf, Inf, 5.0)
#                 
# 
#   scripted_tuple =(1,
#                   1,1,1,
#                   1,1,1,
#     
#                   1, 1, 1,
#                   1, 1, 1,       
#                   1, 1, 1,
#                   1, 1,           1, 1, 1)
#  
#   TC = 700
#   pO2 = [40]
#   bias = 0.0
# 
#   data_set = "OLD_MONO_100"
#   simulations = ["CV"]
#               
#   physical_model_name = "ysz_model_GAS_LoMA_shared"
#   
#   #######################################################
#   ###############  SIM_fitting definition ###############
#   #######################################################
#   #######################################################
#   
#   SIM_fitting = build_SIM_fitting(
#                                     TC=TC,
#                                     pO2=pO2,
#                                     bias=bias,
#                                     data_set=data_set,
#                                     simulations=simulations,
#                                     physical_model_name=physical_model_name,
#                                     #
#                                     prms_names=prms_names,
#                                     #x0=output_prms_lists,
#                                     mask=mask,
#                                     lower_bounds=lower_bounds,
#                                     upper_bounds=upper_bounds,
#                                     #
#                                     print_to_file=false,
#                                     save_dir="../data/EEC/temp/log/", 
#                                     file_name="SIM_fitting_default.txt",
#                                     #
#                                     bboptimize_bool=false, 
#                                     iteration_count=100,
#                                     )
# 
#   pyplot = true
#   plot_each_x_th = 20
#   print_only_result = true
#   
#   # 
#   if only_return_SiM_fitting
#     return SIM_fitting
#   else
#     return meta_run_par_study(only_return_SIM_fitting=false,
#                             pyplot=pyplot,
#                             plot_each_x_th=plot_each_x_th,
#                             print_only_result=print_only_result,
#                             SIM_fitting=SIM_fitting,
#                             scripted_tuple=scripted_tuple,
#                           
#                               ### if true, no script is called! Just direclty run_par_study_script_wrap()
#                               direct_bool = true,
#   
#                             SIM_fitting_mode = true,    #!#!#!#!#!#!#!#!#!#!
#                   
#                             bash_command = "sbatch",
#                             #bash_command = "echo",
#                             #bash_command = "julia",
#                             
#                             #mode = "test_one_prms",
#                             #mode = "only_print",
#                             mode = "go",
#                             
#                             express3_bool = true
#                             ) 
#   end
# end
