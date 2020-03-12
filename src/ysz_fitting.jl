module ysz_fitting
#######################
####### TODO ##########
# [x] better ramp ... starting directly from steadystate :)
# [x] general global search using projection to each variable
# [x] compute EIS exacly on checknodes and therefore remove plenty of "EIS_apply_checknodes"
# [x] put appropriate and finished stuff into "CV_fitting_supporting_stuff"
# [x] spoustet run_new() v ruznych procedurach stejnou funkci "EIS_default_run_new( ... )"
# [o] implement LM algorithm
# [x] sjednotit, co znamena pO2, jeslti jsou to procenta a jak se prenasi do simulace!  a taky T .. Celsia a Kelvina
# [x] srovnat vodivost elektrolytu s experimentem CV i EIS naraz
# [x] vymyslet novy relevantni vektor parametru
# [?] aplikovat masku pro fitting pro scan_2D_recursive
# [ ] snehurka, velke objemy dat, maska a metafile
# ---[x] save_dir preposilany parametrem a odlisit scripted_dir
# ---[x] sneh - zadavani i jinych parametru nez prms
# ---[x] sneh - pouzivat jen jeden soubor pro spouteni sbatch ... v hlavicce #!(..)/julia
# ---[ ] mozna pouzit precompile
# ---[ ] nechat meta_run_par_study() vytvorit skriptovaci soubor (dle poctu parametru, module, pocet experimentu)
# ------[ ] spis bude lepsi skriptovaci soubor nechat stejny a prohanet pres nej jen ARGS bez cisel
# [!] postarat se o modularitu ysz_model_COSI, at se nemusi vytvaret znovu "experiments_COSI" -> JUERGEN
# ---[!] zadavani modulu by melo byt v hlavicce meta_run_par_study()
# ------[!] mozna by se include mel vykonavat na vyssi urovni a do ysz_experiments preposilat uz hotovy modul
# [x] vymyslet lepe zadavani parametru pro ruzne modely s ruznymi parametry 
# [x] opravdu promyslet to objektove programovatni ... CV_sim, EIS_sim ... prms_lists ...
# ---[x] pridat treba dalsi experiment s dalsimi daty, vuci kterym se da srovnavat (kapacitance)
# ---[x] udelat tridu par_study
# [x] PAR_STUDY vyhodnoceni ... soubory do slozky dane studie .. automatizovane
# ---[ ] pri zobrazovani dat udelat clustery dle chyby/prms a pak zobrazit (prms ci Nyquist) 1 reprezentanta z kazde
# ---[ ] automatizovat ukladani souboru vysledku par_study.vyhodnoceni()
# [ ] do experiments.jl pridat obecne zaznamenavani promennych od final_plot
# [x] get rig of shared_prms and shared_add_prms !!!
# [ ] snehurka by rada ./zabal.sh a pocitac zase ./rozbal_par_study.sh
# [!] snehurkove fitovani by slo zrychlit, kdyz bych skriptoval EQ parametry a job by menil jen kineticke? ... chrm ...
# ---[ ] nejak si preposilat steadystate?! 
# ------[x] zacina se ze steady_statu, ktery se pocita rychle
# [ ] simple_run by mohl vracet par_study
# ---[ ] par study by pak mohla mit funcki "uloz me"
# [ ] par_study_struct by mohla nest velikosti poli
# [ ] promyslet velikost Floatu pri importu velkeho mnozstvi dat
# [x] funkci pro par_study, aby umela zmenit sve info (redukovat pO2_list a tak) -> funkce par_study_filter()
# [ ] par_study_plot_the_best() by mela vyhodit jeden figure s dvema subploty... asi ... 
# [ ] prms_names a prms_values dat do jedne tridy
# ---[ ] vlastne jakoby nevim, jestli je to dobry napad?
# [ ] automaticke vybirani vhodne scripted_tuple
#
# ##### interpolace ######
# [o] vizualizace trendu chyby mezi interpolanty
# ---[o] zobrazovat soucty odchylek
#
# #### od Affra ####
# [ ] fitness funkce s maximovou metrikou
#
# 
# #### Fitting process
# [!] prozkoumat experimentalni data (a udelat prislusne procedury)
# [o] non-gas-ads-model
# ---[ ] find best add_prms (from cap. fitting)
# ---[ ] find best prms for both EIS and CV
# ------[x] fing best prms for EIS
# ------[ ] fing best prms for CV
# ------[ ] the TRICK with SA ! .. R0 := R0/SA (exp .... - exp ...)
# ---[o] try perturbate add prms by par_study
# [ ] gas-exp-model
# [ ] gas-LoMA-model
# ---[ ] try to fit CV and EIS at once!
# ------[ ] compute a lots of EIS, compute some CVs and interpolate
#######################


using Printf
using PyPlot
#using Plots
using DataFrames
using LeastSquaresOptim

import Base.string

########  import of model_file  ######
#TODO !!!!


######################################

#include("../examples/ysz_experiments.jl")


include("../src/general_supporting_stuff.jl")
include("../src/import_experimental_data.jl")
include("../src/general_simulation.jl")
include("../src/CV_simulation.jl")
include("../src/EIS_simulation.jl")
include("../src/par_study.jl")
#include("../src/interpolation_fitting.jl")






function get_fitted_all_prms()
  # for which model???
  prms_names=["A0", "R0", "K0", "SA", "SR", "SO", "DGA", "DGR", "DGO", "betaA", "betaR", "betaO", "DD"]
  #prms_values=[19.7, 19.7, 18.6,    1, 1, 1,    0.7, -0.8, -0.3,      0.5, 0.5, 0.5,    5.35e-13]  # fitted to EIS 800, 100, 0.0
    
  return prms_names, prms_values
end




function simple_run(SIM_list::Array{abstract_simulation}; pyplot=0, use_experiment=true, prms_values=[], prms_names=[], 
                        test=false)
  # here starts the true body
  if test
    test_result = 0
  end
  
  for SIM in SIM_list
  
    function recursive_simple_run_call(output_prms, plot_names, plot_values, active_idx)
      if active_idx > size(prms_names,1)

        ###################################################################
        # for each combination of output_prms from prms_lists perform this #######
        prms_names_in=[]
        prms_values_in=[]
        
        append!(prms_names_in, prms_names)
        append!(prms_values_in, output_prms)
      
        
        if check_equal_size(prms_names, prms_values)  
          SIM_raw = typical_run_simulation(SIM, prms_names_in, prms_values_in, pyplot)        
        end    
        if size(plot_names,1) < 1
          plot_prms_string = ""
        else
          plot_prms_string = " $(string(plot_names)) = $(plot_values)"
        end
        SIM_sim = apply_checknodes(SIM, SIM_raw, checknodes)
        if pyplot > 0
            typical_plot_sim(SIM, SIM_sim, plot_prms_string)
        end  
        if use_experiment
          if test
            test_result+=fitnessFunction(SIM, SIM_exp, SIM_sim)
          else
            fitness_error_report(SIM, plot_prms_string, SIM_exp, SIM_sim)
          end
        end
        return
        ###################################################################
        ###################################################################
      
      end
      if size(prms_values[active_idx],1)>1
        for i in prms_values[active_idx]
          recursive_simple_run_call(
            push!(deepcopy(output_prms),i),
            push!(deepcopy(plot_names),prms_names[active_idx]),
            append!(deepcopy(plot_values),i),
            active_idx + 1)
        end
      else
        recursive_simple_run_call(
          push!(deepcopy(output_prms),prms_values[active_idx][1]),
          plot_names,
          plot_values,
          active_idx + 1)
      end
    end

    # here the true body continues
    checknodes =  get_shared_checknodes(SIM)
    if use_experiment
      SIM_exp = apply_checknodes(SIM, import_data_to_DataFrame(SIM), checknodes)
      if (pyplot > 0)
        typical_plot_exp(SIM, SIM_exp)
      end
    end
    recursive_simple_run_call([], Array{String}(undef,(0)), Array{Float64}(undef,(0)), 1)
  end
  
  if test
    return test_result
  else
    return
  end
end


# useful wrap
function simple_run(;TC, pO2, bias=0.0, simulations=[], pyplot=0, use_experiment=true, prms_values=[], prms_names=[], 
                         test=false)
    simple_run(get_SIM_list_rectangle(TC, pO2, bias, simulations); pyplot=pyplot, use_experiment=use_experiment, prms_values=prms_values, prms_names=prms_names, 
                        test=false)
end


###########################################################
###########################################################
#### Fitness_function_interpolation #######################
###########################################################
###########################################################

function findin(a::Array, item)
  findall(x->x==item, a)
end

function remove!(a::Array, item)
    deleteat!(a, findin(a, item))
end

function par_study_get_sim(this_par_study::par_study_struct, SIM_idx, prms_indicies::Union{Array, Tuple})
  this_par_study.SIMs_data[SIM_idx, prms_indicies...  ]
end

function in_interval_filter(array::Array, a, b)
  res = []
  for item in array
    if is_between(item, a, b)
      append!(res, [item])
    end
  end
  return res
end

function par_study_display_fitness_profile(tested_prm_name, int_a, int_b)
  
  function for_each_non_tested_indicies_in_prms_lists(this_par_study, perform_generic)
    function recursive_call(output_set, active_idx)
      if active_idx > size(this_par_study.prms_lists,1)
        perform_generic(output_set)
      else
        if this_par_study.prms_names[active_idx]==tested_prm_name
          recursive_call(
            push!(
              deepcopy(output_set), 
              this_par_study.prms_lists[active_idx]
            ), 
            active_idx + 1
          )
        else
          for parameter in this_par_study.prms_lists[active_idx]
            recursive_call(push!(deepcopy(output_set), parameter), active_idx + 1)
          end
        end
      end
    end
    recursive_call([],1)
    return
  end

  function display_fitness_profile(output_prms)
    println("ted")
    test_par_study = ysz_fitting.par_study_get_subset_of_info(act_par_study, prms_names=act_par_study.prms_names, prms_values=output_prms);
    ysz_fitting.par_study_import_data!(test_par_study, verbose=0);
    
    idx_array = output_prms
    ok_bool = true
    for (SIM_idx, SIM) in enumerate(test_par_study.SIM_list)
      
      idx_array .= 1
      SIM_ref = par_study_get_sim(test_par_study, SIM_idx, idx_array)
      if size(SIM_ref,1) == 0
        ok_bool=false
        println("F")
        return
      end
      
      # preparing trend_tuples
      trend_tuples = initialize_trend_tuples(SIM, SIM_ref)
      
      #typical_plot_exp(SIM, SIM_ref)
      for (prm_idx, value) in enumerate(test_par_study.prms_lists[tested_prm_idx])
        idx_array[tested_prm_idx] = prm_idx
        SIM_test = par_study_get_sim(test_par_study, SIM_idx, idx_array)
        if size(SIM_test, 1) == 0
          ok_bool=false
          println("F")
          return
        end
        #typical_plot_sim(SIM, SIM_test)
        
        trend_tuple = get_trend_tuple(SIM, SIM_ref, SIM_test)
        push!(trend_tuples, [value, trend_tuple...])
      end
      for i in 1:(size(trend_tuples,2)-1)
        trend_tuples[!, Symbol(string(i))] = trend_tuples[!, Symbol(string(i))]./(trend_tuples[!, Symbol(string(i))][end])
        
        figure(8)
        plot(trend_tuples.prm_value, trend_tuples[!, Symbol(string(i))], "x")  
        
        
        #SIM_trend_tuples_holder[SIM_idx], trend_tuples)
      end
      
    end
  end
  
  my_par_study = ysz_fitting.par_study_import_info_from_metafile("GAS_LoMA_$(tested_prm_name)_test")
    
  tested_prm_idx = findin(my_par_study.prms_names, tested_prm_name)
  if size(tested_prm_idx,1) < 1
    println("ERROR: par_study_display_fitness_profile: tested_prm_name not in prms_names ($tested_prm_name in $(my_par_study.prms_names))")
    throw(Exception)
  elseif size(tested_prm_idx, 1) > 1
    println("ERROR: par_study_display_fitness_profile: tested_prm_name not UNIQUE in prms_names ($tested_prm_name in $(my_par_study.prms_names))")
    throw(Exception)
  else
    tested_prm_idx = tested_prm_idx[1]
  end
  
  # TO DELETE !!!!
  act_par_study = ysz_fitting.par_study_get_subset_of_info(my_par_study, 
    simulations=["CV"], pO2=[60], 
    prms_names=[tested_prm_name], 
    prms_values=[in_interval_filter(my_par_study.prms_lists[tested_prm_idx], int_a, int_b)]);
  
  #test_par_study = my_par_study;
        
  SIM_trend_tuples_holder = Array{Any}(undef, size(act_par_study.SIM_list, 1))
  
  PyPlot.ioff()
  for_each_non_tested_indicies_in_prms_lists(act_par_study, display_fitness_profile)
  PyPlot.show()

end





###########################################################
###########################################################
#### Levenberg-Marquardt ##################################
###########################################################
###########################################################


# function LM_optimize(;EIS_opt_bool=false, CV_opt_bool=false, pyplot=false)
#     function prepare_prms(mask, x0, x)
#         prms = zeros(0)
#         xi = 1
#         for i in collect(1 : 1 :length(mask))
#             if convert(Bool,mask[i])
#                 append!(prms, x[xi])
#                 xi += 1
#             else
#                 append!(prms, x0[i])
#             end
#         end
#         return prms
#     end
# 
#    function to_optimize(x) 
#         #err = run_new(print_bool=false, fitting=true, voltametry=true, pyplot=false, voltrate=0.005, sample=8, bound=0.41, 
#         #    prms_in=x)
#         prms = prepare_prms(mask, x0, x)
#         print(" >> mask = ",mask)
#         print(" || prms = ",prms)
#         
#         prms_names_in=[]
#         prms_values_in=[]
#         
#         append!(prms_names_in, prms_names)
#         append!(prms_values_in, prms)
#       
#        
#         err = 0.0
#         err_string = ""
#         if CV_opt_bool
#             CV_sim = ysz_experiments.run_new(
#                         physical_model_name=physical_model_name,
#                         out_df_bool=true, voltammetry=true, sample=8, pyplot=false,
#                         T=TCtoT(TC), pO2=pO2tosim(pO2),
#                         prms_names_in=prms_names_in,
#                         prms_values_in=prms_values_in,
#             )
#             if pyplot
#                 figure(1)
#                 CV_plot(CV_sim)
#             end
#             
#             err +=CV_penalty_factor * CV_fitnessFunction(
#                 CV_apply_checknodes(CV_sim, CV_get_shared_checknodes()), 
#                 CV_exp
#             )
#             err_string = string(err_string," CV ")
#         end
#         if EIS_opt_bool
#             EIS_sim = ysz_experiments.run_new( 
#                     physical_model_name=physical_model_name,
#                     pyplot=false, EIS_IS=true, out_df_bool=true, bias=bias, omega_range=EIS_get_shared_omega_range(),
#                     dx_exp=-9,
#                     T=TCtoT(TC), pO2=pO2tosim(pO2),
#                     prms_names_in=prms_names_in,
#                     prms_values_in=prms_values_in,
#             )
# 
#             
#             if pyplot
#                 figure(2)
#                 Nyquist_plot(EIS_sim)
#             end
#             
#             err += EIS_fitnessFunction(EIS_sim, EIS_exp)
#             err_string = string(err_string," EIS ")
#         end
#         println(" || ",err_string,": err =", err)
#         return [err]
#     end
#     
# 
#     
#     TC=800
#     pO2 = 100
#     bias=0.0
#     physical_model_name="necum"
#     CV_penalty_factor = 10  # for fitnessFunction = factor*CV + EIS 
#     
#     
#     if pyplot
#         PyPlot.close()
#     end
#     
#     if CV_opt_bool
#         CV_exp = CV_apply_checknodes(
#             import_CVtoDataFrame(TC=TC, pO2=pO2),
#             CV_get_shared_checknodes()
#         )
#         if pyplot
#             figure(1)
#             CV_plot(CV_exp, "exp $(CV_experiment_legend(TC, pO2))")
#         end
#     end
#     
#     if EIS_opt_bool
#         EIS_exp = EIS_apply_checknodes(
#             import_EIStoDataFrame(TC=TC,pO2=pO2,bias=bias),
#             EIS_get_shared_checknodes()
#         )
# 
#         if pyplot
#             figure(2)
#             Nyquist_plot(EIS_exp, "exp $(EIS_experiment_legend(TC, pO2, bias))")
#         end
#     end
#             
#     #x0 = zeros(2)
#     #optimize(rosenbrock, x0, LevenbergMarquardt())
#     
#     prms_names=["A0", "R0", "K0", "SA", "SR", "SO", "DGA", "DGR", "DGO", "betaA", "betaR", "betaO", "DD", "nu", "nus", "ms_par"]
#     
#     lower_bounds      = [17, 17, 17,     -5, -5, -5,     -0.9, -0.9, -0.9,    0.0, 0.0, 0.0,     0.1e-13, 0.1, 0.1, 0.01 ]
#     upper_bounds      = [23, 23, 23,      1,  1,  1,      0.9,  0.9,  0.9,    1.0, 1.0, 1.0,     9.0e-13, 0.9, 0.9, 1.5 ]
#     
#     x0  = [19.7, 19.7, 18.6,    1, 1, 1,    0.7, -0.8, -0.3,      0.5, 0.5, 0.5,     5.35e-13,  0.85,  0.21, 0.05]
#     #x0 = [18.673485702096716, 18.598123074869402, 17.50858588747129, 1.0, 1.0, 1.0, 0.8666319406487695, -0.8342275189659124, -0.7728570698687104, 0.5, 0.5, 0.5, 5.35e-13] 
#     
#     
#     mask = [1, 1, 1,     0, 0, 0,     1,1,1,      0,0,0,   1, 1, 1, 1] # determining, which parametr should be fitted
#     
#     x0M = zeros(0)
#     lowM = zeros(0)
#     uppM = zeros(0)
#     for i in collect(1 : 1 : length(mask))
#         if convert(Bool,mask[i])
#             append!(x0M, x0[i])
#             append!(lowM, lower_bounds[i])
#             append!(uppM, upper_bounds[i])
#         end
#     end
#     
#    
#     
#     #to_optimize(x0M)
#     println(optimize(to_optimize, x0M, lower=lowM, upper=uppM, Δ=1000000, f_tol=1.0e-14, g_tol=1.0e-14, LevenbergMarquardt()))
#     
#     ####optimize(to_optimize, x0M, Dogleg())
#     return
# end


###########################################################
###########################################################
#### CLUSTER COMPUTATION ##################################
###########################################################
###########################################################



function meta_run_par_study()  
  function recursive_bash_call(output_prms_lists, active_idx)
    if active_idx > size(scripted_tuple,1)
      scripted_tuple_string = string(scripted_tuple)
      output_prms_lists_string = string(output_prms_lists)
      prms_names_string = string(prms_names)
      simulations_string = string(simulations)
      TC_string = string(TC)
      pO2_string = string(pO2)
      bias_string = string(bias) 
      
      if direct_bool
        return run_par_study_script_wrap(
                    output_prms_lists_string ,
                    save_dir ,
                    name ,
                    scripted_tuple_string ,
                    prms_names_string ,
                    TC_string ,
                    pO2_string ,
                    bias_string ,
                    simulations_string ,
                    mode ,
                    physical_model_name        
        )
      else
        return run(`
                    $bash_command  
                    $run_file_name 
                    $output_prms_lists_string 
                    $save_dir 
                    $name 
                    $scripted_tuple_string 
                    $prms_names_string 
                    $TC_string 
                    $(pO2_string) 
                    $bias_string 
                    $simulations_string 
                    $mode 
                    $physical_model_name 
        `)
      end
    else
      if scripted_tuple[active_idx] == 1
        for i in prms_lists[active_idx]
          recursive_bash_call(push!(deepcopy(output_prms_lists),i), active_idx + 1)
        end
      else
        recursive_bash_call(push!(deepcopy(output_prms_lists),prms_lists[active_idx]), active_idx + 1)
      end
    end
  end
  
  function consistency_check()
    if (size(prms_names,1) != size(scripted_tuple,1) || 
        size(prms_names,1) != size(prms_lists,1))
      return false
    end
    
    for i in 1:size(prms_lists,1)
      if size(prms_lists[i],1) < 1
        return false
      end
    end
    
    return true
  end
  
  #######################################################
  #######################################################
  #######################################################
  ########### par_study definition ######################
  
  name = "GAS_LoMA_ULTRA_level_2"
  
  prms_names = ["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"]
  # BE CAREFUL! data are saved to files only with TWO digits after comma!
  prms_lists = [
    collect(20.0 : 1.0 : 21.0),  
    collect(20.0 : 1.0 : 21.0),  
    collect(20.0 : 1.0 : 21.0), 
    collect(-0.0 : 1.0 : 0.0), 
    collect(-0.0 : 1.0 : 0.0), 
    collect(-0.0 : 1.0 : 0.0),
    # hint: TC = (700, 750, 800, 850)  => DD = ( 1.277, 2.92, 5.35, 9.05)e-13
    [9.05e-13]
  ]
  scripted_tuple = (0, 1, 0, 0, 0, 0, 0)
  
  TC = 850
  pO2 = [60, 80]
  bias = 0.0

  simulations = ["EIS", "CV"]  
  
  save_dir = "../data/par_studies/"
  #######################################################  
  # preparing bash output ###############################
  
  # if true, no script is called! Just direclty run_par_study_script_wrap()
  direct_bool = true
  
  #bash_command = "sbatch"
  #bash_command = "echo"
  bash_command = "julia"
  
  #mode = "test_one_prms"
  #mode = "only_print"
  mode = "go"
  
  run_file_name = "../snehurka/run_ysz_fitting_par_study-prms-.jl"
  #######################################################
  ####################################################### 
  #######################################################
  #######################################################
  
  
  ### be careful, this has no impact on real included model ... yet  ... :)
  physical_model_name = "GAS_LoMA"
  ###
  
  if mode == "test_one_prms"
    prms_lists = [list[Int64(ceil(end/2.0))] for list in prms_lists]
  end
  
  # counter of output files
  nodes_count = 1
  per_node_count = 1
  for (i, byte) in enumerate(scripted_tuple)
    if byte == 1
      nodes_count *= size(prms_lists[i],1)
    else
      per_node_count *= size(prms_lists[i],1)
    end 
  end
  
  
  # saving metafile

  metafile_string = "-----------  METAFILE for par_study  -------------\n"
  metafile_string = string(metafile_string,"name=", name,"\n")
  metafile_string = string(metafile_string,"physical_model_name=", physical_model_name,"\n")
  prms_strings = [string(item) for item in prms_lists]
  for (i,str) in enumerate(prms_strings)
    metafile_string = string(metafile_string, prms_names[i],"_list=",str,"\n")
  end
  
  metafile_string = string(metafile_string,"#### nodes / pernode_count  =  ", nodes_count," / ",per_node_count," ( = ",per_node_count*3.0/60.0,"m)  #### \n")
  metafile_string = string(metafile_string,"prms_names=", prms_names,"\n")
  metafile_string = string(metafile_string,"scripted_tuple=", scripted_tuple,"\n")
  metafile_string = string(metafile_string,"TC=", TC,"\n")
  metafile_string = string(metafile_string,"pO2=", pO2,"\n")
  metafile_string = string(metafile_string,"bias=", bias,"\n")
  metafile_string = string(metafile_string,"simulations=", simulations,"\n")
  #metafile_string = string(metafile_string,"string(SIM_list)=", string(SIM_list),"\n") 
  metafile_string = string(metafile_string,"save_dir=", save_dir,"\n")
  metafile_string = string(metafile_string,"mode=", mode,"\n")
  metafile_string = string(metafile_string,"bash_command=", bash_command,"\n")
  metafile_string = string(metafile_string,"run_file_name=", run_file_name,"\n")
  metafile_string = string(metafile_string,"prms_lists=", prms_lists,"\n")
  #metafile_string = string(metafile_string,"SIM_list=", SIM_list,"\n")
  
  
  
  if mode == "only_print"
    println(metafile_string)
    return
  end
  
  meta_file_path = string(save_dir, name, "/")
  run(`mkdir -p $(meta_file_path)`)
  write("$(meta_file_path)__metafile_par_study.txt", metafile_string)
  
  println()
  println(metafile_string)
  
  #setting jobs
  if !consistency_check()
    println("ERROR: meta_run_par_study(): shape mismatch (!consistency_check())")
    return throw(Exception)
  end 
  recursive_bash_call([], 1)
    
  println("ok :) ")
end
        


end
