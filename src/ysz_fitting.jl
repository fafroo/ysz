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
# [ ] promyslet velikost Floatu pri importu velkeho mnozstvi dat
# [x] funkci pro par_study, aby umela zmenit sve info (redukovat pO2_list a tak) -> funkce par_study_filter()
# [ ] par_study_plot_the_best() by mela vyhodit jeden figure s dvema subploty... asi ... 
# [ ] prms_names a prms_values dat do jedne tridy
# ---[ ] vlastne jakoby nevim, jestli je to dobry napad?
# [ ] automaticke vybirani vhodne scripted_tuple
# [ ] f_range pouzivat jako Array frekvenci, nikoliv jako trojici cisel
# [x] zobrazovani jmena data_setu i v simple_run 
# [x] obecne zavest data_set jako soucast definice experimentu v SIM
# [ ] konzistentne pridat "f_interval" do vsech EIS
# [x] plot_legend universalne pridat vsude
# [ ] add data_set and other features to CAP
# [ ] make I-V simulation
# [ ] zaradit jednoduchy vztah teplotni zavislosti na odporu
# [ ] dodelat v CV simulaci vse potrebne, co je v EIS. I importovani novych souboru.
# [ ] plot_legend nechat jako volitelny parameter typical_plot ... asi
#
# ##### stacionarni krivky
# [ ] zatim to modeluji jako CV s nizkym voltratem
#
# ##### interpolace ######
# [o] vizualizace trendu chyby mezi interpolanty
# ---[x] zobrazovat soucty odchylek
# ------[x] obrazky po normalizaci vypadaji ruznorode
# ---[o] zkusit vykoukat trend z nenormalizovanych dat
#
#
# #### od Affra ####
# [ ] fitness funkce s maximovou metrikou
#
# 
# #### Fitting process
# [x] prozkoumat experimentalni data (a udelat prislusne procedury)
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
# 
# [!] fitting with one dominating SIM.fitness_factor, but also see the others with lower factor
# [ ] get Fminbox() to work ! ... or do it by my own withing the to_optimize() function 
#######################


using Printf
using PyPlot
#using Plots
using DataFrames
using LeastSquaresOptim
#using Optim
using BlackBoxOptim
using LinearAlgebra

import Base.string

########  import of model_file  ######
#TODO !!!!


######################################

#include("../examples/ysz_experiments.jl")


include("../src/general_supporting_stuff.jl")
include("../src/import_experimental_data.jl")

include("../src/simulations/general_simulation.jl")
include("../src/simulations/CV_simulation.jl")
include("../src/simulations/EIS_simulation.jl")
include("../src/simulations/CAP_simulation.jl")

include("../src/par_study.jl")
#include("../src/interpolation_fitting.jl")
include("../src/EEC_module.jl")






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
  res = DataFrame()
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
        SIM_sim = apply_checknodes(SIM, SIM_raw, SIM.checknodes)
        
        ## TODO !!!
        res = SIM_sim
        ###
        
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
    if use_experiment
      SIM_exp = apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
      if (pyplot > 0)
        typical_plot_exp(SIM, SIM_exp)
      end
    end
    recursive_simple_run_call([], Array{String}(undef,(0)), Array{Float64}(undef,(0)), 1)
  end
  
  if test
    return test_result
  else
    return res
  end
end


# useful wrap
function simple_run(;TC, pO2=1.0, bias=0.0, data_set="MONO", simulations=[], pyplot=0, use_experiment=true, prms_values=[], prms_names=[], 
                         test=false)
    simple_run(get_SIM_list_rectangle(TC, pO2, bias, data_set, simulations); pyplot=pyplot, use_experiment=use_experiment, prms_values=prms_values, prms_names=prms_names, 
                        test=false)
end


###########################################################
###########################################################
#### Fitness_function_interpolation #######################
###########################################################
###########################################################




###########################################################
###########################################################
####     D R T     ########################################
###########################################################
###########################################################



function test_DRT(;lambda=0.0, mode="EEC", TC=800, pO2=80, bias=0.0, R_ohm=1, R1=1, C1=0.001, R2=1, C2=0.0001, alpha=1, prms_names=[], prms_values=[], backward_check=true, draw_semicircles=false, plot_option="Nyq DRT Bode RC", f_range=EIS_get_shared_f_range(), data_set="MONO", 
tau_min_fac=10, tau_max_fac=10, tau_range_fac=2,
peak_merge_tol=0.0, show_legend=true, fig_num=EIS_standard_figure_num)
  

  if data_set=="POLY_OCV_test"
    bias=0.0
    data_set_list = ["POLY_$idx" for idx in [1,2,3]]
  else  
    data_set_list = [data_set]
  end
  
  for TC_item in TC, pO2_item in pO2, bias_item in bias, lambda_item in lambda, data_set_item in data_set_list

      
    
    
    # to add .... , tau_min_fac=tau_min_fac, tau_max_fac=tau_max_fac, tau_range_fac=tau_range_fac
    DRT_control = DRT_control_struct(lambda_item, tau_min_fac, tau_max_fac, tau_range_fac, peak_merge_tol)
    
    SIM_list = EIS_simulation(TC_item, pO2_item, bias_item, data_set=data_set_item, DRT_draw_semicircles=draw_semicircles, 
                DRT_control=DRT_control, plot_option=plot_option, fig_num=fig_num, f_range=f_range)
    SIM = SIM_list[1]
    
    #@show SIM
    
    if mode=="EEC"
      EIS_df = EIS_get_RC_CPE_elements(R1, C1, R2, C2, alpha, R_ohm, f_range=f_range)
      typical_plot_sim(SIM, EIS_df, "! EEC ($R1, $C1) ($R2, $C2, $alpha)")
    elseif mode=="sim"
      EIS_df = ysz_fitting.simple_run(SIM_list, pyplot=1, 
        prms_names=prms_names, 
        prms_values=prms_values, use_experiment=false)
    elseif mode=="exp"
      #@show SIM.checknodes
      EIS_df = apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
      if data_set_item[end-1:end]==".z"
        typical_plot_exp(SIM, EIS_df, "!$(data_set_item)") 
      else
        typical_plot_exp(SIM, EIS_df) 
      end
    end

    
    if backward_check
      DRT_actual = get_DRT(EIS_df, DRT_control)
      println("Fitness error = ",fitnessFunction(EIS_simulation(), DRT_actual.EIS_df, EIS_df))
      typical_plot_sim(EIS_simulation(800, 80, 0.0, use_DRT=false, DRT_control=DRT_control, plot_option=plot_option)..., DRT_actual.EIS_df, "! DRT_backward_check")
      
#       DRT_actual = get_DRT(EIS_df, lambda_item, debug_mode=true)
#       println("Fitness error = ",fitnessFunction(EIS_simulation(), DRT_actual.EIS_df, EIS_df))
#       
#       plot_DRT_h(DRT_actual)
#       typical_plot_sim(EIS_simulation(800, 80, 0.0, use_DRT=false)..., DRT_actual.EIS_df, "! DRT_backward_check")
    end
  end
    
  #return DRT_actual
  return
end





###########################################################
###########################################################
#### Working space ########################################
###########################################################
###########################################################

# function EIS_get_and_plot_RC_element(R, C, Rohm=0)
#   EIS_RC = DataFrame( f = [], Z = [])
#   for f in get_shared_checknodes(EIS_simulation(800,100,0.0)...)
#     push!(EIS_RC, (f, Rohm + R/(1 + im*2*pi*f*R*C)))
#   end
#   typical_plot_sim(EIS_simulation(800,100,0.0)..., EIS_RC, " RC_elem")
#   return EIS_RC
# end


###########################################################
###########################################################
#### Levenberg-Marquardt ##################################
###########################################################
###########################################################

function get_fitting_initial_condition()
  # fitted things for 850, 100, 0, MONO
  #x0 = (0.0, 0.0, 0.0, 20.01931310942335, 19.13133370088284, 20.30899006042895, 1, 21, -0.10582329539200844, -0.5612452246493124, -0.004425970436021694, 9.0e-13, 0.85, true, 0.25, 0.5)
  x0 = (0.0, 0.0, 0.0, 20.009046378779413, 19.133842639499733, 20.243362989853896, 1, 21, -0.13315570684220343, -0.5493240083421153, -0.01362036864906086, 9.0e-13, 0.85, true, 0.25, 0.5)
  #x0 =  (0.0, 0.0, 0.0, 20.304854396292857, 19.100695998282248, 20.203452081086848, 29.786448291451634, 22.50852142987179, -0.09991041679025696, -0.47107633036111757, -0.01598123040079823, 9.0e-13, 0.85, true, 0.25, 0.5)
  #x0 =  (0.0, 0.0, 0.0, 20.304854396292857, 19.100695998282248, 20.303452081086848, 25.786448291451634, 22.50852142987179, -0.09991041679025696, -0.47107633036111757, -0.01598123040079823, 9.0e-13, 0.85, true, 0.25, 0.5)
  
  x0 = (0.0, 0.0, 0.0, 19.85559447494341, 20.057302906172676, 19.30709527857658, 5.129458479008138, 25.342869613876285, -0.45254880302146255, -0.6377067007673729, -0.940802850242693, 9.0e-13, 0.85, true, 0.25, 0.5) # Optim ||  EIS: err =0.00442703203277445 (850, 100, 0.0)
  
  #x0 = (0.0, 0.0, 0.0, 19.907580949649912, 20.028132917464422, 19.403302626940977, 6.0097725379445714, 25.284531020122714, -0.46362112792111193, -0.6359175506862039, -0.9403837928159715, 9.0e-13, 0.85, true, 0.25, 0.5) # Optim ||  EIS: err =0.00442703203277445 (850, 100, 0.0)
  
  #x0 = (0.0, 0.0, 0.0, 19.7936, 20.2862, 18.7722, 5.0, 25.0, -0.508173, -0.648986, -0.795695, 9.0e-13, 0.85, true, 0.25, 0.5)
  
  x0 = (1.0, 0.0, 1.0, 20.85559447494341, 20.057302906172676, 20.30709527857658, 20.129458479008138, 20.342869613876285, -0.25254880302146255, -0.1377067007673729, -0.140802850242693, 9.0e-13, 0.85, true, 0.25, 0.5) # Optim ||  EIS: err =0.00442703203277445 (850, 100, 0.0)
  
  x0 = (1.0, 0.0, 1.0, 24.68638793247588, 16.929296123207436, 22.170921800982974, 18.553478252136934, 19.897488751668668, -0.5769408452228384, 0.12810885586649326, -0.33908733417798587, 9.0e-13, 0.85, true, 0.25, 0.5) # err = 0.03701 .. cosi pro (850, [60, 80, 100], 0.0) MONO
  
  
  x0 = (1.0, 0.0, 1.0, 23.9595, 23.4001, 21.1791, 5.0, 25.0, -0.78899, -0.522923, 0.591665, 9.0e-13, 0.85, true, 0.25, 0.5) # err =  0.025134653 .. cosi pro (850, [60, 80, 100], 0.0) MONO
  
  x0 = (1.0, 0.0, 1.0, 26.429724575986242, 21.599841316855407, 21.084721296249796, -6.663086428479204, 23.342499176410396, -2.123456036083893, -0.6371129649203712, 0.6658969910747834, 9.0e-13, 0.85, true, 0.25, 0.5) # err =  0.013466 .. cosi pro (850, [60, 80, 100], 0.0) MONO
  
  x0 = (1.0, 0.0, 1.0, 20.93358829626075, 21.402618660603686, 21.136328351766778, 5.47009852248677, 23.37319636187158, -0.4140996223040165, -0.6454107557668001, 0.6964208008004982, 9.0e-13, 0.85, true, 0.25, 0.5) # intermediate ... to delete !
  
  x0 = (1.0, 0.0, 1.0, 20.93358829626075, 20.402618660603686, 20.136328351766778, 20.47009852248677, 20.37319636187158, -0.4140996223040165, -0.6454107557668001, 0.6964208008004982, 9.0e-13, 0.85, true, 0.25, 0.5) # intermediate ... to delete !
  
  return x0
end

mutable struct SIM_fitting_struct
  TC
  pO2
  bias
  data_set
  simulations
  #
  SIM_list
  #
  prms_names::Array{String}
  x0
  mask
  lower_bounds
  upper_bounds
  fitted_prms
  #
  BBO_bool::Bool
  iteration_count::Int64
  #
  print_to_file::Bool
  save_dir
  file_name
  
  SIM_fitting_struct() = new()
end

function run_SIM_fitting_script_wrap(;
                    TC_string="850", 
                    pO2_string="100", 
                    bias_string="0.0", 
                    data_set="MONO",
                    simulations_stirng="[\"EIS\"]",
                    #
                    x0_string = string(get_fitting_initial_condition()),
                    prms_names_string=string(["kappaA", "kappaR", "kappaO", 
                                "rA", "rR", "rO",         "rB", "rC",     
                                "DGA", "DGR", "DGO",     
                                "DD", "nu", "separate_vacancy",       "sites_Om0", "sites_Om1"  ]),
                    mask_string        =string((0, 0, 0,
                                    1, 1, 1,        0, 0,
                                    1, 1, 1,
                                    0, 0,     0,        0, 0    )),
                    lower_bounds_string=string((0.0, 0.0, 0.0,
                                  15.5, 15.9, 15.7,        5, 5,       
                                  -0.8, -0.8, -0.8,     
                                  [90]*1.0e-14, collect(0.85 : 0.05 : 0.85), true,       1/4, 1/2    )),
                    upper_bounds_string=string((0.0, 0.0, 0.0,
                                  25.5, 25.9, 25.7,        25, 25,        
                                  0.8, 0.8, 0.8,     
                                  [90]*1.0e-14, collect(0.85 : 0.05 : 0.85), true,       1/4, 1/2    )),
                    #
                    BBO_bool_string=string(false), 
                    iteration_count_string="100",   
                    #
                    print_to_file=string(true), 
                    save_dir="../data/EEC/temp/log/", 
                    file_name="default.txt",
                    #
                    pyplot_string=string(false),  
                    plot_each_x_th_string=string(50),
                    print_only_result_string=string(false),
                    )
  
  SIM_fitting = SIM_fitting_struct()
  #
  TC = eval(Meta.parse(TC_string))
  pO2 = eval(Meta.parse(pO2_string))
  bias = eval(Meta.parse(bias_string))
  simulations = eval(Meta.parse(simulations_string))
  #
  SIM_fitting.TC = TC
  SIM_fitting.pO2 = pO2
  SIM_fitting.bias = bias
  SIM_fitting.simulations = simulations
  #
  SIM_fitting.SIM_list = get_SIM_list_rectangle(TC, pO2, bias, data_set, simulations)
  #
  SIM_fitting.x0 = eval(Meta.parse(x0_string))
  SIM_fitting.prms_names = eval(Meta.parse(prms_names_string))
  SIM_fitting.mask = eval(Meta.parse(mask_string))
  SIM_fitting.lower_bounds = eval(Meta.parse(lower_bounds_string))
  SIM_fitting.upper_bounds = eval(Meta.parse(upper_bounds_string))
  #
  SIM_fitting.BBO_bool = eval(Meta.parse(BBO_bool_string))
  SIM_fitting.iteration_count = eval(Meta.parse(iteration_count_string))
  SIM_fitting.print_to_file = eval(Meta.parse(print_to_file_string))
  
  run_SIM_fitting(SIM_fitting, 
                  pyplot=eval(Meta.parse(pyplot_string)),
                  plot_each_x_th=eval(Meta.parse(plot_each_x_th_string)),
                  print_only_result=eval(Meta.parse(print_only_result_string))
                  )
  return
end

function build_SIM_fitting(;TC=850, pO2=100, bias=0.0, data_set="MONO", simulations=["EIS"],  
                      #
                      x0 = get_fitting_initial_condition(),
                      prms_names=["kappaA", "kappaR", "kappaO", 
                                  "rA", "rR", "rO",         "rB", "rC",     
                                  "DGA", "DGR", "DGO",     
                                  "DD", "nu", "separate_vacancy",       "sites_Om0", "sites_Om1"  ],
                      mask          =(0, 0, 0,
                                      1, 1, 1,        0, 0,
                                      1, 1, 1,
                                      0, 0,     0,        0, 0    ),
                      lower_bounds=(0.0, 0.0, 0.0,
                                    15.5, 15.9, 15.7,        5, 5,       
                                    -0.8, -0.8, -0.8,     
                                    [90]*1.0e-14, collect(0.85 : 0.05 : 0.85), true,       1/4, 1/2    ),
                      upper_bounds=(0.0, 0.0, 0.0,
                                    25.5, 25.9, 25.7,        25, 25,        
                                    0.8, 0.8, 0.8,     
                                    [90]*1.0e-14, collect(0.85 : 0.05 : 0.85), true,       1/4, 1/2    ),
                      #
                      BBO_bool=false, iteration_count=100,   
                      #
                      print_to_file=true, save_dir="../data/EEC/temp/log/", file_name="default.txt",
                      )
  SIM_fitting = SIM_fitting_struct()
  
  SIM_fitting.TC = TC
  SIM_fitting.pO2 = pO2
  SIM_fitting.bias = bias
  SIM_fitting.data_set = data_set
  SIM_fitting.simulations = simulations

  #
  SIM_fitting.prms_names = prms_names
  @show typeof(x0)
  SIM_fitting.x0 = x0
  SIM_fitting.mask = mask
  SIM_fitting.lower_bounds = lower_bounds
  SIM_fitting.upper_bounds = upper_bounds
  #
  SIM_fitting.BBO_bool = false
  SIM_fitting.iteration_count = iteration_count
  #
  SIM_fitting.print_to_file = print_to_file
  SIM_fitting.save_dir = save_dir
  SIM_fitting.file_name = file_name

  return SIM_fitting
end

function run_par_study_script_wrap(
                    prms_lists_string=("[20, 20, 20, 0.0, 0.0, 0.0]"),
                    save_dir="./kadinec/",
                    name="default_par_study_name",
                    scripted_tuple_string=("(1, 0, 0, 0, 0, 0)"),
                    prms_names_string="[\"A0\", \"R0\", \"K0\", \"DGA\", \"DGR\", \"DGO\"]", 
                    TC_string="800",
                    pO2_string="100",
                    bias_string="0.0",
                    simulations_string="[\"EIS\"]",
                    mode="test",
                    physical_model_name="nothing")
  
  
  actual_par_study = par_study_struct()
  
  actual_par_study.physical_model_name = physical_model_name
  actual_par_study.name = name
  actual_par_study.save_dir = save_dir
  actual_par_study.prms_lists = eval(Meta.parse(prms_lists_string))
  actual_par_study.scripted_tuple = eval(Meta.parse(scripted_tuple_string))
  actual_par_study.prms_names = eval(Meta.parse(prms_names_string))
  #
  TC = eval(Meta.parse(TC_string))
  pO2 = eval(Meta.parse(pO2_string))
  bias = eval(Meta.parse(bias_string))
  simulations = eval(Meta.parse(simulations_string))
  #
  actual_par_study.SIM_list = get_SIM_list_rectangle(TC, pO2, bias, simulations)
  
  run_par_study(    actual_par_study,
                    save_dir=save_dir, 
                    save_file_bool=true,
                    mode=mode)
end

function run_SIM_fitting_example()
  SIM_fitting = build_SIM_fitting()
  
  fit = run_SIM_fitting(SIM_fitting, pyplot=false, print_only_result=false)
  return fit
end


function run_SIM_fitting(SIM_fitting::SIM_fitting_struct;
                      pyplot=false,  plot_each_x_th=50,
                      print_only_result=false,
                      )
  function plot_error_projection(prms_values, prms_names, error)
    fig_num = 333
    prms_length = length(prms_names)
    plot_edge = ceil(sqrt(prms_length))
    
    if projection_plot_maximum < 0 && error < 20
      projection_plot_maximum = error
    end
    
    if error < projection_plot_minimum
      projection_plot_minimum = error
    end
    
    figure(fig_num)
    for i in 1:prms_length
      subplot(plot_edge, plot_edge, i)
      xlabel(prms_names[i])
      ylabel("err")
      ylim(projection_plot_minimum, projection_plot_maximum)
      plot(prms_values[i], error, "o")
    end
    
    return
  end
  
  function prepare_prms(mask, x0, x)
      prms = []
      xi = 1
      for i in collect(1 : 1 :length(mask))
          if convert(Bool,mask[i])
              append!(prms, x[xi])
              xi += 1
          else
              append!(prms, x0[i])
          end
      end
      return Tuple(prms)
  end
  
  function check_x_in(x, low, upp)
    for (i, item) in enumerate(x)
      if item < low[i] || upp[i] < item
        return false
      end
    end
    return true
  end
  
  function to_optimize(x)
      if !(check_x_in(x, lowM, uppM))
        print_only_result || (output_buffer *= "    OUT OF THE BOUNDS   \n")
        print_only_result || (print_and_delete_buffer(output_buffer))
        return 10000
      end
      iteration_counter += 1
      
      prms_values = prepare_prms(SIM_fitting.mask, SIM_fitting.x0, x)
#       print(" >> SIM_fitting.mask = ",SIM_fitting.mask)
#       print(" || prms = ",prms_values)

      print_only_result || (output_buffer *= "> prms=$(prms_values)\n")
      err = 0.0
      SIM_err = 0.0
      err_string = ""

      
      for (i, SIM) in enumerate(SIM_fitting.SIM_list)        
        try
          if check_equal_size(SIM_fitting.prms_names, prms_values)  
            SIM_sim = apply_checknodes(SIM, 
                      typical_run_simulation(SIM, SIM_fitting.prms_names, prms_values, pyplot ? 1 : 0),
                      SIM.checknodes)
          end    
                
          if pyplot && mod(iteration_counter, plot_each_x_th) == 1
            typical_plot_sim(SIM, SIM_sim)
          end
          
          SIM_err = SIM.fitness_factor * fitnessFunction(SIM, SIM_exp[i], SIM_sim)
          err += SIM_err
        catch e
          if e isa InterruptException
            rethrow(e)
          else
            print(e)
            err += 1000*norm(x0M.-x) 
          end
        end 
        
        err_string = " "*string(SIM)
        print_only_result || (output_buffer *= "$(err_string)=$(SIM_err) || ")
      end
      print_only_result || (output_buffer *= "err=$(err)\n\n")
      
      
      if pyplot
        plot_error_projection(
          take_only_masked(SIM_fitting.mask, prms_values), 
          take_only_masked(SIM_fitting.mask, SIM_fitting.prms_names), 
          err)
      end
      
      print_only_result || (print_and_delete_buffer(output_buffer))
      
#       println(" || ",err_string,": err =", err)

      return err
  end

  function take_only_masked(mask, array)
    output = []
    for (i, item) in enumerate(mask)
      if item == 1
        append!(output, [array[i]])
      end
    end
    return output
  end
  
  function print_and_delete_buffer(buffer="string")
    if SIM_fitting.print_to_file
      open(SIM_fitting.save_dir*SIM_fitting.file_name, "a") do f
        write(f, buffer)
      end
    else
      print(buffer)
    end
  end

  
  
  physical_model_name="necum"
  iteration_counter = 0
  
  if pyplot
      PyPlot.close()
  end
  
  SIM_exp = Array{DataFrame}(undef, length(SIM_fitting.SIM_list))
  for (i, SIM) in enumerate(SIM_fitting.SIM_list)
    SIM.plot_legend=false
    SIM_exp[i] = apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
    if (pyplot > 0)
      typical_plot_exp(SIM, SIM_exp[i])
    end
  end

  output_buffer = ""  
  if SIM_fitting.print_to_file
    run(`mkdir -p $(SIM_fitting.save_dir)`)
  end
  #######################################
  #######################################
  #######################################
  #######################################
  # control panel
  #SIM_fitting.x0 = zeros(2)
  #optimize(rosenbrock, SIM_fitting.x0, LevenbergMarquardt())
  

  ######################################            

  if length(SIM_fitting.x0) != length(SIM_fitting.prms_names)
    output_buffer *= ("ERROR: length(SIM_fitting.x0) != length(SIM_fitting.prms_names)")
    print_and_delete_buffer(output_buffer)
    return throw(Exception)
  end

  
  #######################################
  #######################################
  #######################################

  #######################################
  #######################################
  #######################################
  #######################################
  
  projection_plot_maximum = -1
  projection_plot_minimum = Inf
  
  x0M = zeros(0)
  lowM = zeros(0)
  uppM = zeros(0)
  for i in collect(1 : 1 : length(SIM_fitting.mask))
      if convert(Bool,SIM_fitting.mask[i])
          append!(x0M, SIM_fitting.x0[i])
          append!(lowM, SIM_fitting.lower_bounds[i])
          append!(uppM, SIM_fitting.upper_bounds[i])
      end
  end
  
  function get_SearchRange(lowM, uppM)
    output = Array{Tuple{Float64, Float64}}(undef, 0)
    for i in 1:length(lowM)
      append!(output, [(lowM[i], uppM[i])])
    end
    return output
  end
  
  function rosenbrock2d(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
  end

  #bboptimize(rosenbrock2d, SearchRange = [(-5.0, 5.0), (-2.0, 2.0)])
  
  # for BlackBoxOptim
  SearchRange = get_SearchRange(lowM, uppM)
  #@show SearchRange
  
  if SIM_fitting.BBO_bool
    fit= bboptimize(to_optimize, SearchRange = SearchRange)
  else
#     fit = optimize(
#       to_optimize,
#       x0M,
#       #rosenbrock2d,
#       #[0., 0.],
#       #NelderMead(),
#       ParticleSwarm(lower=lowM, upper=uppM, n_particles = 10) 
#       )
      
    fit = optimize(to_optimize, x0M, NelderMead(),
                     Optim.Options(
                      iterations = SIM_fitting.iteration_count,
                      
                      )
                     )
     
#    fit = optimize(to_optimize, lowM, uppM, x0M, Fminbox(NelderMead()), 
#                     Optim.Options(
#                      iterations = 10,
#                      
#                       )
#                     )
  end
  
  SIM_fitting.prms_values = prepare_prms(SIM_fitting.mask, SIM_fitting.x0, fit.minimizer)
  
  if print_only_result
    output_buffer *= "\n"
    output_buffer *= "SIM_list = $(string.(SIM_fitting.SIM_list))\n"
    output_buffer *= "data_set = $(string.(SIM_fitting.data_set))\n"
    output_buffer *= "prms_names = $(SIM_fitting.prms_names)\n"
    output_buffer *= "x0 = $(SIM_fitting.x0)\n"
    output_buffer *= "lower_bounds = $(SIM_fitting.lower_bounds)\n"
    output_buffer *= "upper_bounds = $(SIM_fitting.upper_bounds)\n"
    output_buffer *= "prms_values = $(SIM_fitting.prms_values)\n"
    output_buffer *= "mask = $(SIM_fitting.mask)\n"
    output_buffer *= "err=$(fit.minimum)\n"
    print_and_delete_buffer(output_buffer)
  end
  
  #println(optimize(to_optimize, x0M, lower=lowM, upper=uppM, Δ=1000000, f_tol=1.0e-14, g_tol=1.0e-14, LevenbergMarquardt()))
  
  ####optimize(to_optimize, x0M, Dogleg())
  return fit  
end





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
  
  express3_bool = true

  if express3_bool
    run_file_name = "../snehurka/run_EX3_ysz_fitting_par_study-prms-.jl"
  else
    run_file_name = "../snehurka/run_ysz_fitting_par_study-prms-.jl"
  end
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
