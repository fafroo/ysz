module ysz_fitting
#######################
####### TODO ##########
# [x] better ramp ... starting directly from steadystate :)
# [x] general global search using projection to each variable
# [x] compute EIS exacly on checknodes and therefore remove plenty of "EIS_apply_checknodes"
# [ ] put appropriate and finished stuff into "CV_fitting_supporting_stuff"
# [ ] spoustet run_new() v ruznych procedurach stejnou funkci "EIS_default_run_new( ... )"
# [x] implement LM algorithm
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
# ------[!] mozna by se include mel vykonavat na vyssi urovni a do ysz_experiments preposilat uz hotovou tridu
# [x] vymyslet lepe zadavani parametru pro ruzne modely s ruznymi parametry 
# [o] opravdu promyslet to objektove programovatni ... CV_sim, EIS_sim ... prms_lists ...
# ---[o] pridat treba dalsi experiment s dalsimi daty, vuci kterym se da srovnavat (kapacitance)
# ---[ ] udelat tridu par_study
# [ ] PAR_STUDY vyhodnoceni ... soubory do slozky dane studie .. automatizovane
# ---[ ] pri zobrazovani dat udelat clustery dle chyby/prms a pak zobrazit (prms ci Nyquist) 1 reprezentanta z kazde
# ---[ ] automatizovat ukladani souboru vysledku par_study.vyhodnoceni()
# [ ] do experiments.jl pridat obecne zaznamenavani promennych od final_plot
# [x] get rig od shared_prms and shared_add_prms !!!
# [ ] snehurka by rada ./zabal.sh a pocitac zase ./rozbal_par_study.sh
# [!] snehurkove fitovani by slo zrychlit, kdyz bych skriptoval EQ parametry a job by menil jen kineticke? ... chrm ...
# ---[!!!!!!] nejak si preposilat steadystate?! 
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
using DataFrames
using LeastSquaresOptim

import Base.string

########  import of model_file  ######
#TODO !!!!


######################################

#include("../examples/ysz_experiments.jl")

include("../src/import_experimental_data.jl")
include("../src/incommon_simulation.jl")
include("../src/CV_simulation.jl")
include("../src/EIS_simulation.jl")
include("../src/CV_fitting_supporting_stuff.jl")
include("../src/EIS_fitting_supporting_stuff.jl")






function get_fitted_all_prms()
  prms_names=["A0", "R0", "K0", "SA", "SR", "SO", "DGA", "DGR", "DGO", "betaA", "betaR", "betaO", "DD"]
  prms_values=[19.7, 19.7, 18.6,    1, 1, 1,    0.7, -0.8, -0.3,      0.5, 0.5, 0.5,    5.35e-13]  # fitted to EIS 800, 100, 0.0
    
  return prms_names, prms_values
end


####
####


function for_each_prms_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else
      for parameter in prms_lists[active_idx]
        recursive_call(push!(deepcopy(output_set),parameter), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end

function for_each_indicies_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else
      for i in 1:size(prms_lists[active_idx],1)
        recursive_call(push!(deepcopy(output_set),i), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end

function get_prms_from_indicies(prms_lists, tuple)
  output_array = []
  for (i,list) in enumerate(prms_lists)
    append!(output_array, list[tuple[i]])
  end
  return output_array
end


# #############################################################
# ############################################################
# ###################   TO ANOTHER FILE ->>> par_study.jl ####

mutable struct par_study_struct
  physical_model_name::String
  simulation_list::Array{abstract_simulation}
  prms_names::Array{Any}
  prms_lists::Array{Any}
  scripted_tuple::Tuple
  simulation_data::Array{Any}
  simulation_err::DataFrame
  #
  name::String
  #
  par_study_struct() = new()
end

# function par_study()
#   output = []
#   for TC_item in TC
#     for pO2_item in pO2
#       for bias_item in bias
#         this = EIS_simulation()
#         
#         this.TC = TC_item
#         this.pO2 = pO2_item
#         this.bias = bias_item
#         this.dx_exp = dx_exp
#         this.omega_range = omega_range
#         this.fig_size = fig_size
#         
#         push!(output, this)
#       end
#     end
#   end
#   return output
# end


#############################################################
#############################################################
#############################################################







############## import large amount of data ##########
function import_par_study_from_metafile(;save_dir="../snehurka/data/", name="", metafile_name="__metafile_par_study.txt", CV_bool=false, EIS_bool=false, verbose=false)
  prms_lists = []
  prms_names = []
  scripted_tuple = []
  TC=-600
  pO2=-600
  EIS_bias=-600
  
  file_lines = readlines(string(save_dir,name,"/",metafile_name))
  for str in file_lines
    if occursin('=',str)
      (token_name,token_value)=split(str,"=")
      token_name=="prms_lists" && (prms_lists = eval(Meta.parse(token_value)))
      token_name=="prms_names" && (prms_names = eval(Meta.parse(token_value)))
      token_name=="scripted_tuple" && (scripted_tuple = eval(Meta.parse(token_value)))
      token_name=="TC" && (TC = eval(Meta.parse(token_value)))
      token_name=="pO2" && (pO2 = eval(Meta.parse(token_value)))
      token_name=="EIS_bias" && (EIS_bias = eval(Meta.parse(token_value)))
    end
  end

  (
    import_par_study(save_dir=save_dir, name=name, prms_lists=prms_lists, prms_names=prms_names, scripted_tuple=scripted_tuple, CV_bool=CV_bool,  EIS_bool=EIS_bool, verbose=verbose), 
    prms_lists, 
    scripted_tuple, 
    prms_names, 
    TC, 
    pO2, 
    EIS_bias
  )
end


function import_par_study(;save_dir="../snehurka/data/", name="", prms_lists=[13, 13, 0.10, 0.10, 0.4, 0.0], prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), scripted_tuple=(1,1,0,0,0,0), EIS_bool=false, CV_bool=false, verbose=false)
  
  
  save_dir_forwarded = string(save_dir, name, "/")
  if CV_bool
    CV_hypercube = Array{DataFrame}(undef, Tuple([size(list,1) for list in prms_lists]))
    CV_good_count = 0
    CV_bad_count = 0
  end
  
  if EIS_bool
    EIS_hypercube = Array{DataFrame}(undef, Tuple([size(list,1) for list in prms_lists]))
    EIS_good_count = 0
    EIS_bad_count = 0
  end

  counter=1
  for list in prms_lists
    counter*= size(list,1)
  end
  prmsDIV100 = counter/100
  
  counter=0
  println("Loading prms_lists")
  println("|====================================================================================================|")

  function perform_import(prms_indicies)
    
    counter += 1
    if counter >= prmsDIV100 - 0.00000001
      print("o")
      counter=1
    end
    
    if CV_bool
      try
        CV_hypercube[prms_indicies...] = CV_load_file_prms(save_dir=save_dir_forwarded, prms=get_prms_from_indicies(prms_lists, prms_indicies), prms_names=prms_names, scripted_tuple=scripted_tuple, throw_exception=true, verbose=verbose)
        CV_good_count += 1
      catch e
        if e isa InterruptException
          rethrow(e)
        else
          CV_hypercube[prms_indicies...] = DataFrame()
          CV_bad_count += 1
        end
      end
    end
    if EIS_bool
      try
        EIS_hypercube[prms_indicies...] = EIS_load_file_prms(save_dir=save_dir_forwarded, prms=get_prms_from_indicies(prms_lists, prms_indicies), prms_names=prms_names, scripted_tuple=scripted_tuple, throw_exception=true,verbose=verbose)
        EIS_good_count += 1
      catch e
        if e isa InterruptException
          rethrow(e)
        else
          EIS_hypercube[prms_indicies...] = DataFrame()
          EIS_bad_count += 1
        end
      end
    end
  end
  
  print("|")
  for_each_indicies_in_prms_lists(prms_lists, perform_import)
  println("|")
  
  if CV_bool && EIS_bool
    println("import_par_study: CV_good / all = ",CV_good_count, " / ", CV_good_count + CV_bad_count)
    println("import_par_study: EIS_good / all = ",EIS_good_count, " / ", EIS_good_count + EIS_bad_count)
    return (CV_hypercube, EIS_hypercube)
  end
  if CV_bool
    println("import_par_study: CV_good / all = ",CV_good_count, " / ", CV_good_count + CV_bad_count)
    return CV_hypercube
  end
  if EIS_bool
    println("import_par_study: EIS_good / all = ",EIS_good_count, " / ", EIS_good_count + EIS_bad_count)
    return EIS_hypercube
  end
end



function get_error_projection_to_prms(sim::EIS_simulation, EIS_hypercube, prms_lists, TC, pO2, EIS_bias,  prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), omega_range=EIS_get_shared_omega_range())
  #scripted_tuple = (1,1,0,0,0,0)
#prms_lists=([13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0], [13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0], [-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7], [-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7], 0.4, [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
  #name = "prvni_vrh_lin_ads"
  
  checknodes = EIS_get_checknodes_geometrical(omega_range...)
  
  EIS_exp = EIS_apply_checknodes(import_EIStoDataFrame(TC=TC, pO2=pO2, bias=EIS_bias), checknodes)
  
  error_df = DataFrame(error = [], prms_indicies = [])
  
  counter=1
  for list in prms_lists
    counter*= size(list,1)
  end
  prmsDIV100 = counter/100
  
  counter=0
  println("performing error projection")
  println("|====================================================================================================|")
  
  function perform_error_projection(prms_indicies)
    
    counter += 1
    if counter >= prmsDIV100 - 0.00000001
      print("o")
      counter=1
    end
    
    if size(EIS_hypercube[prms_indicies...],1) > 0
      #@show EIS_hypercube[prms_indicies...]
      #EIS_sim = EIS_apply_checknodes(EIS_hypercube[prms_indicies...], EIS_get_shared_checknodes())
      #close()
      #Nyquist_plot(EIS_exp)
      #Nyquist_plot(EIS_sim)      
      #println("OK")
      push!(error_df, 
        (EIS_fitnessFunction(
          EIS_exp,
          EIS_apply_checknodes(EIS_hypercube[prms_indicies...], checknodes)
          )
          , 
          prms_indicies
        )
      )
    end
  end

  print("|")
  for_each_indicies_in_prms_lists(prms_lists, perform_error_projection)
  println("|")
  
  sort!(error_df, :error)
  

  error_df
end



###################################################


function display_err_projection(error_df, prms_lists, prms_names, count)
  # TODO
  # [ ] mention the par_study name in the title (and checknodes)
  
  range_to_display = 1:count
  figure(4)
  suptitle("The best $count")
  for i in 1:6
    subplot(230+i)
    title(prms_names[i])
    plot([prms_lists[i][item[i]] for item in error_df.prms_indicies[range_to_display]], error_df.error[range_to_display], "x")
  end
end

function display_the_best(err_df, EIS_hypercube, prms_lists, TC, pO2, EIS_bias; from=1, till=10)
  # TODO
  # [ ] mention the par_study name in the title ( and checknodes)
  
  figure(3)
  EIS_exp = EIS_apply_checknodes(import_EIStoDataFrame(TC=TC, pO2=pO2, bias=EIS_bias),EIS_get_shared_checknodes())
  Nyquist_plot(EIS_exp, "exp $(EIS_experiment_legend(TC, pO2, EIS_bias))")
  for i in from:till
    Nyquist_plot(EIS_hypercube[err_df.prms_indicies[i]...], "$(i): $(get_prms_from_indicies(prms_lists, err_df.prms_indicies[i]))")
  end
end









#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
# function new_common_DD()
#   4.35e-13
# end

function check_equal_size(list_1, list_2)
  if size(list_1,1) == size(list_2,1)
    return true
  else
    println("ERROR: check_equal_size: shape mismatch $(list_1) and $(list_2)")
    return throw(Exception)
  end
end

function simple_run(SIM_list::Array{abstract_simulation}; pyplot=0, use_experiment=true, prms_values=[], prms_names=[], 
                        make_EIS_hypercube=false, test=false)
  if test
    test_result = 0
  end
  
  ## TODO !!!!!!!
  
#   if make_EIS_hypercube
#     EIS_hypercube = Array{DataFrame}(undef, Tuple([size(list,1) for list in prms_lists]))
#   end
  
  for SIM in SIM_list
    function recursive_simple_run_call(output_prms, plot_names, plot_values, active_idx)
      if active_idx > size(prms_names,1)
    
        # perform the standard code #####################
        prms_names_in=[]
        prms_values_in=[]
        
        append!(prms_names_in, prms_names)
        append!(prms_values_in, output_prms)
      
      
        if check_equal_size(prms_names, prms_values)  
          SIM_df = typical_run_simulation(SIM, prms_names_in, prms_values_in, pyplot)        
        end    
        #@show EIS_df
        if size(plot_names,1) < 1
          plot_prms_string = ""
        else
          plot_prms_string = " $(string(plot_names)) = $(plot_values)"
        end
        SIM_sim = apply_checknodes(SIM, SIM_df, checknodes)
#             if make_EIS_hypercube
#               EIS_hypercube[output_prms...]
#             end
        if pyplot > 0
            typical_plot_sim(SIM, SIM_sim, plot_prms_string)
        end  
        if use_experiment
          if test
            test_result+=fitnessFunction(SIM, SIM_exp, SIM_sim)
          else
            fitting_report(SIM, plot_prms_string, SIM_exp, SIM_sim)
          end
        end
        return
        #####################################  
      
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
    
    # here starts the body
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
  end
end



#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################

# # # function get_dfs_to_interpolate(
# # #                             EIS_bool=true, pyplot=false,
# # #                            TC=800, pO2=100, EIS_bias=0.0,
# # #                             #rprm1=range(15, stop=20, length=3),
# # #                             #rprm2=range(15, stop=20, length=3),
# # #                             #rprm1=range(18.4375, stop=18.59375, length=3),
# # #                             #rprm2=range(20.59375, stop=21, length=3),
# # #                             rprm1=range(19.3, stop=21.0, length=2),
# # #                             #rprm2=range(15.0, stop=23, length=7),
# # #                             rprm2=range(-3.0, stop=-1.6, length=2),
# # #                             nx=100, ny=100,
# # #                             wp=[1, 1, 0, 0, 0, 0], #TODO !!!
# # #                             depth_remaining=1,
# # #                             approx_min=Inf,
# # #                             approx_min_position=[0, 0],
# # #                             recursive=false,
# # #                             recursive_string="",
# # #                             )
# # # 
# # #     println("nx = ",nx, " ..... ny = ",ny)
# # #     #PyPlot.close()
# # #     #PyPlot.close()
# # #     
# # #     ######################################################
# # #     #rprm1 = collect(-5.0 : 1.0 : 5.0)
# # #     #rprm1 = collect(-3. : 0.5 : 3.0)
# # #     #rprm1=range(15, stop=20, length=3)
# # #     
# # #     #rprm2 = collect(-3. : 0.5 : 3.0)
# # #     #rprm2=range(15, stop=20, length=3)
# # #     ######################################################
# # #     
# # #     wpn = ["A0","R0","DGA","DGR","betaR","SR"]
# # #     #A0 = 21.71975544711280
# # #     #R0 = 19.53
# # #     A0 = 19.50
# # #     R0 = 19.85
# # #     DGA = 0.0905748
# # #     DGR = -0.708014
# # #     beta = 0.6074566741435283
# # #     A = -0.356
# # #     
# # #     (A0, R0, DGA, DGR, beta, A) = (21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1)
# # # 
# # #     println(recursive_string,"computing 2D scan... ")
# # #     lrprm1 = size(rprm1,1)
# # #     lrprm2 = size(rprm2,1)
# # #     
# # #     global_min = Inf
# # #     global_min_position = [0,0]
# # #     
# # #     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
# # #     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
# # #     if EIS_bool
# # #         println(recursive_string,"Computing EISs")
# # #     else
# # #         println(recursive_string,"Computing CVs")
# # #     end
# # #     @time(
# # #         for i1 in 1:lrprm1
# # #             for i2 in 1:lrprm2
# # #                 try
# # #                     print(recursive_string,rprm1[i1], " / ", rprm2[i2])
# # #                     
# # #                     ##########################################
# # #                     R0 = rprm1[i1]
# # #                     A = rprm2[i2]
# # #                     ##########################################
# # #                     
# # #                    if EIS_bool 
# # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # #                             pyplot=false, EIS_IS=true, out_df_bool=true, 
# # #                             tref=0,
# # #                             prms_in=[A0, R0, DGA, DGR, beta, A],  
# # #                             nu_in=0.9, pO2_in=1.0, DD_in=9.5658146540360312e-10
# # #                         )
# # #                     else
# # #                         
# # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # #                             out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
# # #                             prms_in=[A0, R0, DGA, DGR, beta, A]
# # #                         )
# # #                     end
# # #                     ok_matrix[i1, i2] = 1
# # #                     println("  .. ok :)")
# # #                 catch e
# # #                     if e isa InterruptException
# # #                         rethrow(e)
# # #                     else
# # #                         ok_matrix[i1, i2] = 0
# # #                         println("  :(   ")
# # #                         continue
# # #                     end
# # #                 end
# # #             end
# # #         end
# # #     )
# # #     
# # #     #PyPlot.figure(figsize=(6,4))
# # #     #xlabel("parameter")
# # #     #zlabel("error")
# # #     #title("Fitness Function")
# # # 
# # #     if pyplot && !(recursive)
# # #         figure(figsize=(3,3))
# # #     end
# # #     #print_ok_matrix(ok_matrix)
# # #     
# # #     println(recursive_string,"Preparing for interpolation ... ")
# # #         for i1 in 1:lrprm1-1
# # #             for i2 in 1:lrprm2-1
# # #                 if all_four_nodes_ok(ok_matrix,i1,i2)
# # #                     print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
# # #                     
# # #                     if EIS_bool
# # #                             EIS_list = [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]]
# # #                             EIS_with_checknodes_list = []
# # #                             checknodes =  EIS_get_shared_checknodes()
# # #                             for i in 1:size(EIS_list,1)
# # #                                 push!(EIS_with_checknodes_list, DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_list[i], checknodes)))
# # #                             end
# # #                             return EIS_with_checknodes_list
# # #                     else
# # #                         # TODO !!!
# # #                         FF=CV_get_FF_interp_2D(
# # #                             [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
# # #                             rprm1[i1:i1+1],
# # #                             rprm2[i2:i2+1];
# # #                             nx=nx, ny=ny,
# # #                             TC=TC, pO2=pO2
# # #                         )
# # #                     end
# # #                 end 
# # #             end
# # #         end
# # # end

function EIS_view_interpolation_at(Q_list, x, y)
    EIS_intrp = DataFrame(
        f = Q_list[1].f,
        Z = EIS_biliComb(
            Q_list[1],
            Q_list[2],
            Q_list[3],
            Q_list[4],
            x,y).Z
    )
    Nyquist_plot(EIS_intrp, string("x/y = ",x,"/",y))
end

function plot_all(df_list, prms_list)
    PyPlot.figure(figsize=(6,4))
    #PyPlot.figure(figsize=(10,8))
    for i in 1:size(df_list,1)
        #println(i)
        CV_plot(
            df_list[i],
            string("prms = ",prms_list[i])
        )
    end
    
    CV_orig = import_CVtoDataFrame(;TC=800,pO2=100)
    #checknodes = CV_get_checknodes(0.01,0.95,-0.95,-0.01,0.05)
    checknodes = CV_get_shared_checknodes()
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_orig, checknodes))
    
    CV_plot(CV_exp, CV_experiment_legend)
end

function plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, TC, pO2)
    PyPlot.figure(figsize=(6,4))
    
    
    for i1 in 1:size(rprm1,1)
        for i2 in 1:size(rprm2,1)
            if Bool(ok_matrix[i1,i2])
                #println(" ploting i1 / i2 : ",i1, " / ", i2)
                CV_plot(
                    CV_matrix[i1,i2],
                    string("prms = ", rprm1[i1], " / ", rprm2[i2])
                )
            end
        end
    end

    CV_orig = import_CVtoDataFrame(;TC=TC,pO2=pO2)
    #checknodes = CV_get_checknodes(0.01,0.95,-0.95,-0.01,0.05)
    checknodes = CV_get_shared_checknodes()
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_orig, checknodes))
    
    CV_plot(CV_exp, "exp $(CV_experiment_legend(TC, pO2))")
end



function plot_EIS_matrix(EIS_matrix, ok_matrix, rprm1, rprm2, TC, pO2, EIS_bias=0.0)
    PyPlot.figure(figsize=(6,4))
    
    checknodes =  EIS_get_shared_checknodes()  
    
    for i1 in 1:size(rprm1,1)
        for i2 in 1:size(rprm2,1)
            if Bool(ok_matrix[i1,i2])
                #println(" ploting i1 / i2 : ",i1, " / ", i2)
                Nyquist_plot(
                    EIS_apply_checknodes(EIS_matrix[i1, i2], checknodes),
                    string("prms = ", rprm1[i1], " / ", rprm2[i2])
                )
            end
        end
    end

    EIS_raw = import_EIStoDataFrame(;TC=TC,pO2=pO2,bias=EIS_bias)
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    Nyquist_plot(EIS_exp, "exp $(EIS_experiment_legend(TC, pO2, EIS_bias))")
end


#####################
## 2D optimization ##
#####################

mutable struct FF_2D
    x_matrix
    y_matrix
    err_matrix
end



function surfplot_FF_2D(FF::FF_2D; my_title="Objective function")
    surf(FF.x_matrix, FF.y_matrix, FF.err_matrix, cstride=1, #cmap=ColorMap("gray"), 
        alpha=0.8);
    xlabel("prm1")
    ylabel("prm2")
    title(my_title)
    zlabel("error")
end


function CV_get_FF_interp_2D(CV_list, prm1_list, prm2_list; nx=20, ny=20, TC=800, pO2=100)
    # ordered for prms A, B as [(A1,B1), (A1,B2), (A2,B1), (A2, B2)]
    # prms_list ordered as [A1 B1 A2 B2]
    CV_raw = import_CVtoDataFrame(;TC=TC,pO2=pO2)
    checknodes = CV_get_checknodes(0.05,0.95,-0.95,-0.05,0.01)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_raw, checknodes))
    
    CV_with_checknodes_list = []
    
    for i in 1:size(CV_list,1)
        push!(CV_with_checknodes_list, DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[i], checknodes)))
    end
    
    # bilinear interpolation
    rx = range(0.0, stop=1.0, length=nx)
    ry = range(0.0, stop=1.0, length=ny)
    
    x_matrix = zeros(nx,ny)
    y_matrix = zeros(nx,ny)
    err_matrix = zeros(nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            x_matrix[ix, iy] = rx[ix]*(prm1_list[2]-prm1_list[1]).+prm1_list[1]
            y_matrix[ix, iy] = ry[iy]*(prm2_list[2]-prm2_list[1]).+prm2_list[1]
            err_matrix[ix, iy] = CV_fitnessFunction(
                    CV_exp,
                    DataFrame(
                        U = checknodes[:,1], 
                        I = biliComb_CVs(
                            CV_with_checknodes_list[1],
                            CV_with_checknodes_list[2],
                            CV_with_checknodes_list[3],
                            CV_with_checknodes_list[4],
                            rx[ix],ry[iy]).I
                    )
                )
        end
    end
    
    return FF_2D(x_matrix, y_matrix, err_matrix)
end



function EIS_get_FF_interp_2D(EIS_list, prm1_list, prm2_list; nx=20, ny=20, TC=800, pO2=100, EIS_bias=0.0)
    # ordered for prms A, B as [(A1,B1), (A1,B2), (A2,B1), (A2, B2)]
    # prms_list ordered as [A1 B1 A2 B2]
    EIS_raw = import_EIStoDataFrame(;TC=TC,pO2=pO2,bias=EIS_bias)
    checknodes =  EIS_get_shared_checknodes()
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    EIS_with_checknodes_list = []
    
    for i in 1:size(EIS_list,1)
        push!(EIS_with_checknodes_list, DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_list[i], checknodes)))
    end
    
    # bilinear interpolation
    rx = range(0.0, stop=1.0, length=nx)
    ry = range(0.0, stop=1.0, length=ny)
    
    x_matrix = zeros(nx,ny)
    y_matrix = zeros(nx,ny)
    err_matrix = zeros(nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            x_matrix[ix, iy] = rx[ix]*(prm1_list[2]-prm1_list[1]).+prm1_list[1]
            y_matrix[ix, iy] = ry[iy]*(prm2_list[2]-prm2_list[1]).+prm2_list[1]
            err_matrix[ix, iy] = EIS_fitnessFunction(
                    EIS_exp,
                    DataFrame(
                        f = checknodes[:], 
                        Z = EIS_biliComb(
                            EIS_with_checknodes_list[1],
                            EIS_with_checknodes_list[2],
                            EIS_with_checknodes_list[3],
                            EIS_with_checknodes_list[4],
                            rx[ix],ry[iy]).Z
                    )
                )
        end
    end
    
    return FF_2D(x_matrix, y_matrix, err_matrix)
end


function all_four_nodes_ok(ok_matrix, x, y)
    for i in 0:1
        for j in 0:1
            if ok_matrix[x+i, y + j] == 0
                return false
            end
        end
    end
    return true
end

function print_ok_matrix(ok_matrix)
println("... printing ok_matrix ...")
    for i2 in 1:size(ok_matrix,2)
        for i1 in 1:size(ok_matrix,1)
        
            print("|", ok_matrix[i1, i2])
        end
        print("|\n")
    end
end

function find_mins(FF, dx; debug_print_bool=false)
    hard_min = Inf
    count_of_mins = 0
    values = []
    positions = []
    curv = []
    
    if size(FF.err_matrix,1) < 3 || size(FF.err_matrix,2) < 3
        println("ERROR: is_min_there: some dimension of FF.err_matrix is smaller than 3")
        return throw(Exception)
    end
    
    for i1 in 2:size(FF.err_matrix,1)-1
        for i2 in 2:size(FF.err_matrix,2)-1
            pot_min = FF.err_matrix[i1, i2]
            if pot_min < hard_min
                hard_min = pot_min
            end
            
            if debug_print_bool
                println(i1, " ", i2, " > ",pot_min)
                
                println(pot_min < FF.err_matrix[i1-1, i2], " ",
                    pot_min < FF.err_matrix[i1+1, i2], " ",
                    pot_min < FF.err_matrix[i1, i2 - 1], " ",
                    pot_min < FF.err_matrix[i1, i2 + 1])
            end
                
            if (pot_min < FF.err_matrix[i1 - 1, i2 - 1] &&
                pot_min < FF.err_matrix[i1 - 1, i2] &&
                pot_min < FF.err_matrix[i1 - 1, i2 + 1] &&
                pot_min < FF.err_matrix[i1, i2 - 1] &&
                pot_min < FF.err_matrix[i1, i2 + 1] &&
                pot_min < FF.err_matrix[i1 + 1, i2 - 1] &&
                pot_min < FF.err_matrix[i1 + 1, i2] &&
                pot_min < FF.err_matrix[i1 + 1, i2 + 1]
                )
                
                c1 = (FF.err_matrix[i1 - 1, i2] - 2*pot_min + FF.err_matrix[i1 + 1, i2])/(dx[1]*dx[1])
                c2 = (FF.err_matrix[i1, i2 - 1] - 2*pot_min + FF.err_matrix[i1, i2 + 1])/(dx[2]*dx[2])
                c12= (FF.err_matrix[i1 - 1, i2 - 1] + FF.err_matrix[i1 + 1, i2 + 1] 
                    - FF.err_matrix[i1 + 1, i2 - 1] - FF.err_matrix[i1 - 1, i2 + 1])/(4*dx[1]*dx[2])
                count_of_mins += 1
                append!(values,pot_min)
                push!(positions, [FF.x_matrix[i1, i2], FF.y_matrix[i1, i2]])
                push!(curv,[c1, c2,c12])
            end
        end
    end
    
    return count_of_mins, values, positions, curv, hard_min
end

function pickup_min(min_count, min_values, min_positions)
    value = Inf;
    position = [0, 0]
    for i in 1:min_count
        if min_values[i] < value
            value = min_values[i]
            position = min_positions[i]
        end
    end
    if value == Inf
        println("ERROR: pickup_min: value = -1")
        throw(Exception)
    else
        return value, position
    end
end

function scan_2D_recursive(;pyplot=false, 
                            EIS_bool=true,
                            TC=800, pO2=100, EIS_bias=0.0,
                            
                            wp=[1, 1, 0, 0, 0, 0], #TODO !!!
                            rprm1=range(0.1, stop=0.9, length=10),
                           #rprm2=range(15.0, stop=23, length=7),
                            rprm2=range(0.8, stop=1.0, length=4),
                            
                            nx=30, ny=30,

                            depth_remaining=1,
                            approx_min=Inf,
                            approx_min_position=[0, 0],
                            recursive=false,
                            recursive_string="",
                            )

    println("nx = ",nx, " ..... ny = ",ny)
    #PyPlot.close()
    #PyPlot.close()
    
    ######################################################
    #rprm1 = collect(-5.0 : 1.0 : 5.0)
    #rprm1 = collect(-3. : 0.5 : 3.0)
    #rprm1=range(15, stop=20, length=3)
    
    #rprm2 = collect(-3. : 0.5 : 3.0)
    #rprm2=range(15, stop=20, length=3)
    ######################################################
    
    wpn = ["A0","R0","DGA","DGR","betaR","SR"]
    
   
    println(recursive_string,"computing 2D scan... ")
    lrprm1 = size(rprm1,1)
    lrprm2 = size(rprm2,1)
    
    global_min = Inf
    global_min_position = [0,0]
    
    CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
    ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
    if EIS_bool
        println(recursive_string,"Computing EISs")
    else
        println(recursive_string,"Computing CVs")
    end
    @time(
        for i1 in 1:lrprm1
            for i2 in 1:lrprm2
                try
                    print(recursive_string,rprm1[i1], " / ", rprm2[i2])
                    
                    ##########################################
                    beta = rprm1[i1]
                    A = rprm2[i2]
                    ##########################################
                    
                   if EIS_bool
                        CV_matrix[i1,i2] = ysz_experiments.run_new(
                            pyplot=false, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias, omega_range=EIS_get_shared_omega_range(),
                            dx_exp=-8,
                            T=TCtoT(TC), pO2=pO2tosim(pO2),
                            prms_in=[A0, R0, DGA, DGR, beta, A],
                            add_prms_in=(DD, nu, nus, ms_par)
                        )
                    else
                        CV_matrix[i1,i2] = ysz_experiments.run_new(
                            out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
                            prms_in=[A0, R0, DGA, DGR, beta, A],
                            add_prms_in=(DD, nu, nus, ms_par)
                        )
                    end
                    ok_matrix[i1, i2] = 1
                    println("  .. ok :)")
                catch e
                    if e isa InterruptException
                        rethrow(e)
                    else
                        ok_matrix[i1, i2] = 0
                        println("  :(   ")
                        continue
                    end
                end
            end
        end
    )
    
    #PyPlot.figure(figsize=(6,4))
    #xlabel("parameter")
    #zlabel("error")
    #title("Fitness Function")

    if pyplot && !(recursive)
        figure(figsize=(3,3))
    end
    #print_ok_matrix(ok_matrix)
    
    println(recursive_string,"Computing interpolation ... ")
    @time(
        for i1 in 1:lrprm1-1
            for i2 in 1:lrprm2-1
                if all_four_nodes_ok(ok_matrix,i1,i2)
                    print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
                    
                    if EIS_bool
                         FF=EIS_get_FF_interp_2D(
                            [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
                            rprm1[i1:i1+1],
                            rprm2[i2:i2+1];
                            nx=nx, ny=ny,
                            TC=TC, pO2=pO2, EIS_bias=EIS_bias
                        )                       
                    else
                        FF=CV_get_FF_interp_2D(
                            [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
                            rprm1[i1:i1+1],
                            rprm2[i2:i2+1];
                            nx=nx, ny=ny,
                            TC=TC, pO2=pO2
                        )
                    end
                    
                    min_count, min_values, min_positions, min_curv, hard_min = find_mins(
                        FF, 
                        [(rprm1[i1+1] - rprm1[i1])/nx,
                         (rprm2[i2+1] - rprm2[i2])/ny]
                    )
                    print("... h_min = ",hard_min)
                    
                    
                    if min_count > 0
                        println(" ... min_c ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
                        #println(" ... min ", min_count)#
                        
                        approx_min, approx_min_position = pickup_min(min_count, min_values, min_positions)
                        
                        if depth_remaining > 0
                            println(recursive_string,">> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
                            new_min, new_min_position = scan_2D_recursive(;pyplot=true, 
                                TC=TC, pO2=pO2,
                                rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
                                rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
                                approx_min=approx_min,
                                approx_min_position=approx_min_position,
                                recursive=true,
                                depth_remaining=depth_remaining - 1,
                                recursive_string=string(recursive_string, "<", depth_remaining - 1, ">"),
                                nx = Int32(round(nx/1.2)), ny = Int32(round(ny/1.2))
                            )
                            if new_min < global_min
                                global_min = new_min
                                global_min_position =new_min_position
                            end
                        else
                            if approx_min < global_min
                                global_min = approx_min
                                global_min_position = approx_min_position
                            end
                            if pyplot
                                surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
                            end
                        end
                    else
                        if pyplot
                            surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
                        end
                        println()
                    end
                end 
            end
        end
    )
    
    
    
    if pyplot && !(recursive)
        if EIS_bool 
            plot_EIS_matrix(CV_matrix, ok_matrix, rprm1, rprm2, TC, pO2, EIS_bias)
        else
            plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, TC, pO2)
        end
    end
    return global_min, global_min_position
end


# # # # function scan_2D_recursive(;pyplot=false, 
# # # #                             TC=800, pO2=100,
# # # #                             #rprm1=range(15, stop=20, length=3),
# # # #                             #rprm2=range(15, stop=20, length=3),
# # # #                             #rprm1=range(18.4375, stop=18.59375, length=3),
# # # #                             #rprm2=range(20.59375, stop=21, length=3),
# # # #                             rprm1=range(10.0, stop=28., length=4),
# # # #                             #rprm2=range(15.0, stop=23, length=7),
# # # #                             rprm2=range(-2.5, stop=1.0, length=4),
# # # #                             nx=100, ny=100,
# # # #                             wp=[1, 1, 0, 0, 0, 0], #TODO !!!
# # # #                             depth_remaining=1,
# # # #                             approx_min=Inf,
# # # #                             approx_min_position=[0, 0],
# # # #                             recursive=false,
# # # #                             recursive_string="",
# # # #                             )
# # # # 
# # # #     println("nx = ",nx, " ..... ny = ",ny)
# # # #     #PyPlot.close()
# # # #     #PyPlot.close()
# # # #     
# # # #     ######################################################
# # # #     #rprm1 = collect(-5.0 : 1.0 : 5.0)
# # # #     #rprm1 = collect(-3. : 0.5 : 3.0)
# # # #     #rprm1=range(15, stop=20, length=3)
# # # #     
# # # #     #rprm2 = collect(-3. : 0.5 : 3.0)
# # # #     #rprm2=range(15, stop=20, length=3)
# # # #     ######################################################
# # # #     
# # # #     wpn = ["A0","R0","DGA","DGR","betaR","SR"]
# # # #     #A0 = 21.71975544711280
# # # #     #R0 = 19.53
# # # #     A0 = 19.50
# # # #     R0 = 19.85
# # # #     DGA = 0.0905748
# # # #     DGR = -0.708014
# # # #     beta = 0.6074566741435283
# # # #     A = -0.356
# # # # 
# # # #     println(recursive_string,"computing 2D scan... ")
# # # #     lrprm1 = size(rprm1,1)
# # # #     lrprm2 = size(rprm2,1)
# # # #     
# # # #     global_min = Inf
# # # #     global_min_position = [0,0]
# # # #     
# # # #     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
# # # #     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
# # # #     println(recursive_string,"Computing CVs")
# # # #     @time(
# # # #         for i1 in 1:lrprm1
# # # #             for i2 in 1:lrprm2
# # # #                 try
# # # #                     print(recursive_string,rprm1[i1], " / ", rprm2[i2])
# # # #                     
# # # #                     ##########################################
# # # #                     R0 = rprm1[i1]
# # # #                     A = rprm2[i2]
# # # #                     ##########################################
# # # #                     
# # # #                    if true
# # # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # # #                             out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
# # # #                             prms_in=[A0, R0, DGA, DGR, beta, A]
# # # #                         )
# # # #                     else
# # # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # # #                             pyplot=false, EIS_IS=true, out_df_bool=false, 
# # # #                             tref=0,
# # # #                             prms_in=[A0, R0, DGA, DGR, beta, A],  
# # # #                             nu_in=0.9, pO2_in=1.0, DD_in=9.5658146540360312e-10
# # # #                         )
# # # #                     end
# # # #                     ok_matrix[i1, i2] = 1
# # # #                     println("  .. ok :)")
# # # #                 catch e
# # # #                     if e isa InterruptException
# # # #                         rethrow(e)
# # # #                     else
# # # #                         ok_matrix[i1, i2] = 0
# # # #                         println("  :(   ")
# # # #                         continue
# # # #                     end
# # # #                 end
# # # #             end
# # # #         end
# # # #     )
# # # #     
# # # #     #PyPlot.figure(figsize=(6,4))
# # # #     #xlabel("parameter")
# # # #     #zlabel("error")
# # # #     #title("Fitness Function")
# # # # 
# # # #     if pyplot && !(recursive)
# # # #         figure(figsize=(3,3))
# # # #     end
# # # #     #print_ok_matrix(ok_matrix)
# # # #     
# # # #     println(recursive_string,"Computing interpolation ... ")
# # # #     @time(
# # # #         for i1 in 1:lrprm1-1
# # # #             for i2 in 1:lrprm2-1
# # # #                 if all_four_nodes_ok(ok_matrix,i1,i2)
# # # #                     print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
# # # #                     
# # # #                     FF=get_FF_interp_2D(
# # # #                         [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
# # # #                         rprm1[i1:i1+1],
# # # #                         rprm2[i2:i2+1];
# # # #                         nx=nx, ny=ny,
# # # #                         TC=TC, pO2=pO2
# # # #                     )
# # # #                     
# # # #                     min_count, min_values, min_positions, min_curv, hard_min = find_mins(
# # # #                         FF, 
# # # #                         [(rprm1[i1+1] - rprm1[i1])/nx,
# # # #                          (rprm1[i2+1] - rprm1[i2])/ny]
# # # #                     )
# # # #                     print("... h_min = ",hard_min)
# # # #                     
# # # #                     
# # # #                     if min_count > 0
# # # #                         println(" ... min_c ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
# # # #                         #println(" ... min ", min_count)#
# # # #                         
# # # #                         approx_min, approx_min_position = pickup_min(min_count, min_values, min_positions)
# # # #                         
# # # #                         if depth_remaining > 0
# # # #                             println(recursive_string,">> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
# # # #                             new_min, new_min_position = scan_2D_recursive(;pyplot=true, 
# # # #                                 TC=TC, pO2=pO2,
# # # #                                 rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
# # # #                                 rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
# # # #                                 approx_min=approx_min,
# # # #                                 approx_min_position=approx_min_position,
# # # #                                 recursive=true,
# # # #                                 depth_remaining=depth_remaining - 1,
# # # #                                 recursive_string=string(recursive_string, "<", depth_remaining - 1, ">"),
# # # #                                 nx = Int32(round(nx/1.2)), ny = Int32(round(ny/1.2))
# # # #                             )
# # # #                             if new_min < global_min
# # # #                                 global_min = new_min
# # # #                                 global_min_position =new_min_position
# # # #                             end
# # # #                         else
# # # #                             if approx_min < global_min
# # # #                                 global_min = approx_min
# # # #                                 global_min_position = approx_min_position
# # # #                             end
# # # #                             if pyplot
# # # #                                 surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
# # # #                             end
# # # #                         end
# # # #                     else
# # # #                         if pyplot
# # # #                             surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
# # # #                         end
# # # #                         println()
# # # #                     end
# # # #                 end 
# # # #             end
# # # #         end
# # # #     )
# # # #     
# # # #     
# # # #     
# # # #     if pyplot && !(recursive)
# # # #         plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, TC, pO2)
# # # #     end
# # # #     return global_min, global_min_position
# # # # end



# function scan_2D(;pyplot=false, 
#                 TC=800, pO2=100,
# #                rprm1=range(15, stop=20, length=3),
# #                rprm2=range(15, stop=20, length=3),
#                 #rprm1=range(18.4375, stop=18.59375, length=3),
#                 #rprm2=range(20.59375, stop=21, length=3),
#                 rprm1=range(18.4375, stop=18.59375, length=3),
#                 rprm2=range(20.59375, stop=21, length=3),
#                 wp=[1, 1, 0, 0, 0, 0], #TODO !!!
#                 depth_remaining=2,
#                 recursive=false,
#                 nx=100, ny=100
#                 )
# 
#     
#     #PyPlot.close()
#     #PyPlot.close()
#     
#     ######################################################
#     #rprm1 = collect(-5.0 : 1.0 : 5.0)
#     #rprm1 = collect(-3. : 0.5 : 3.0)
#     #rprm1=range(15, stop=20, length=3)
#     
#     #rprm2 = collect(-3. : 0.5 : 3.0)
#     #rprm2=range(15, stop=20, length=3)
#     ######################################################
#     
#     wpn = ["A0","R0","DGA","DGR","betaR","SR"]
#     A0 = 21.71975544711280
#     R0 = 19.53
#     DGA = 0.0905748
#     DGR = -0.708014
#     beta = 0.6074566741435283
#     A = -0.356
# 
#     println("computing 2D scan... ")
#     lrprm1 = size(rprm1,1)
#     lrprm2 = size(rprm2,1)
#     
#     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
#     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
#     println("Computing CVs")
#     @time(
#         for i1 in 1:lrprm1
#             for i2 in 1:lrprm2
#                 try
#                     print(rprm1[i1], " / ", rprm2[i2])
#                     
#                     ##########################################
#                     A0 = rprm1[i1]
#                     R0 = rprm2[i2]
#                     ##########################################
#                     
#                     CV_matrix[i1,i2] = ysz_experiments.run_new(
#                         out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
#                         prms_in=[A0, R0, DGA, DGR, beta, A]
#                         )
#                     ok_matrix[i1, i2] = 1
#                     println("  .. ok :)")
#                 catch e
#                     if e isa InterruptException
#                         rethrow(e)
#                     else
#                         ok_matrix[i1, i2] = 0
#                         println("  :(   ")
#                         continue
#                     end
#                 end
#             end
#         end
#     )
#     
#     #PyPlot.figure(figsize=(6,4))
#     #xlabel("parameter")
#     #zlabel("error")
#     #title("Fitness Function")
# 
#     if pyplot && !(recursive)
#         figure(figsize=(3,3))
#     end
#     print_ok_matrix(ok_matrix)
#     
#     println("Computing interpolation ... ")
#     @time(
#         for i1 in 1:lrprm1-1
#             for i2 in 1:lrprm2-1
#                 if all_four_nodes_ok(ok_matrix,i1,i2)
#                     print(" intrp i1 / i2 : ",i1, " / ", i2)
#                     
#                     FF=get_FF_interp_2D(
#                         [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
#                         rprm1[i1:i1+1],
#                         rprm2[i2:i2+1];
#                         nx=nx, ny=ny,
#                         TC=TC, pO2=pO2
#                     )
#                     
#                     min_count, min_values, min_positions, min_curv = find_mins(
#                         FF.err_matrix, 
#                         [(rprm1[i1+1] - rprm1[i1])/nx,
#                          (rprm1[i2+1] - rprm1[i2])/ny]
#                     )
#                     
# 
#                     
#                     if pyplot
#                         println("jsem tu")
#                         surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
#                     end
#                     
#                     if min_count > 0
#                         print(" ... min ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
#                         
#                         if depth_remaining > 0
#                             println("\n >>>>>>>> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
#                             scan_2D(;pyplot=true, 
#                                 TC=TC, pO2=pO2,
#                                 rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
#                                 rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
#                                 depth_remaining=depth_remaining - 1,
#                                 recursive=true
#                                 )
#                         end
#                     end
#                     
#                     print("\n")
#                 end 
#             end
#         end
#     )
#     
#     
#     
#     if pyplot && !(recursive)
#         plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, TC, pO2)
#     end
#     return
# end



function del_pyplot()
    for i in 1:10
        PyPlot.close()
    end
end










#####################
## 1D optimization ##
#####################
function get_FF_interp_1D(CV_list, prms_list; n=20)
    CV_raw = import_CVtoDataFrame(;TC=800,pO2=100)
    checknodes = CV_get_checknodes(0.05,0.95,-0.95,-0.05,0.01)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_raw, checknodes))

    CV1 = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[1],checknodes))
    CV2 = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[end],checknodes))
    node_list = range(0.0, stop=1.0, length=n)
    err_list_interp = [CV_fitnessFunction(
            CV_exp,
            DataFrame(
                U = checknodes[:,1], 
                I = (linComb_CVs(CV1, CV2, i)).I
            )
        )
        for i in node_list
    ]
    prms_list_interp = node_list.*(prms_list[end]-prms_list[1]).+prms_list[1]
    FF_out = DataFrame(prms = prms_list_interp, err = err_list_interp)
end

function comparison()
    CV_list = []
    prms_list = []
    PyPlot.close()
    
    #for A in collect(-0.36: 0.1 : -0.32)
    #for A in [-0.345, -0.33]
    for A in [-0.356, -0.354]
    #for A in [-0.1, -0.6]
    #for A in collect(-0.1 : -0.01 : -0.6)
    #for A in [-0.3893]
    #for R0 in [21.4, 21.8]
    #for R0 in collect(21.4: 0.02 : 21.8)
    for R0 in [20.53]
        println(A)
        ww=DataFrame()
        push!(CV_list,ysz_experiments.run_new(
            out_df_bool=true, voltammetry=true, sample=10,
            prms_in=[21.71975544711280, R0, 0.0905748, -0.708014, 0.6074566741435283, A]
            )
        )
        push!(prms_list,A)
    end
    end
    println("jdeme dal")
    plot_all(CV_list, prms_list)
    plot_error(CV_list, prms_list)
    return
end


function scan_1D()
    PyPlot.close()
    PyPlot.close()
    #rA = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.2]    
    #rA = collect(-0.4 : 0.02 : -0.3)
    
    #rprms = [19, 20, 21, 22]
    rprms = collect(-0.5 : 0.1 : 0.0)
    
    
    A0 = 21.71975544711280
    R0 = 20.606423236896422
    DGA = 0.0905748
    DGR = -0.708014
    beta = 0.6074566741435283
    A = -0.356
    
    println("computing ... ")
    CV_list = []
    ok_prms_list = []
    for i in 1:size(rprms,1)
        try
        print(rprms[i])
        
        ########
        A = rprms[i]
        ########
        
        push!(CV_list,ysz_experiments.run_new(
            out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
            prms_in=[A0, R0, DGA, DGR, beta, A]
            )
        )
        push!(ok_prms_list,rprms[i])
        println("  .. ok")
        catch e
            if e isa InterruptException
                rethrow(e)
            else
                println()
                continue
            end
        end
    end

    PyPlot.figure(figsize=(6,4))
    xlabel("parameter")
    ylabel("error")
    title("Fitness Function")
    FF = DataFrame()
    for i in 1:(size(ok_prms_list,1)-1)
        FF = get_FF_interp_1D([CV_list[i], CV_list[i+1]],[ok_prms_list[i], ok_prms_list[i+1]])
        plot(FF.prms, FF.err, linewidth=3.0, label=string("FF ",i))
        plot(FF.prms[1],FF.err[1], "kx", markersize=12.0, ) 
        #plot(FF.prms, FF.err)
    end
    plot(FF.prms[end],FF.err[end], "kx", markersize=12.0)
    legend(loc="best")
    grid()
    
    plot_all(CV_list, ok_prms_list)
    
    return
end



####################################
####################################
####################################
######## Levenberg-Maquard ##########
####################################

function LM_optimize(;EIS_opt_bool=false, CV_opt_bool=false, pyplot=false)
    function prepare_prms(mask, x0, x)
        prms = zeros(0)
        xi = 1
        for i in collect(1 : 1 :length(mask))
            if convert(Bool,mask[i])
                append!(prms, x[xi])
                xi += 1
            else
                append!(prms, x0[i])
            end
        end
        return prms
    end

   function to_optimize(x) 
        #err = run_new(print_bool=false, fitting=true, voltametry=true, pyplot=false, voltrate=0.005, sample=8, bound=0.41, 
        #    prms_in=x)
        prms = prepare_prms(mask, x0, x)
        print(" >> mask = ",mask)
        print(" || prms = ",prms)
        
        prms_names_in=[]
        prms_values_in=[]
        
        append!(prms_names_in, prms_names)
        append!(prms_values_in, prms)
      
       
        err = 0.0
        err_string = ""
        if CV_opt_bool
            CV_sim = ysz_experiments.run_new(
                        physical_model_name=physical_model_name,
                        out_df_bool=true, voltammetry=true, sample=8, pyplot=false,
                        T=TCtoT(TC), pO2=pO2tosim(pO2),
                        prms_names_in=prms_names_in,
                        prms_values_in=prms_values_in,
            )
            if pyplot
                figure(1)
                CV_plot(CV_sim)
            end
            
            err +=CV_penalty_factor * CV_fitnessFunction(
                CV_apply_checknodes(CV_sim, CV_get_shared_checknodes()), 
                CV_exp
            )
            err_string = string(err_string," CV ")
        end
        if EIS_opt_bool
            EIS_sim = ysz_experiments.run_new( 
                    physical_model_name=physical_model_name,
                    pyplot=false, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias, omega_range=EIS_get_shared_omega_range(),
                    dx_exp=-9,
                    T=TCtoT(TC), pO2=pO2tosim(pO2),
                    prms_names_in=prms_names_in,
                    prms_values_in=prms_values_in,
            )

            
            if pyplot
                figure(2)
                Nyquist_plot(EIS_sim)
            end
            
            err += EIS_fitnessFunction(EIS_sim, EIS_exp)
            err_string = string(err_string," EIS ")
        end
        println(" || ",err_string,": err =", err)
        return [err]
    end
    

    
    TC=800
    pO2 = 100
    EIS_bias=0.0
    physical_model_name="necum"
    CV_penalty_factor = 10  # for fitnessFunction = factor*CV + EIS 
    
    
    if pyplot
        PyPlot.close()
    end
    
    if CV_opt_bool
        CV_exp = CV_apply_checknodes(
            import_CVtoDataFrame(TC=TC, pO2=pO2),
            CV_get_shared_checknodes()
        )
        if pyplot
            figure(1)
            CV_plot(CV_exp, "exp $(CV_experiment_legend(TC, pO2))")
        end
    end
    
    if EIS_opt_bool
        EIS_exp = EIS_apply_checknodes(
            import_EIStoDataFrame(TC=TC,pO2=pO2,bias=EIS_bias),
            EIS_get_shared_checknodes()
        )

        if pyplot
            figure(2)
            Nyquist_plot(EIS_exp, "exp $(EIS_experiment_legend(TC, pO2, EIS_bias))")
        end
    end
            
    #x0 = zeros(2)
    #optimize(rosenbrock, x0, LevenbergMarquardt())
    
    prms_names=["A0", "R0", "K0", "SA", "SR", "SO", "DGA", "DGR", "DGO", "betaA", "betaR", "betaO", "DD", "nu", "nus", "ms_par"]
    
    lower_bounds      = [17, 17, 17,     -5, -5, -5,     -0.9, -0.9, -0.9,    0.0, 0.0, 0.0,     0.1e-13, 0.1, 0.1, 0.01 ]
    upper_bounds      = [23, 23, 23,      1,  1,  1,      0.9,  0.9,  0.9,    1.0, 1.0, 1.0,     9.0e-13, 0.9, 0.9, 1.5 ]
    
    x0  = [19.7, 19.7, 18.6,    1, 1, 1,    0.7, -0.8, -0.3,      0.5, 0.5, 0.5,     5.35e-13,  0.85,  0.21, 0.05]
    #x0 = [18.673485702096716, 18.598123074869402, 17.50858588747129, 1.0, 1.0, 1.0, 0.8666319406487695, -0.8342275189659124, -0.7728570698687104, 0.5, 0.5, 0.5, 5.35e-13] 
    
    
    mask = [1, 1, 1,     0, 0, 0,     1,1,1,      0,0,0,   1, 1, 1, 1] # determining, which parametr should be fitted
    
    x0M = zeros(0)
    lowM = zeros(0)
    uppM = zeros(0)
    for i in collect(1 : 1 : length(mask))
        if convert(Bool,mask[i])
            append!(x0M, x0[i])
            append!(lowM, lower_bounds[i])
            append!(uppM, upper_bounds[i])
        end
    end
    
   
    
    #to_optimize(x0M)
    println(optimize(to_optimize, x0M, lower=lowM, upper=uppM, Δ=1000000, f_tol=1.0e-14, g_tol=1.0e-14, LevenbergMarquardt()))
    
    ####optimize(to_optimize, x0M, Dogleg())
    return
end


###############################
###############################
#### SNEHURKA COMPUTATION #####
###############################

function run_par_study(actual_par_study::par_study_struct; save_dir="../snehurka/data/", save_file_bool=false, mode="test")
     
     
  function consistency_check()
      
      
      # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # (ale az v pripade, ze to bude obecne napsane)
      
      return true
  end
  SIM_length = size(actual_par_study.simulation_list,1)
  
  SIM_err_counter = zeros(Int64,SIM_length)
  SIM_good_counter = zeros(Int64,SIM_length)
  all_counter::Int64 = 0
  
  save_file_path = "$(save_dir)$(actual_par_study.name)/"

  

  #println(" >> number of sets of parameters: ",
  #    length(A0_list)*length(R0_list)*length(DGA_list)*length(DGR_list)*length(beta_list)*length(A_list))
    
  function perform_par_study(actual_prms)
    prms = actual_prms
    prms_names = actual_par_study.prms_names
    
    prms_names_in=[]
    prms_values_in=[]
    
    append!(prms_names_in, prms_names)
    append!(prms_values_in, prms)
    
    string_prms_names = ""
    string_prms = ""
    overall_count = 1
    for i in 1:size(prms,1)
      string_prms_names = "$(string_prms_names)$(prms_names[i]), "
      string_prms = "$(string_prms)$(prms[i]), "
      overall_count *= size(actual_par_study.prms_lists[i],1)
    end
    string_prms_names = string_prms_names[1:end-2]
    string_prms = string_prms[1:end-2]

    # printing info
    all_counter = all_counter + 1
    println(
      string(" <><><><><><><><><><><><> all_counter <><><> calculating ",all_counter," of ", overall_count)
    )
    print("parameters: ($(string_prms_names)) = ($(string_prms))")
    
    for (i,SIM) in enumerate(actual_par_study.simulation_list)
        try

            SIM_sim = typical_run_simulation(SIM, prms_names_in, prms_values_in)
            
            if save_file_bool
              save_file_prms(SIM, SIM_sim, "$(save_file_path)$(string(SIM))/", prms_values_in, prms_names_in, actual_par_study.scripted_tuple)           
            end
            SIM_good_counter[i] += 1
            print(" $(string(SIM)) ok :) ")
          catch e
            if e isa InterruptException
                rethrow(e)
            else
                println(e)
                SIM_err_counter[i] += 1
                print("<<<<< $(string(SIM))  FAIL! >>>>>")
            end
        end        
    end
    println()    
  end
      
  for_each_prms_in_prms_lists(actual_par_study.prms_lists, perform_par_study)

  println(string("<<<<< SIM_good_count / SIM_err_count  >>> ", SIM_good_counter,"  /  ", SIM_err_counter ," >>>>>>>"))

end



function run_par_study_script_wrap(
                    prms_lists_string=("[20, 20, 20, 0.0, 0.0, 0.0]"),
                    save_dir="./kadinec/",
                    name="default_par_study_name",
                    scripted_tuple_string=("(1, 0, 0, 0, 0, 0)"),
                    prms_names_string="[\"A0\", \"R0\", \"K0\", \"DGA\", \"DGR\", \"DGO\"]", 
                    TC_string="800",
                    pO2_string="100",
                    EIS_bias_string="0.0",
                    CV_bool_string="true", 
                    EIS_bool_string="true", 
                    mode="test",
                    physical_model_name="nothing")
  
  
  actual_par_study = par_study_struct()
  
  actual_par_study.physical_model_name = physical_model_name
  actual_par_study.name = name
  actual_par_study.prms_lists = eval(Meta.parse(prms_lists_string))
  actual_par_study.scripted_tuple = eval(Meta.parse(scripted_tuple_string))
  actual_par_study.prms_names = eval(Meta.parse(prms_names_string))
  #
  TC = eval(Meta.parse(TC_string))
  pO2 = eval(Meta.parse(pO2_string))
  EIS_bias = eval(Meta.parse(EIS_bias_string))
  CV_bool = eval(Meta.parse(CV_bool_string))
  EIS_bool = eval(Meta.parse(EIS_bool_string)) 
  #
  if CV_bool && EIS_bool
    actual_par_study.simulation_list = [
      CV_simulation(TC, pO2)... ,
      EIS_simulation(TC, pO2, EIS_bias)...
    ]
  elseif CV_bool
    actual_par_study.simulation_list = [
      CV_simulation(TC, pO2)...
    ]
  elseif EIS_bool
    actual_par_study.simulation_list = [
      EIS_simulation(TC, pO2, EIS_bias)...
    ]
  else
    actual_par_study.simulation_list = Array{abstract_simulation}[]
  end
  
  run_par_study(    actual_par_study,
                    save_dir=save_dir, 
                    save_file_bool=true,
                    mode=mode)
end


function meta_run_par_study()  
  function recursive_bash_call(output_prms_lists, active_idx)
    if active_idx > size(scripted_tuple,1)
      scripted_tuple_string = string(scripted_tuple)
      output_prms_lists_string = string(output_prms_lists)
      prms_names_string = string(prms_names)
      CV_bool_string = string(CV_bool)
      EIS_bool_string = string(EIS_bool)
      TC_string = string(TC)
      pO2_string = string(pO2)
      EIS_bias_string = string(EIS_bias)
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
                  $EIS_bias_string 
                  $CV_bool_string 
                  $EIS_bool_string 
                  $mode 
                  $physical_model_name 
      `)
    end
    if scripted_tuple[active_idx] == 1
      for i in prms_lists[active_idx]
        recursive_bash_call(push!(deepcopy(output_prms_lists),i), active_idx + 1)
      end
    else
      recursive_bash_call(push!(deepcopy(output_prms_lists),prms_lists[active_idx]), active_idx + 1)
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
  # par_study definition ####################################
  
  physical_model_name = "ysz_model_GAS_LoMA"
  name = "GAS_LoMA_ULTRA_level_2"
  
  prms_names = ["A0", "R0", "K0", "DGA", "DGR", "DGO", "DD"]
  prms_lists = [
    collect(18.5 : 0.5 : 21.5),  
    collect(18.5 : 0.5 : 21.5),  
    collect(18.5 : 0.5 : 21.5), 
    collect(-0.7 : 0.35 : 0.7), 
    collect(-0.7 : 0.35 : 0.7), 
    collect(-0.7 : 0.35 : 0.7),
    # hint: TC = (700, 750, 800, 850)  => DD = ( 1.277, 2.92, 5.35, 9.05)e-13
    [9.05e-13]
  ]
  scripted_tuple = (1, 1, 0, 0, 0, 0, 0)
  
  TC = 850
  pO2 = [60, 80]
  EIS_bias = 0.0
  
  EIS_bool = true
  CV_bool = true
  
  # TODO !!!
#   simulation_list = [
#     CV_simulation(TC, pO2)... ,
#     EIS_simulation(TC, pO2, EIS_bias)...
#   ] 
  
  #######################################################
  
  # preparing bash output ###############################
  #bash_command = "sbatch"
  #bash_command = "echo"
  bash_command = "julia"
  
  save_dir = "../snehurka/data/"

  
  mode = "test_one_prms"
  #mode = "only_print"
  #mode = "go"
  
  run_file_name = "../snehurka/run_ysz_fitting_par_study-prms-.jl"
  #######################################################
  #######################################################

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

  metafile_string = "METAFILE for par_study\n"
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
  metafile_string = string(metafile_string,"EIS_bias=", EIS_bias,"\n")
  metafile_string = string(metafile_string,"CV_bool=", CV_bool,"\n")
  metafile_string = string(metafile_string,"EIS_bool=", EIS_bool,"\n")
  #metafile_string = string(metafile_string,"string(simulation_list)=", string(simulation_list),"\n") 
  metafile_string = string(metafile_string,"name=", name,"\n")
  metafile_string = string(metafile_string,"save_dir=", save_dir,"\n")
  metafile_string = string(metafile_string,"mode=", mode,"\n")
  metafile_string = string(metafile_string,"bash_command=", bash_command,"\n")
  metafile_string = string(metafile_string,"run_file_name=", run_file_name,"\n")
  metafile_string = string(metafile_string,"prms_lists=", prms_lists,"\n")
  #metafile_string = string(metafile_string,"simulation_list=", simulation_list,"\n")
  
  
  
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
