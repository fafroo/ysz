using CSV
using DataFrames
using PyPlot

############################################################
mutable struct par_study_struct
  physical_model_name::String
  SIM_list::Array{abstract_simulation}
  prms_names::Array{Any}
  prms_lists::Array{Any}
  scripted_tuple::Tuple
  SIMs_data::Array{Any}
  SIMs_err::DataFrame
  #
  name::String
  save_dir::String
  #
  par_study_struct() = new()
end

# function par_study()
# end


#############################################################


############## import large amount of data ##########
function par_study_import_info_from_metafile(name; save_dir="../data/par_studies/", metafile_name="__metafile_par_study.txt")

  actual_par_study = par_study_struct()
  
  actual_par_study.save_dir = save_dir
  
  TC = Inf
  pO2 = Inf
  bias = Inf
  
  simulations = []
  
  file_lines = readlines(string(save_dir,name,"/",metafile_name))
  for str in file_lines
    if occursin('=',str)
      (token_name,token_value)=split(str,"=")
      token_name=="physical_model_name" && (actual_par_study.physical_model_name = token_value)
      token_name=="prms_names" && (actual_par_study.prms_names = eval(Meta.parse(token_value)))
      token_name=="prms_lists" && (actual_par_study.prms_lists = eval(Meta.parse(token_value)))
      token_name=="scripted_tuple" && (actual_par_study.scripted_tuple = eval(Meta.parse(token_value)))
      token_name=="name" && (actual_par_study.name = token_value)
      #
      token_name=="TC" && (TC = eval(Meta.parse(token_value)))
      token_name=="pO2" && (pO2 = eval(Meta.parse(token_value)))
      token_name=="bias" && (bias = eval(Meta.parse(token_value)))
      token_name=="simulations" && (simulations = eval(Meta.parse(token_value)))
    end
  end
  
  if (TC == Inf) || (pO2 == Inf) || (bias == Inf)
    println("ERROR: par_study_import_info_from_metafile: load error of TC, pO2 or bias!")
    throw(Exception)
  else
    actual_par_study.SIM_list = get_SIM_list_rectangle(TC, pO2, bias, simulations)
  end

  return actual_par_study
end


function par_study_import_data!(actual_par_study; verbose=false)
  
  
  save_dir_forwarded = string(actual_par_study.save_dir, actual_par_study.name, "/")
  
  SIMs_length = size(actual_par_study.SIM_list,1)
  SIMs_dimensions = Tuple(append!([SIMs_length], [size(list,1) for list in actual_par_study.prms_lists]))
  
  SIMs_data = Array{DataFrame}(undef, SIMs_dimensions)
  SIM_good_count = zeros(Int64, size(actual_par_study.SIM_list,1))
  SIM_bad_count = zeros(Int64, size(actual_par_study.SIM_list,1))
  
  
  empty_df = DataFrame()
  
  counter=1
  circle_counter = 0
  for list in actual_par_study.prms_lists
    counter*= size(list,1)
  end
  prmsDIV100 = counter/100
  
  counter=0
  println("Loading prms_lists")
  println("|====================================================================================================|")

  function perform_import(prms_indicies)
    
    counter += 1
    while (circle_counter + 1)*prmsDIV100 <= counter + 0.00000001
      circle_counter+=1
      print("o")
    end
    
    for (i, SIM) in enumerate(actual_par_study.SIM_list)
      try
        SIMs_data[i, prms_indicies...] = load_file_prms(
            SIM, 
            save_dir=string(save_dir_forwarded,string(SIM),"/"), 
            prms=get_prms_from_indicies(actual_par_study.prms_lists, prms_indicies), 
            prms_names=actual_par_study.prms_names, 
            scripted_tuple=actual_par_study.scripted_tuple, 
            throw_exception=true, 
            verbose=verbose)
        SIM_good_count[i] += 1
      catch e
        if e isa InterruptException
          rethrow(e)
        else
          SIMs_data[i, prms_indicies...] = empty_df
          SIM_bad_count[i] += 1
        end
      end
    end
  end
  
  print("|")
  for_each_indicies_in_prms_lists(actual_par_study.prms_lists, perform_import)
  println("|")
  
  actual_par_study.SIMs_data = SIMs_data
  println(string(actual_par_study.SIM_list))
  println("par_study_import_data: SIM_good / all = ",SIM_good_count, " / ", SIM_good_count .+ SIM_bad_count)
  return
end


function par_study_err_evaluation!(actual_par_study; data_already_imported=false, verbose=false)
  save_dir_forwarded = string(actual_par_study.save_dir, actual_par_study.name, "/")
 
  # compute length of error DataFrame
  product_size_prms=1
  for list in actual_par_study.prms_lists
    product_size_prms*= size(list,1)
  end
  prmsDIV100 = product_size_prms/100
  SIMs_length = size(actual_par_study.SIM_list,1)
  
  # SIMs_err stores errors in columns (1, 2, 3, 4, .. SIMs_length, prms, total_error)
  SIMs_err = DataFrame()
  for i in 1:SIMs_length
    SIMs_err[!, Symbol(string(i))] = Array{Float32}(undef, product_size_prms)
  end
  SIMs_err[!, Symbol("prms_indicies")] = Array{Array{Int16}}(undef, product_size_prms)
  SIMs_err[!, Symbol("error_total")] = Array{Float32}(undef, product_size_prms)

  # experimental data to compare
  SIM_list_exp = []
  for (i, SIM) in enumerate(actual_par_study.SIM_list)
    push!(SIM_list_exp, apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes))
  end  

  SIM_good_count = zeros(Int64, size(actual_par_study.SIM_list,1))
  SIM_bad_count = zeros(Int64, size(actual_par_study.SIM_list,1))

  all_errors_sum::Float64 = 0
  aux::Float64 = 0
  
  counter::Int64 = 0
  circle_counter::Int64 = 0
  println("Evaluating error")
  println("|====================================================================================================|")

  function perform_error_evaluate(prms_indicies)
    
    counter += 1
    while (circle_counter + 1)*prmsDIV100 <= counter + 0.00000001
      circle_counter+=1
      print("o")
    end
    SIMs_err[!, SIMs_length + 1][counter] = prms_indicies

    all_errors_sum = 0
      
    for (i, SIM) in enumerate(actual_par_study.SIM_list)
      if data_already_imported        
        # load SIM from par_study_struct
        if size(actual_par_study.SIMs_data[i, prms_indicies...],1) > 0
          SIM_sim = actual_par_study.SIMs_data[i, prms_indicies...]
          
          # compute fitness function
          error = fitnessFunction(SIM, 
                SIM_list_exp[i],
                apply_checknodes(SIM, SIM_sim, SIM.checknodes)
              )
          all_errors_sum += error * SIM.fitness_factor
          SIMs_err[!, i][counter] = error
          
          SIM_good_count[i] += 1
        else
          all_errors_sum += Inf
          SIMs_err[!, i][counter] = Inf
          SIM_bad_count[i] += 1
        end
      else
        # try load SIM from a file
        try
          SIM_sim = load_file_prms(
            SIM, 
            save_dir=string(save_dir_forwarded, string(SIM), "/"), 
            prms=get_prms_from_indicies(actual_par_study.prms_lists, prms_indicies), 
            prms_names=actual_par_study.prms_names, 
            scripted_tuple=actual_par_study.scripted_tuple, 
            throw_exception=true, 
            verbose=verbose)
          
          # compute fitness function
          error = fitnessFunction(SIM, 
                SIM_list_exp[i],
                apply_checknodes(SIM, SIM_sim, SIM.checknodes)
              )
          all_errors_sum += error * SIM.fitness_factor
          SIMs_err[!, i][counter] = error

          SIM_good_count[i] += 1
        catch e
          if e isa InterruptException
            rethrow(e)
          else
            all_errors_sum += Inf
            SIMs_err[!, i][counter] = Inf
            SIM_bad_count[i] += 1
          end
        end
      end
      
    end
    
    SIMs_err[!, SIMs_length + 2][counter] = all_errors_sum/SIMs_length
    
    
  end
  
  print("|")
  for_each_indicies_in_prms_lists(actual_par_study.prms_lists, perform_error_evaluate)
  println("|")
  
  sort!(SIMs_err, SIMs_length + 2)
  
  actual_par_study.SIMs_err = SIMs_err
  
  println("par_study_err_evaluation: SIM_good / all = ",SIM_good_count, " / ", SIM_good_count .+ SIM_bad_count)
  return
end





# function par_study_plot_err_projections_total(actual_par_study, count=300)
#   prms_names = actual_par_study.prms_names
#   prms_lists = actual_par_study.prms_lists
#   error_df = actual_par_study.SIMs_err
#   
#   range_to_display = 1:count
#   figure(4)
#   suptitle("The best $count of $(actual_par_study.name) w.r.t. all SIMs")
#   for i in 1:6
#     subplot(230+i)
#     title(prms_names[i])
#     plot([prms_lists[i][item[i]] for item in error_df.prms_indicies[range_to_display]], error_df.error_total[range_to_display], "x")
#   end
# end

### supporting stuff 

function get_show_tuple(actual_par_study; TC=Nothing, pO2=Nothing, bias=Nothing, simulations=Nothing)  
  SIMs_length = size(actual_par_study.SIM_list,1)
  show_array = Array{Int8}(undef,SIMs_length)
  show_array .= 1
  
  (TC == Nothing ? TC_filter = false : TC_filter = true)
  (pO2 == Nothing ? pO2_filter = false : pO2_filter = true)
  (bias == Nothing ? bias_filter = false : bias_filter = true)
  (simulations == Nothing ? simulations_filter = false : simulations_filter = true)
  for (i, SIM) in enumerate(actual_par_study.SIM_list)
    if TC_filter && !(SIM.TC in TC)
      show_array[i] = 0
      continue
    end
    if pO2_filter && !(SIM.pO2 in pO2) 
      show_array[i] = 0
      continue
    end
    if bias_filter && !(SIM.bias in bias)
      show_array[i] = 0
      continue
    end
    if simulations_filter && !(SIM.name in simulations)
      show_array[i] = 0
      continue
    end
  end
  return Tuple(show_array)
end

function get_SIMs_err_for_show_tuple(actual_par_study, show_tuple)
  SIMs_length = size(actual_par_study.SIM_list,1) 
  if size(show_tuple,1) != SIMs_length
    println("ERROR: get_SIMs_err_for_show_tuple: show_tuple has incorrect dimension!")
    throw(Exception)
  end
  err_sum = 0
  active_SIM_count = 0
  for i in 1:SIMs_length
    if show_tuple[i] == 1
      active_SIM_count += 1
    end
  end
  if active_SIM_count == 0
    println("ERROR: get_SIMs_err_for_show_tuple: show_tuple has no \"1\"!")
    throw(Exception)    
  end

  product_size_prms=1
  for list in actual_par_study.prms_lists
    product_size_prms*= size(list,1)
  end
  temp_error_df = DataFrame(
    prms_indicies = Array{Array{Int16}}(undef, product_size_prms),
    error_total = Array{Float32}(undef, product_size_prms)
    )


  temp_error_df.prms_indicies .= actual_par_study.SIMs_err.prms_indicies
  for i in 1:product_size_prms
    err_sum = 0
    for (j, SIM) in enumerate(actual_par_study.SIM_list)
      if show_tuple[j] == 1
        err_sum += actual_par_study.SIMs_err[!, j][i] * SIM.fitness_factor
      end
    end
    temp_error_df.error_total[i] = err_sum/active_SIM_count
  end
  sort!(temp_error_df, :error_total)
  return temp_error_df
end


### API functions ######################################
function par_study_plot_err_projections(actual_par_study, show_tuple; count=300)
#   if show_tuple == Nothing
#     show_tuple = get_show_tuple(actual_par_study)
#   end
  prms_names = actual_par_study.prms_names
  prms_lists = actual_par_study.prms_lists
  error_df = get_SIMs_err_for_show_tuple(actual_par_study, show_tuple)
  
  
  if count>size(error_df,1)
    println("Warning: par_study_plot_err_projections: count > number of data!")
    count = size(error_df,1)
  end
  
  range_to_display = 1:count
  figure(14)
  suptitle("The best $count of par study ($(actual_par_study.name)) w.r.t. $(show_tuple)")
  for i in 1:6
    subplot(230+i)
    ylabel("error")
    xlabel(prms_names[i])
    plot([prms_lists[i][item[i]] for item in error_df.prms_indicies[range_to_display]], error_df.error_total[range_to_display], "x")
  end
end

function par_study_plot_err_projections(actual_par_study; count=300,
                                            TC=Nothing, pO2=Nothing, bias=Nothing, simulations=Nothing)
  par_study_plot_err_projections(
    actual_par_study, get_show_tuple(actual_par_study, TC=TC, pO2=pO2, bias=bias, simulations=simulations),
    count=count
    )
end




function par_study_plot_the_best(actual_par_study, show_tuple; from=1, till=1, data_already_imported=false, plot_all_SIMs=false)
  par_study_dir = string(actual_par_study.save_dir, actual_par_study.name, "/")
  error_df = get_SIMs_err_for_show_tuple(actual_par_study, show_tuple)
  for (j, SIM) in enumerate(actual_par_study.SIM_list)
    if (show_tuple[j] == 1) || plot_all_SIMs
      typical_plot_exp(SIM, get_experiment(SIM))
      for i in from:till
        if data_already_imported
          SIM_sim = actual_par_study.SIMs_data[err_df.prms_indicies[i]...]
        else
          SIM_sim = load_file_prms(
              SIM, 
              save_dir=string(par_study_dir, string(SIM), "/"), 
              prms=get_prms_from_indicies(actual_par_study.prms_lists, error_df.prms_indicies[i]), 
              prms_names=actual_par_study.prms_names, 
              scripted_tuple=actual_par_study.scripted_tuple, 
              throw_exception=true, 
              verbose=false)
        end
        typical_plot_sim(SIM, SIM_sim, "$(i): $(get_prms_from_indicies(actual_par_study.prms_lists, error_df.prms_indicies[i]))")
      end
    end
  end
end


function par_study_plot_the_best(actual_par_study; from=1, till=1, plot_all_SIMs=false, data_already_imported=false,
                                  TC=Nothing, pO2=Nothing, bias=Nothing, simulations=Nothing)  
  par_study_plot_the_best(
    actual_par_study, get_show_tuple(actual_par_study, TC=TC, pO2=pO2, bias=bias, simulations=simulations),
    from=from, till=till, plot_all_SIMs=plot_all_SIMs, data_already_imported=data_already_imported
  )
end



function run_par_study(actual_par_study::par_study_struct; save_dir="../data/par_studies/", save_file_bool=false, mode="test", verbose=true)
     
     
  function consistency_check()
      
      
      # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # (ale az v pripade, ze to bude obecne napsane)
      
      return true
  end
  SIM_length = size(actual_par_study.SIM_list,1)
  
  SIM_err_counter = zeros(Int64,SIM_length)
  SIM_good_counter = zeros(Int64,SIM_length)
  all_counter::Int64 = 0
  
  #include("../src/models/$(actual_par_study.physical_model_name).jl")
  save_file_path = "$(save_dir)$(actual_par_study.name)/"

  #println(" >> number of sets of parameters: ",
  #    length(A0_list)*length(R0_list)*length(DGA_list)*length(DGR_list)*length(beta_list)*length(A_list))
  if verbose
    println("Start running par_study ...")
  end
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
    if verbose
      println(
        string(" <><><><><><><><><><><><> all_counter <><><> calculating ",all_counter," of ", overall_count)
      )
      println("parameters: ($(string_prms_names)) = ($(string_prms))")
    end
    
    for (i,SIM) in enumerate(actual_par_study.SIM_list)
        try

            SIM_sim = typical_run_simulation(SIM, prms_names_in, prms_values_in)
            
            if save_file_bool
              save_file_prms(SIM, SIM_sim, "$(save_file_path)$(string(SIM))/", prms_values_in, prms_names_in, actual_par_study.scripted_tuple)           
            end
            SIM_good_counter[i] += 1
            if verbose
              print("  $(string(SIM)) ok")
            end
          catch e
            if e isa InterruptException
                rethrow(e)
            else
                println(e)
                SIM_err_counter[i] += 1
                if verbose
                  print("<<<<< $(string(SIM))  FAIL! >>>>>")
                end
            end
        end        
    end
    println()    
  end
      
  for_each_prms_in_prms_lists(actual_par_study.prms_lists, perform_par_study)
  if verbose
    println(string("<<< run_par_study(): SIM_good_count / all_count  >>> ", SIM_good_counter,"  /  ", SIM_err_counter .+ SIM_good_counter ," >>>>>>>\n"))
  end
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
