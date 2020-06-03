using Printf
using PyPlot
#using Plots
using DataFrames
#using LeastSquaresOptim
using Optim
using CSV

using NNLS
#using LinearAlgebra
using LsqFit


# structure ############
# input:
# # # EEC_structure
# # # alpha_range
# # # list of initial guesses (maybe from a file???)
# # # data_set_string + experimental conditions
# 
# EEC_data struct 
#   function, 
#   list of parameters (strings), 
#   string representing EEC_structure
#   actual_fitted_parmeters
# build EEC (EEC_structure::string)
#   from string "L-RC-RCPE-RCPE" and return its representation (as a function) and list of parameters as strings
# import 
#   via implemented function
# perform fit!(alpha_range, list_of_initial_guesses)
#   use of LM .. .probably
# export_results
# 
#
#

L_units = 1.0e-6 # inductance is in micro-Henry


mutable struct EEC_data_struct
  structure::String
  #Z::Function # Z(omega, prms)::Complex
  prms_names::Array{String}
  prms_values::Array{Float32}
  #
  standard_deviations::Array{Float32}
  lower_limits_for_fitting::Array{Float32}
  upper_limits_for_fitting::Array{Float32}
  #
  EEC_data_struct() = new()
end

function EEC_data_struct(structure::String)
  this = EEC_data_struct()
  this.structure = structure
  return this
end


# function get_omega_and_prms_tuple_string(prms_names)
#   output = "(omega, "
#   for item in prms_names
#     output = "$(output)$item, "
#   end
#   output = output[1:end-1]
#   output = "$(output))"
#   return output
# end
  
function get_EEC(EEC_structure)

  
  EEC_actual = EEC_data_struct(EEC_structure)
  EEC_structure_splitted = split(EEC_actual.structure, "-")
  prms_names = Array{String}(undef, 0)
  
  impedance_string = "0"
  
  for (i, token) in enumerate(EEC_structure_splitted)
    if token=="R"
      append!(prms_names, ["R$(i)"])
    elseif token=="L"
      append!(prms_names, ["L$(i)"])
    elseif token=="RCPE"
      append!(prms_names, ["R$(i)", "C$(i)", "alpha$(i)"])
    elseif token=="RC"
      append!(prms_names, ["R$(i)", "C$(i)"])
    elseif token=="RL"
      append!(prms_names, ["R$(i)", "L$(i)"])
    else
      println("ERROR: build_EEC: Incorrect input of EEC_structure $(EEC_structure)")
      throw(Exception)
    end
  end
  
  EEC_actual.prms_names = prms_names
  return  EEC_actual
end

function g(EEC::EEC_data_struct, name::String)
  #@show name
  return EEC.prms_values[findall(x -> x==name, EEC.prms_names)[1]]
end

function evaluate_EEC(EEC::EEC_data_struct, omega)

  total_Z = 0
  for (i, token) in enumerate(split(EEC.structure, "-"))
    if token=="R"
      total_Z += g(EEC, "R$(i)")
    elseif token=="L"
      total_Z += im*omega*g(EEC, "L$(i)")*L_units
    elseif token=="RCPE"
      total_Z += g(EEC, "R$(i)")    /( 1 + (g(EEC, "R$(i)") * g(EEC, "C$(i)") * im*omega)^(g(EEC, "alpha$(i)")) )
    elseif token=="RC"
      total_Z += g(EEC, "R$(i)")    /( 1 + im * omega* g(EEC, "R$(i)") * g(EEC, "C$(i)") )
    elseif token=="RL"
      total_Z += im*omega* g(EEC, "L$(i)") * g(EEC, "R$(i)")       /( g(EEC, "R$(i)") + im * omega* g(EEC, "L$(i)") )
    else
      println("ERROR: evaluate_EEC: Incorrect input of EEC_structure $(EEC_structure)")
      throw(Exception)
    end
  end
  return total_Z
end

function get_EIS_from_EEC(EEC_actual::EEC_data_struct; f_range=[])
  EIS_out = DataFrame( f = [], Z = [])
  #@show EEC_actual.prms_names
  #@show EEC_actual.prms_values
  for f in f_range
    push!(
      EIS_out, 
      (f, evaluate_EEC(EEC_actual, 2*pi*f))
    )
  end
  return EIS_out
end



function HF_LF_check(prms_values)
  if prms_values[3]*prms_values[4] < prms_values[6]*prms_values[7]
    return true
  else
    return false
  end
end

function HF_LF_correction(EEC)
  if HF_LF_check(EEC.prms_values)
    aux = EEC.prms_values[6]
    EEC.prms_values[6] = EEC.prms_values[3]
    EEC.prms_values[3] = aux
    
    aux = EEC.prms_values[7]
    EEC.prms_values[7] = EEC.prms_values[4]
    EEC.prms_values[4] = aux
    
    aux = EEC.prms_values[8]
    EEC.prms_values[8] = EEC.prms_values[5]
    EEC.prms_values[5] = aux
  end
end

#
#     mask can be generated via string input, e.g. every "L" should be constant >>
#
function EEC_find_fit!(EEC_actual::EEC_data_struct, EIS_exp::DataFrame; mask=Nothing, alpha_low=0.0, alpha_upp=1.0, with_errors=false)
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
  
  function prepare_masked_stuff(mask, x0, lower_bounds, upper_bounds)
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
    return x0M, lowM, uppM
  end
  
  function get_EIS_value_from_gf(EIS_df, gf_list)
    output = []
    for gf in gf_list
      if mod(gf,2) == 0
        append!(output, real(EIS_df.Z[div(gf,2)]))
      else
        append!(output, imag(EIS_df.Z[div(gf,2)]))
      end
    end
    return output
  end

  x_previous = Nothing
  dramatic_tolerance = 1.0e-4
  function check_dramatic_change(x)
    for (i, item) in enumerate(x)
      if abs(x[i] - x_previous[i]) > dramatic_tolerance
        return true
      end
    end
    return false
  end
  
  function to_optimize(x)
    if !(check_x_in(x, lowM, uppM))
      #println("    OUT OF THE BOUNDS   \n")
      return 10000
    end
    
#     if !(HF_LF_check(prepare_prms(mask, x0, x)))
#       println("     HF_LF_check is violated !!! ")
#       return 10000
#     end
    
    EEC_actual.prms_values = prepare_prms(mask, x0, x)
    EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
    EIS_EEC_plot = get_EIS_from_EEC(EEC_actual, f_range=EIS_get_checknodes_geometrical((1, 10000, 10)...))
    SIM = EIS_simulation(800, 100, 0, use_DRT=false, plot_option="Bode Nyq", plot_legend=false)[1]
#     if check_dramatic_change(x)
#       typical_plot_sim(SIM, EIS_EEC_plot)
#     end
    err = fitnessFunction(SIM, EIS_EEC, EIS_exp)
#    println("~~~~~ LM e = $(err)\nx = $(x)")
    return err
  end
  
  function model(gf, x)    
    EEC_actual.prms_values = prepare_prms(mask, x0, x)
    EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
    SIM = EIS_simulation(800, 100, 0, use_DRT=false, plot_option="Bode Nyq", plot_legend=false)[1]
    #typical_plot_sim(SIM, EIS_EEC)
    #println("e = $(fitnessFunction(SIM, EIS_EEC, EIS_exp))\nx = $(x)")
    return get_EIS_value_from_gf(EIS_EEC, gf)
  end  
  
  prms_length = length(EEC_actual.prms_values)
  if mask == Nothing
    mask = Array{Int16}(undef, prms_length)
    mask .= 1
  end
  
  
  x0 = copy(EEC_actual.prms_values)
  
  lower_bounds = Array{Float32}(undef, prms_length)
  upper_bounds = Array{Float32}(undef, prms_length)
  lower_bounds_threshold = 0.000005
  for (i, name) in enumerate(EEC_actual.prms_names)
      if occursin("alpha", name)
        lower_bounds[i] = max(alpha_low, lower_bounds_threshold)
        upper_bounds[i] = alpha_upp
      else
        lower_bounds[i] = lower_bounds_threshold
        upper_bounds[i] = Inf
      end
  end
  
  for (i, name) in enumerate(EEC_actual.prms_names)
    lower_bounds[i] = max(lower_bounds[i], EEC_actual.lower_limits_for_fitting[i])
    upper_bounds[i] = min(upper_bounds[i], EEC_actual.upper_limits_for_fitting[i])
  end
  
  #######################################
  #######################################
  #######################################
  #######################################
  x0M, lowM, uppM = prepare_masked_stuff(mask, x0, lower_bounds, upper_bounds)
  #@show lowM
  #@show uppM
 
  x_previous = x0M
 
  begin # LeastSquaresOptim
#   fit_LSO = optimize(to_optimize, x0M, lower=lowM, upper=uppM, 
#     #Δ=1000000, 
#     x_tol=1.0e-14,
#     f_tol=1.0e-14, 
#     g_tol=1.0e-14, 
#     #autodiff=:central,
#     LevenbergMarquardt(),
#     #Dogleg(),
#     )
  end
 
  if with_errors
    begin # LsqFit
      gf_range = collect( 2 : 1 : 1 + 2*length(EIS_exp.f))
      y_data = get_EIS_value_from_gf(EIS_exp, gf_range)   
      fit = curve_fit(model, gf_range, y_data, x0M, lower=lowM, upper=uppM)
      EEC_actual.prms_values = prepare_prms(mask, x0, fit.param)
    end  
  else
    begin # Optim
      fit_O = optimize(to_optimize, x0M, #lower=lowM, upper=uppM, 
        #Δ=1000000, 
        #x_tol=1.0e-14,
        #f_tol=1.0e-14, 
        #g_tol=1.0e-14, 
        #autodiff=:central,
        #LevenbergMarquardt(),
        #Dogleg(),
        )
      EEC_actual.prms_values = prepare_prms(mask, x0, fit_O.minimizer)  
    end
  end

  
  error0 = deepcopy(x0)
  error0 .= -Inf;
  try
    EEC_actual.standard_deviations = prepare_prms(mask, error0, stderror(fit))
  catch
    #println("SO EIN PECH!")
    EEC_actual.standard_deviations = deepcopy(EEC_actual.prms_values)
  end
  return
end













function get_left_right_width_of_EIS(EIS_df, N_for_sum=2)
  real_sum_left = 0
  real_sum_right = 0  
  for i in 1:N_for_sum
    real_sum_left += real(EIS_df.Z[end-i+1])
    real_sum_right += real(EIS_df.Z[i])
  end
  left = real_sum_left/N_for_sum
  right = real_sum_right/N_for_sum
  width = right - left
  return left, right, width
end

function EIS_data_preprocessing(EIS_df)
  
  #lowest frequency cut off
  lowest_freq_idx = -1
  positive_counter = 0
  for (i, Z) in enumerate(EIS_df.Z)
    if imag(Z) < 0
      positive_counter += 1
    else
      positive_counter = 0
    end
    if positive_counter == 10
      lowest_freq_idx = i - 9
      break
    end
  end
  
  # intersection with x axis
  x_intersection_freq_idx = -1
  negative_counter = 0
  for i in (lowest_freq_idx + 5):length(EIS_df.f)
    if imag(EIS_df.Z[i]) > 0
      negative_counter += 1
    else  
      negative_counter = 0
    end
    if negative_counter == 8
      x_intersection_freq_idx = i - 7
      break
    end
  end
  
  #@show x_intersection_freq_idx
  #@show lowest_freq_idx
  
  if x_intersection_freq_idx == -1
    return DataFrame(f = EIS_df.f[lowest_freq_idx:end], Z = EIS_df.Z[lowest_freq_idx:end])
  else
    # inductance cut off
    accepted_inductance_real_axis_threshold = 0.05*real(EIS_df.Z[lowest_freq_idx]) + 0.95*real(EIS_df.Z[x_intersection_freq_idx])
    highest_freq_idx = -1
    for i in (x_intersection_freq_idx + 1 ):length(EIS_df.f)
      if real(EIS_df.Z[i]) > accepted_inductance_real_axis_threshold
        highest_freq_idx = i
        break
      end
    end
    if highest_freq_idx == -1
      return DataFrame(f = EIS_df.f[lowest_freq_idx:end], Z = EIS_df.Z[lowest_freq_idx:end])
    end

    return DataFrame(f = EIS_df.f[lowest_freq_idx:highest_freq_idx], Z = EIS_df.Z[lowest_freq_idx:highest_freq_idx])
  end
end

function get_init_values(EIS_exp)
  
  # R | L | RCPE | RCPE structure REQUIERED!
  output = Array{Float64}(undef, 8)
  
  # resistors
  smaller_circle_ratio = 0.3  # with respect to full width
  left, right, width = get_left_right_width_of_EIS(EIS_exp)
  
  
  # R_ohm
  output[1] = left
  
  # L
  output[2] = 1.54
#   if imag(EIS_exp.Z[end]) > 0
#     @show  imag(EIS_exp.Z[end])
#     output[2] = imag(EIS_exp.Z[end])/(2*pi*EIS_exp.f[end])*(1/L_units)
#   else
#     output[2] = 1.0
#   end
  
  # R3, R4
  output[3] = width*smaller_circle_ratio
  output[6] = width*(1 - smaller_circle_ratio)
  
  # alphas
  output[5] = 0.8
  output[8] = 0.9
  
  # C3, C4
  output[4] = 1.0
  output[7] = 0.5
  
  #@show imag(EIS_exp.Z[end]) 
  #@show left, right, width
  #@show output
  return output
  #return [2.0, 6.2, 0.5 , 0.001, 1.0,    0.6, 0.01, 0.8]
end

function set_fitting_limits_to_EEC_from_EIS_exp!(EEC::EEC_data_struct, EIS_exp)

  left, right, width = get_left_right_width_of_EIS(EIS_exp)
    
  prms_length = length(EEC.prms_names)
  lower_limits = Array{Float32}(undef, prms_length)
  upper_limits = Array{Float32}(undef, prms_length)
  
  # R | L | RCPE | RCPE structure REQUIERED!
  
  # R_ohm
  lower_limits[1] = left/2
  upper_limits[1] = right
  
  # L2
  lower_limits[2] = (imag(EIS_exp.Z[end])/(2*pi*EIS_exp.f[end])*(1/L_units))/10
  upper_limits[2] = (imag(EIS_exp.Z[end])/(2*pi*EIS_exp.f[end])*(1/L_units))*10

  # R3
  lower_limits[3] = -Inf
  upper_limits[3] = width*2
  
  # R4
  lower_limits[6] = -Inf
  upper_limits[6] = width*2
  
  # C3
  lower_limits[4] = -Inf
  upper_limits[4] = 10
  
  # C4
  lower_limits[7] = -Inf
  upper_limits[7] = 10
  
  # alpha3
  lower_limits[5] = -Inf
  upper_limits[5] = Inf
  
  # alpha4
  lower_limits[8] = -Inf
  upper_limits[8] = Inf
  
  EEC.lower_limits_for_fitting = lower_limits
  EEC.upper_limits_for_fitting = upper_limits
  return
end


function view_EEC(;
      EEC_structure="RL-R-L-RCPE-RCPE", 
      prms_values=[1, 0.001,      1.7, 0,        1. , 0.001, 1.0,    1.0, 0.01, 0.8], 
      f_range=(0.01, 10000, 1.2),
      print_bool=false, pyplot=1, plot_legend=true, use_DRT=true
      )
  
  EEC = get_EEC(EEC_structure)
  
  function recursive_view_EEC_call(output_prms, plot_names, plot_values, active_idx)
      if active_idx > size(prms_names,1)

        ###################################################################
        # for each combination of output_prms from prms_lists perform this #######
        prms_names_in=[]
        prms_values_in=[]
        
        append!(prms_names_in, prms_names)
        append!(prms_values_in, output_prms)
      
        
        if check_equal_size(prms_names, prms_values)  
          EEC.prms_values = output_prms
          EIS_EEC = get_EIS_from_EEC(EEC, f_range=EIS_get_checknodes_geometrical(f_range...))   
        end    
        if size(plot_names,1) < 1
          plot_prms_string = ""
        else
          plot_prms_string = " $(string(plot_names)) = $(plot_values)"
        end
        if print_bool
          #@show EEC.prms_names
          #@show EEC.prms_values
          
          println(get_string_EEC_prms(EEC, plot_errors=false))
        end
        if pyplot > 0
            ysz_fitting.typical_plot_exp(EIS_simulation(800, 100, 0.0, plot_legend=plot_legend, use_DRT=use_DRT)[1], EIS_EEC, "!EEC"*plot_prms_string) 
        end 
        
        
        return
        ###################################################################
        ###################################################################
      
      end
      if size(prms_values[active_idx],1)>1
        for i in prms_values[active_idx]
          recursive_view_EEC_call(
            push!(deepcopy(output_prms),i),
            push!(deepcopy(plot_names),prms_names[active_idx]),
            append!(deepcopy(plot_values),i),
            active_idx + 1)
        end
      else
        recursive_view_EEC_call(
          push!(deepcopy(output_prms),prms_values[active_idx][1]),
          plot_names,
          plot_values,
          active_idx + 1)
      end
    end
  
  prms_names = EEC.prms_names
  recursive_view_EEC_call([], Array{String}(undef,(0)), Array{Float64}(undef,(0)), 1)
  
  return EEC
end

function get_string_EEC_prms(EEC::EEC_data_struct; plot_errors=false)
  output = ""
  for (i, item) in enumerate(EEC.prms_values)
    output *= "$(@sprintf("%8s", EEC.prms_names[i]))  =  $(@sprintf("%12.8f", item))    "
    plot_errors ? (output *="(rel_err: $(@sprintf("%12.8f", EEC.standard_deviations[i]/item)))\n") : output *= "\n"
  end
  return output
end









mutable struct EEC_data_holder_struct
  TC
  pO2
  bias
  data_set
  
  prms_names
  # indexing of data: [TC, pO2, bias, data_set, prms_values]
  data::Array{Any}
  
  EEC_data_holder_struct() = new()
end

function EEC_data_holder_struct(TC, pO2, bias, data_set, prms_names)
  this = EEC_data_holder_struct()
  this.TC = TC
  this.pO2 = pO2
  this.bias = bias
  this.data_set = data_set
  this.prms_names = prms_names
  
  this.data = Array{Any}(undef, (length(TC), length(pO2), length(bias), length(data_set), length(prms_names)))
  return this
end


function make_array_from_string(input)
  if typeof(input) == String
    return [input]
  else
    return input
  end
end

        
function save_EEC_prms_to_file(TC, pO2, bias, data_set, prms_values, saving_destination)
  df_out = DataFrame(TC = TC, pO2 = pO2, bias = bias, data_set = data_set,
                    R1      = [prms_values[1]],
                    L2      = [prms_values[2]],
                    R3      = [prms_values[3]],
                    C3      = [prms_values[4]],
                    alpha3  = [prms_values[5]],
                    R4      = [prms_values[6]],
                    C4      = [prms_values[7]],
                    alpha4  = [prms_values[8]]
                    )
  CSV.write(saving_destination, df_out, delim="\t", append=true)
end

function run_EEC_fitting(;TC=800, pO2=80, bias=0.0, data_set="MONO_110",
                        f_interval=Nothing, succes_fit_threshold = 0.002,
                        fixed_prms_names=[], fixed_prms_values=[],
                        init_values=Nothing, alpha_low=0, alpha_upp=1, #fitting_mask=Nothing,
                        plot_bool=false, plot_legend=true, plot_initial_guess=false, plot_fit=true,
                        with_errors=false,
                        use_DRT=false,
                        save_file_bool=true, save_to_folder="../data/EEC/", file_name="default_output.txt")
  ####
  ####  TODO:
  ####  [?] initial values handling (from file, probably)
  ####  [x] R_ohm estimation (+ R2, R3 estimation?)
  ####  [x] swap RCPEs in the case Hf and Lf interchanged
  ####  [x] check, if fitting converged
  ####  [x] output format
  ####  [x] view_EEC ... recursive call
  ####  [x] define default SIM for ploting Nyquists! and turn of the legend!
  ####  [x] f_interval ... to crop the Nyquist
  ####  [ ] sjednotit plot options, abych mohl zapinat DRT a menit ploty v klidecku!
  ####  [x] upravit hranice parametrickeho prostoru, aby nehledal kraviny!
  ####  [ ] fixed_prms_names ... bug for *3 a *4  ... elements can be interchanged !!!
  ####  [ ] maybe try several automatic initial guesses and pick up the best one?
  
  
  ########## Questions
  ####  [x] Jak se mam zachovat, kdyz fit neni dobry? Mam proste preskocit, nic nezapisovat do souboru a jet dal? Nebo zapsat a warning?
  ####  [x] Format zapisu do souboru -> nemel by to spis byt *.csv soubor s cisly oddelenymi carkami ci tabulatory?
  ####  [x] da se orezat Nyquist pro nizke frekvence, at nejde do zapornych cisel. Chceme?
  ####  [x] teoreticky by se plotovani obrazku mohlo vzdy vypnout, kdyz by clovek chtel ukladat do souboru. Ale myslim, ze to neni nutne
  ####  [x] nazor na to, ze pocitam prumer systemu jako 1.2 cm, pricemz to je jne prumer elektrody, ale ellyt ma 2.5 cm v prumeru
  ####  [x] da se udelat nejake moudre orezani, ktere nedovoli vice hodnotam jit doprava oproti pruseciku s realnou osou
  ####  [ ] nechcete udelat treba animace nebo dalsi obrazky? (vyhledove)
  ####  [ ] chcete vypisovat L2 v mikro-Henry nebo v Henry? stejne jako zadavani a kdekoliv
  
  function succesful_fit(EIS_EEC, EIS_exp)
    if fitnessFunction(EIS_simulation(), EIS_EEC, EIS_exp) > succes_fit_threshold
      return false
    else
      return true
    end
  end

  if save_file_bool
    mkpath(save_to_folder)
  end
  
  EEC_actual = get_EEC("R-L-RCPE-RCPE")
      
  data_set = make_array_from_string(data_set)  
  EEC_data_holder = EEC_data_holder_struct(TC, pO2, bias, data_set, EEC_actual.prms_names)
  
  warning_buffer = ""
  if save_file_bool 
    saving_destination = save_to_folder*file_name
    df_out = DataFrame(TC = [], pO2 = [], bias = [], data_set = [],
                    R1      = [],
                    L2      = [],
                    R3      = [],
                    C3      = [],
                    alpha3  = [],
                    R4      = [],
                    C4      = [],
                    alpha4  = []
                    )
    CSV.write(saving_destination, df_out, delim="\t")
  end
  
  cycle_number = 0
  for   (TC_idx, TC_item) in enumerate(TC), 
        (pO2_idx, pO2_item) in enumerate(pO2), 
        (bias_idx, bias_item) in enumerate(bias), 
        (data_set_idx, data_set_item) in enumerate(data_set)
    cycle_number += 1

  
    SIM_list = ysz_fitting.EIS_simulation(TC_item, pO2_item, bias_item, 
                  use_DRT=use_DRT, plot_option="Bode Nyq DRT RC", plot_legend=plot_legend, data_set=data_set_item)
    SIM = SIM_list[1]    
    
    
    
    
    
    try
      EIS_exp = ysz_fitting.import_data_to_DataFrame(SIM)
    catch 
      warning_buffer *= " ->->-> ERROR: file for (TC, pO2, bias, data_set_item) = ($(TC_item), $(pO2_item), $bias_item, $data_set_item) ... NOT FOUND! \n"
      missing_array = Array{Any}(undef, 8)
      missing_array .= Missing
      
      EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, :] = deepcopy(missing_array)
      if save_file_bool
        save_EEC_prms_to_file(TC_item, pO2_item, bias_item, data_set_item, missing_array, save_to_folder*file_name)
      end
      continue
    end

    
    
    
    
    
    
    
    if f_interval!=Nothing
      if f_interval == "auto"
        #typical_plot_exp(SIM, EIS_exp, "! before")
        EIS_exp = EIS_data_preprocessing(EIS_exp)
        #typical_plot_exp(SIM, EIS_exp, "! after")
      else
        EIS_exp = EIS_crop_to_f_interval(EIS_exp, f_interval)
      end
    end
    
    function get_fitting_mask_and_apply_fixed_prms_values(EEC, init_values, fixed_prms_names, fixed_prms_values)
      if length(fixed_prms_names) != length(fixed_prms_values)
        println("ERROR: length(fixed_prms_names) != length(fixed_prms_values) ==>> $(length(fixed_prms_names)) != $(length(fixed_prms_values))")
        return throw(Exception)
      else
        output = []
        skip_bool = false
        for (i, original_name) in enumerate(EEC.prms_names)
          skip_bool = false
          for (j, name) in enumerate(fixed_prms_names)
            if name == original_name
              append!(output, 0)
              init_values[i] = fixed_prms_values[j]
              skip_bool = true
              break
            end
          end
          skip_bool || append!(output, 1)
        end
        if length(fixed_prms_names) + sum(output) != length(EEC.prms_names)
          println("ERROR: some of fixed_prms_names $(fixed_prms_names) not found!")
          return throw(Exception)
        end
        #@show output
        return output
      end
    end
    
#     if fitting_mask == Nothing
#       mask = Array{Int16}(undef, length(EEC_actual.prms_names))
#       mask .= 1
#     else
#       mask = fitting_mask
#     end
    set_fitting_limits_to_EEC_from_EIS_exp!(EEC_actual, EIS_exp)

    
    # test!
    #   EEC_actual.prms_values = [0.00001, 2.5, 1, 0.001, 2, 0.005, 0.4 ]
    #   EIS_exp = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
    #
    
    if init_values == Nothing || cycle_number > 1
      init_values = get_init_values(EIS_exp)
    end
    
    mask = get_fitting_mask_and_apply_fixed_prms_values(EEC_actual, init_values, fixed_prms_names, fixed_prms_values)
    
    begin
      begin
        plot_bool && ysz_fitting.typical_plot_exp(SIM, EIS_exp)
        
        EEC_actual.prms_values = init_values        
        
        if !save_file_bool
          println("prms_names = $(EEC_actual.prms_names)")
          println("init_values = $(EEC_actual.prms_values)")
        end
        
        EIS_EEC_pre = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
        if plot_initial_guess
          plot_bool && ysz_fitting.typical_plot_sim(SIM, EIS_EEC_pre, "!EEC initial guess")
        end
        
        EEC_find_fit!(EEC_actual, EIS_exp, mask=mask, alpha_low=alpha_low, alpha_upp=alpha_upp, with_errors=with_errors)
        
        if !save_file_bool
          #println("FITTED_values = $(EEC_actual.prms_values)   <<<<<<<<<<<<<<<<<<<<<< ")
        end        
        
        EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
        HF_LF_correction(EEC_actual)
        if !(succesful_fit(EIS_EEC, EIS_exp))
          warning_buffer *="WARNING: (TC, pO2, bias, data_set_item) = ($(TC_item), $(pO2_item), $bias_item, $data_set_item) maybe NOT CONVERGED !!! // error $(fitnessFunction(EIS_simulation(), EIS_EEC, EIS_exp)) > defined threshold $(succes_fit_threshold) //\n"
          # TODO ... optionally save picture of the result for humanous check
          # or this reports can go to another error_log_file.txt
        end

        
        output_string = "========== TC, pO2, bias: ($TC_item, $(pO2_item), $bias_item) --- data_set_item = $data_set_item) =========="
        output_string *= "\n$(get_string_EEC_prms(EEC_actual))"
        if save_file_bool
          save_EEC_prms_to_file(TC_item, pO2_item, bias_item, data_set_item, EEC_actual.prms_values, save_to_folder*file_name)
        else
          println(output_string)
          println(">>> error >>> ", fitnessFunction(EIS_simulation(), EIS_EEC, EIS_exp))
        end
        
        EIS_EEC = get_EIS_from_EEC(EEC_actual, f_range=EIS_exp.f)
        if plot_fit
          plot_bool && ysz_fitting.typical_plot_sim(SIM, EIS_EEC, "!EEC fit")
        end
        
        EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, :] = deepcopy(EEC_actual.prms_values)
        
      end
    end
  end
  
  println()
  println("          =================================")
  println("          =================================")
  println("          =========== MESSAGES ============")
  if length(warning_buffer) > 0 
    println(warning_buffer)
  else
    println(" no warnings :) ")
  end
  
  return EEC_data_holder
#   T_holder = TC_holder .+ 273.15
#   
#   if length(bias_holder) > 1 || length(T_holder) > 1
#     if bias_holder[1] == bias_holder[2]
#       figure(42)
#       title("sigma VS bias")
#       plot(1000 ./ (T_holder), log.(sigma_holder .* T_holder), label="bias "*string(bias_holder[1]))
#       legend(loc="best")
#       xlabel("1000/T")
#       ylabel("log T sigma")
#     else
#       figure(44)
#       title("sigma VS \"TC\"")
#       plot(bias_holder, sigma_holder, label="TC "*string(TC_holder[1]))
#       legend(loc="best")
#       xlabel("bias")
#       ylabel("sigma")
#     end
#   end
#   
  
end

function save(EEC_data_holder, folder, name)
  # TODO!
end

function load(EEC_data_holder)
  # TODO!
end


function plot_EEC_data_general(EEC_data_holder; 
                                y_name="C3",
                                x_name="bias",
                                #
                                TC_interval = [-Inf, Inf],
                                pO2_interval = [-Inf, Inf],
                                bias_interval = [-Inf, Inf],
                                data_set = Nothing,                              
                                #
                                fig_num=102, 
                                plot_legend=true)
  
  function make_range_list(prm_name, interval)
    output_list = []
    if prm_name == "data_set"
      for (i, item) in enumerate(getfield(EEC_data_holder, Symbol(prm_name)))
        if item in interval
          append!(output_list, i)
        end
      end
    else
      for (i, item) in enumerate(getfield(EEC_data_holder, Symbol(prm_name)))
        if interval[1] <= item && item <= interval[2]
          append!(output_list, i)
        end
      end
    end
    return output_list
  end
  
  function add_legend_contribution(idx, plot_idx)
    if length(plot_idx) == 1 && length(range_list[idx]) > 1
      return "$(identifier_list[idx])=$(getfield(EEC_data_holder, Symbol(identifier_list[idx]))[plot_idx]) "
    else
      return ""
    end
  end
  
  x_name_idx = -1
  identifier_list = ["TC", "pO2", "bias", "data_set"]
  range_list = Array{Any}(undef, 4)
  for i in 1:4
    if x_name == identifier_list[i]
      x_name_idx = i
    end
  end
  y_name_idx = findall(x -> x==y_name, EEC_data_holder.prms_names)[1]
  
  
  range_list[1] = make_range_list("TC", TC_interval)
  range_list[2] = make_range_list("pO2", pO2_interval)
  range_list[3] = make_range_list("bias", bias_interval)
  if data_set != Nothing 
    data_set = make_array_from_string(data_set)
  else
    data_set = make_array_from_string(EEC_data_holder.data_set)
  end
  range_list[4] = make_range_list("data_set", data_set)

  figure(fig_num)
  xlabel(x_name)
  ylabel(y_name)
  title(y_name*"  vs  "*x_name)
  
  for TC_idx in (x_name_idx != 1 ? range_list[1] : [range_list[1]]),
      pO2_idx in (x_name_idx != 2 ? range_list[2] : [range_list[2]]),
      bias_idx in (x_name_idx != 3 ? range_list[3] : [range_list[3]]),
      data_set_idx in (x_name_idx != 4 ? range_list[4] : [range_list[4]])
    
    to_plot = EEC_data_holder.data[TC_idx, pO2_idx, bias_idx, data_set_idx, y_name_idx]
    
    #@show to_plot
    
    plot(
      getfield(
        EEC_data_holder, 
        Symbol(identifier_list[x_name_idx])
      )[range_list[x_name_idx]], 
      to_plot,
      label= "$(add_legend_contribution(1, TC_idx))$(add_legend_contribution(2, pO2_idx))$(add_legend_contribution(3, bias_idx))$(add_legend_contribution(4, data_set_idx))"
    )
    
    if plot_legend
      legend(loc="best")
    end
  end
  
  PyPlot.show()
  return
end
