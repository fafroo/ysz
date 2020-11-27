using CSV
using DataFrames
using PyPlot



const EIS_standard_figure_num = 6
const EIS_default_physical_model_name = "ysz_model_GAS_LoMA_Temperature"

include("../DRT.jl")


######################################


mutable struct EIS_simulation <: abstract_simulation
  TC::Float64
  pO2::Float64
  bias::Float64
  data_set::String
  #
  physical_model_name::String
  #
  dx_exp::Float64
  f_range::Tuple
  f_interval::Tuple
  fig_size::Tuple
  fig_num::Int16
  #
  checknodes::Any
  fitness_factor::Float64
  #
  use_TDS::Int16
  tref::Float64
  #
  use_DRT::Bool
  DRT_control::DRT_control_struct
  DRT_draw_semicircles::Bool
  #
  plot_option::String
  plot_legend::Bool
  #
  name::String
  ID::Int16
  
  EIS_simulation() = new()
end

function string(SIM::EIS_simulation)
  return "$(SIM.name)_TC_$(SIM.TC)_pO2_$(SIM.pO2)_bias_$(SIM.bias)"
end

function EIS_simulation(TC, pO2, bias=0.0; data_set="MONO_110", physical_model_name=EIS_default_physical_model_name, dx_exp=-9, f_range=EIS_get_shared_f_range(), f_interval=(-Inf, Inf), fig_size=(9, 6), fig_num=-1, fitness_factor=1.0, use_TDS=0, tref=0, use_DRT=true, DRT_control=DRT_control_struct(), DRT_draw_semicircles=false, plot_option="Nyq Bode DRT RC leg", plot_legend=false)
  output = Array{abstract_simulation}(undef,0)
  for TC_item in TC
    for pO2_item in pO2
      for bias_item in bias
        this = EIS_simulation()
        #
        this.TC = TC_item
        this.pO2 = pO2_item
        this.bias = bias_item
        this.data_set = data_set
        #
        this.physical_model_name=physical_model_name
        #
        this.dx_exp = dx_exp
        this.f_range = f_range
        this.f_interval = f_interval
        this.fig_size = fig_size
        this.fig_num = (fig_num >= 0 ? fig_num : EIS_standard_figure_num)
        #
        this.checknodes = EIS_get_checknodes_geometrical(f_range...)
        this.fitness_factor = fitness_factor
        #
        this.use_TDS = use_TDS
        this.tref = tref
        #
        this.use_DRT = use_DRT
        this.DRT_control = DRT_control
        this.DRT_draw_semicircles = DRT_draw_semicircles
        #
        this.plot_option = plot_option
        this.plot_legend = plot_legend
        #
        this.name = "EIS"
        this.ID = 2
        
        push!(output, this)
      end
    end
  end
  return output
end

#######################################
# potentially changable stuff #########

function EIS_get_shared_f_range()
    # experimental range is 0.1 - 65000 Hz
    # fs = (w0, w1, w_fac)
    return f_range = (1.1, 9500, 1.2)
end
  
function get_shared_checknodes(SIM::EIS_simulation)
    return EIS_get_checknodes_geometrical(EIS_get_shared_f_range()...)
end

function setting_legend(SIM::EIS_simulation; latex=true)
  if latex
    return "\$\\theta=$(SIM.TC)\$°C \$\\mathrm{O}_2=$(SIM.pO2)\\% \$ \$\\mathrm{bias}=$(SIM.bias)\$ $(SIM.data_set)"
  else
    return "TC=$(SIM.TC)°C pO2=$(SIM.pO2)% bias=$(SIM.bias)V $(SIM.data_set)"
  end
end

function typical_plot_general(SIM::EIS_simulation, EIS_df, my_label, additional_string="", to_standard_figure=true; marker_style="x-")
  if to_standard_figure
    figure(SIM.fig_num, figsize=SIM.fig_size)
  end

  if  (occursin("Rf", SIM.plot_option)) || 
      (occursin("RC", SIM.plot_option)) ||
      (occursin("DRT", SIM.plot_option)) ||
      (occursin("Rtau", SIM.plot_option))
    num_rows = 300
  else
    num_rows = 200
  end
  
  if occursin("leg", SIM.plot_option) || SIM.plot_legend
    sLEG = subplot(num_rows + 23)
    plot([], marker_style, label = my_label)
    PyPlot.axis("off") 
    legend(loc="best")
  end
  
  if occursin("Nyq", SIM.plot_option)
    s1 = subplot(num_rows + 21)
    title("Nyquist plot")
    xlabel("Re\$(Z) \\ [\\Omega]\$")
    ylabel("-Im\$(Z) \\ [\\Omega]\$")
    plot(real(EIS_df.Z), -imag(EIS_df.Z), marker_style, label = my_label)
#     if !(my_label == "") && SIM.plot_legend
#         legend(loc="best")
#     end
    grid(true)
    s1.set_aspect(1.0)
  end
  


  if occursin("Bode", SIM.plot_option)
    s2 = subplot(num_rows + 22)
    title("Bode plot: Re vs log f \$(\$Hz\$)\$")
    ylabel("Re\$(Z) \\ [\\Omega]\$")
    plot(log10.(EIS_df.f), real(EIS_df.Z), marker_style, label = my_label)
      
    s3 = subplot(num_rows + 24)
    title("Bode plot: - Im vs log f \$(\$Hz\$)\$")
    ylabel("-Im\$(Z) \\ [\\Omega]\$")
    plot(log10.(EIS_df.f), -imag(EIS_df.Z), marker_style, label = my_label)
  end
  
  if SIM.use_DRT
    DRT_actual = get_DRT(EIS_df, SIM.DRT_control)
    
    if occursin("DRT", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot_DRT_h(DRT_actual, false)
    elseif occursin("Rtau", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot_DRT_Rtau(DRT_actual, false)
    elseif occursin("Rf", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot_DRT_Rf(DRT_actual, false)
    end
    
    if occursin("RC", SIM.plot_option) 
      s4 = subplot(num_rows + 26)
      plot_DRT_RC(DRT_actual, false)
    end
    
    if SIM.DRT_draw_semicircles
      R_ohm_auxilary = DRT_actual.R_ohm
      for i in 1:size(DRT_actual.peaks_df, 1)
        EIS_RC = EIS_get_RC_CPE_elements(DRT_actual.peaks_df.R[i], DRT_actual.peaks_df.C[i], 0, 1, 1, R_ohm_auxilary, f_range=SIM.f_range)
        SIM_aux = deepcopy(SIM)
        SIM_aux.use_DRT = false
        typical_plot_general(SIM_aux, EIS_RC, "", "", true, marker_style=marker_style)
        R_ohm_auxilary += DRT_actual.peaks_df.R[i]
      end
    end
  else    
    if occursin("DRT", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot([])
    elseif occursin("Rtau", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot([])
    elseif occursin("Rf", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot([])
    end
    
    if occursin("RC", SIM.plot_option)
      s4 = subplot(num_rows + 26)
      plot([])
    end
  end
  

  plt.subplots_adjust(bottom=0.07, top=0.95)  
end

function typical_plot_sim(SIM::EIS_simulation, EIS_df, additional_string="", to_standard_figure=true)
  if length(additional_string)>0 && additional_string[1]=='!'
    my_label = additional_string[2:end]
  else
    my_label = "sim $(setting_legend(SIM))$(additional_string)"
  end
  typical_plot_general(SIM, EIS_df, my_label, additional_string, to_standard_figure, marker_style="x-")
end

function typical_plot_exp(SIM::EIS_simulation, EIS_df, additional_string="", to_standard_figure=true)
  if length(additional_string)>0 && additional_string[1]=='!'
    my_label = additional_string[2:end]
  else
    my_label = "exp $(setting_legend(SIM))$(additional_string)"
  end
  typical_plot_general(SIM, EIS_df, my_label, additional_string, to_standard_figure, marker_style="x:")
end

function fitness_error_report(SIM::EIS_simulation, plot_prms_string, EIS_exp, EIS_sim)
  println("EIS fitness error $(setting_legend(SIM, latex=false)) $(plot_prms_string)  => ", fitnessFunction(SIM, EIS_sim, EIS_exp))
end


#######################################
# relatively static part ##############

function apply_checknodes(SIM::EIS_simulation, EIS_in, checknodes)
    DataFrame( f = checknodes, Z = EIS_get_Z_values(EIS_in, checknodes))
end

function load_file_prms(SIM::EIS_simulation; save_dir, prms, prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), scripted_tuple=(0, 0, 0, 0, 0, 0), throw_exception=true, verbose=true)

    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="EIS", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    
    df_out = DataFrame()
    try
      df = CSV.read(string(out_path, out_name))
      df_out = DataFrame(f = df.f, Z = (df.ReZ .+ df.ImZ .*im))
    catch e
      if e isa InterruptException
        rethrow(e)
      else
        if verbose
          println("file not found: $(out_path)$(out_name)")
        end
        if throw_exception
          rethrow(e)
        end
      end
    end
    return df_out
end

function save_file_prms(SIM::EIS_simulation, df_out, save_dir, prms, prms_names, scripted_tuple; mode="")
    if mode=="exp"
      prefix = "$(string(SIM))_experiment"
      scripted_tuple = Array{Int16}(undef, length(prms_names))
      scripted_tuple .= 0
      prms = []
      prms_names = []
    elseif mode=="sim"
      prefix = "$(string(SIM))_"
      scripted_tuple = Array{Int16}(undef, length(prms_names))
      scripted_tuple .= 0
    else
      prefix="EIS"
    end
    
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix=prefix, prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    mkpath(out_path)

    CSV.write(
        string(out_path,out_name),
        DataFrame(f = df_out.f, ReZ = real(df_out.Z), ImZ = imag(df_out.Z))
    )
    return
end


function typical_run_simulation(SIM::EIS_simulation, prms_names_in, prms_values_in, pyplot::Int=0)
  EIS_df = ysz_experiments.run_new(
      pyplot=(pyplot == 2 ? true : false), 
      EIS_TDS=(SIM.use_TDS > 0 ? true : false), tref=SIM.tref,
      EIS_IS=(SIM.use_TDS > 0 ? false : true), out_df_bool=true, bias=SIM.bias, f_range=SIM.f_range,
      dx_exp=SIM.dx_exp,      
      T=TCtoT(SIM.TC), pO2=pO2tosim(SIM.pO2), data_set=SIM.data_set,
      prms_names_in=prms_names_in,
      prms_values_in=prms_values_in,
      physical_model_name=SIM.physical_model_name
  )
end

function import_data_to_DataFrame(SIM::EIS_simulation)
  import_EIStoDataFrame(TC=SIM.TC, pO2=SIM.pO2, bias=SIM.bias, data_set=SIM.data_set)
end

function EIS_test_checknodes_range(f_range=EIS_get_shared_f_range())
  Nyquist_plot(apply_checknodes(EIS_simulation(), import_EIStoDataFrame(TC=800,pO2=100,bias=0.0),EIS_get_checknodes_geometrical(f_range...)))
end


function EIS_view_experimental_data(;TC, pO2, bias, data_set="MONO_110", plot_option="Nyq Bode Rtau RC", f_interval=Nothing, use_checknodes=false, plot_legend=true, fig_num=12, L=0.0)  
  
  EIS_exp = DataFrame()
  for TC_item in TC
    for pO2_item in pO2
      for bias_item in bias, data_set_item in (typeof(data_set)==String ? [data_set] : data_set)
        if use_checknodes
          SIM = EIS_simulation()
          checknodes =  get_shared_checknodes(SIM)
          EIS_exp = apply_checknodes(SIM, import_EIStoDataFrame(TC=TC_item, pO2=pO2_item, bias=bias_item, data_set=data_set_item), checknodes)
        else
          EIS_exp = import_EIStoDataFrame(TC=TC_item, pO2=pO2_item, bias=bias_item, data_set=data_set_item)
        end
        figure(fig_num)
        EIS_exp = f_interval_preprocessing(EIS_exp, f_interval)        
        EIS_exp = EIS_add_inductance(EIS_exp, L)
        typical_plot_exp(EIS_simulation(TC_item, pO2_item, bias_item, data_set=data_set_item, plot_option=plot_option, plot_legend=plot_legend)..., EIS_exp, "", false)
      end
    end
  end
  return EIS_exp
end

# function EIS_capacitance_test()
#   
# end
  
function EIS_get_Z_values(Q_raw, checknodes)
    Zx_values = []
    #println(Q_raw[!, 1])
    #println(Q_raw[!, 2])

    start_i = 1
    for row = 1:size(checknodes,1)
        Ux = checknodes[row]
        xk = -1
        i = start_i
        while i < length(Q_raw[!, 1])
            #println(Ux, " // ", i, " // [", Q_raw[!, 1][i], ", ", Q_raw[!, 1][i+1], "]")
            if is_between(Ux, Q_raw[!, 1][i], Q_raw[!, 1][i+1]) 
                xk = i
                start_i = i
                break
            end
            i += 1
        end
        
        if xk == -1
            println("Error: EIS_get_Z_values: Zx ",Ux," is not in Q_raw[!, 1]: ",Q_raw[!, 1])
            return Exception()
        end
        Ik = Q_raw[!, 2][xk]
        Il = Q_raw[!, 2][xk+1]
        Ix = Ik + ((Il - Ik)/(Q_raw[!, 1][xk+1] - Q_raw[!, 1][xk])) * (Ux - Q_raw[!, 1][xk])
        
        Zx_values = [Zx_values; Ix]
    end
    
    return Zx_values
end

function linComb(SIM::EIS_simulation, Q1,Q2,lambda)
    if Q1[!, 1] == Q2[!, 1]
        return DataFrame(f = Q1.f, Z = Q1.Z.*(1-lambda) .+ Q2.Z.*(lambda))
    else
        println("ERROR: linComb_EIS: shape or value mismatch:\n",Q1.Z,"\n",Q2.Z)
        return Exception()
    end
end

function biliComb(SIM::EIS_simulation, Q11,Q12,Q21,Q22,x,y)
    if Q11[!, 1] == Q12[!, 1] && Q11[!, 1] == Q21[!, 1] && Q11[!, 1] ==Q22[!, 1]
        return DataFrame(f = Q11[!, 1], 
             Z = Q11[!, 2].*(1-x).*(1-y) .+ Q21[!, 2].*x.*(1-y) .+ Q12[!, 2].*(1-x).*y .+ Q22[!, 2].*x.*y
            )
    else
        println("ERROR: biliComb_EIS: shape mismatch or different *.f values")
        return Exception()
    end
end

function real_Bode_get_normalization_factor(a, b, A, B=1)
  # a, A is LF value of experimental and the reference Bode, respectively.
  # b, B is HF value of experimental and the reference Bode, respectively.
  # the following formula comes from the problem:
  #
  #               min_alpha ( (alpha*a - A)^2 + (alpha*b - B)^2 )
  return (a*A + b*B)/(a*a + b*b)
end

function fitnessFunction(SIM::EIS_simulation, sim_EIS::DataFrame, exp_EIS::DataFrame; relative_error=false)
        err = 0.0
        
        min_imag = Inf     
        
        sum_err_real = 0.0
        sum_err_imag = 0.0        
        
        if  exp_EIS.f == sim_EIS.f
                
                for row = 1:size(exp_EIS,1)
                    if (act_imag = imag(exp_EIS.Z[row])) < min_imag                      
                      min_imag = act_imag
                    end            
                    
                    sum_err_real += (
                                      (real(exp_EIS.Z[row]) - real(sim_EIS.Z[row]))
                                    )^2
                    sum_err_imag += (                        
                                      (act_imag - imag(sim_EIS.Z[row]))
                                    )^2 

#                     diff_real = abs(imag(exp_EIS.Z[row]) - imag(sim_EIS.Z[row]))
#                     diff_imag = abs(real(exp_EIS.Z[row]) - real(sim_EIS.Z[row]))
#                     
#                     diff_real > max_err_real ? max_err_real = diff_real : true
#                     diff_imag > max_err_imag ? max_err_imag = diff_imag : true
#                     
#                     err +=( diff_real^2 + diff_imag^2)
                end
                
                # normalization                
                real_factor = (
                                  real_Bode_get_normalization_factor(
                                                      real(exp_EIS.Z[1]), 
                                                      real(exp_EIS.Z[end]),
                                                      6   # LF value of the real reference EIS
                                                      )
                              )
                imag_factor = 1/(min_imag)    # EIS normalized such that the highest peak has height = 1
                #
                sum_err_imag *= imag_factor^2
                sum_err_real *= real_factor^2                              
                #
                err = sum_err_imag + sum_err_real 
        else
                println("ERROR: EIS_fitnesFunction: shape mismatch or different *.f values")
                return Exception()
        end
        # returns average error per checknode
        return sqrt(err)/Float32(size(exp_EIS,1))
#         return sqrt(err)/Float32(size(exp_EIS,1)) + max_err_real + max_err_imag
end

# function monotony_test_function(SIM::EIS_simulation, exp_EIS::DataFrame, sim_EIS::DataFrame)
#         err = 0.0
#         if  exp_EIS.f == sim_EIS.f
#                 
#                 for row = 1:size(exp_EIS,1)
#                         err +=( (real(exp_EIS.Z[row]) - real(sim_EIS.Z[row]))^2 
#                                  + (imag(exp_EIS.Z[row]) - imag(sim_EIS.Z[row]))^2)
#                 end
#         else
#                 println("ERROR: EIS_fitnesFunction: shape mismatch or different *.f values")
#                 return Exception()
#         end
#         # returns average error per checknode
#         return sqrt(err)/Float32(size(exp_EIS,1))
# end


function EIS_get_checknodes_short()
    # frequency nodes to compare
    # must be sorted upwards
    checknodes = (
        1.0e+0,
        1.0e+1,
        1.0e+2,
        1.0e+3
    )
end

function EIS_get_checknodes_geometrical(start_n, end_n, n_fac)
    # frequency nodes to compare
    # must be sorted upwards
    w_list = zeros(0)
    w = start_n
    if start_n == end_n
            append!(w_list, start_n)
            return w_list
    end
    if n_fac <= 1 || start_n > end_n
        println("ERROR: get_checknodes: n_fac <= 1 || start_n > end_n")
        return throw(Exception)
    end
    while w < end_n
        append!(w_list, w)
        w *= n_fac
    end    
    return w_list
end

function initialize_trend_tuples(SIM::EIS_simulation, EIS_ref::DataFrame)
  trend_tuples = DataFrame(prm_value=[])
  for i in 1:size(EIS_ref, 1)
    trend_tuples[!, Symbol(string(i))] = Array{Complex}(undef, 0)
  end
  return trend_tuples
end

function get_trend_tuple(SIM::EIS_simulation, EIS_ref::DataFrame, EIS_test::DataFrame)
  # plot signed deviation from referent simulation_curve
  trend_tuple = deepcopy(EIS_test.Z - EIS_ref.Z)
  return trend_tuple
#   return [
#     real(trend_tuple),
#     imag(trend_tuple),
#     abs(trend_tuple),
#     angle(trend_tuple)
#     ]
end

# function plot_trend_tuples(SIM::EIS_simulation, trend_tuples)
#   
# end



function EIS_crop_to_f_interval(EIS_exp, f_interval)
  return filter(row -> (f_interval[1] < row.f) && (row.f < f_interval[2]), EIS_exp)
end

function EIS_data_preprocessing(EIS_df)
  
  function get_lowest_freq_idx(EIS_df; find_at_least_negative=10)
    lowest_freq_idx = -1
    negative_counter = 0
    for (i, Z) in enumerate(EIS_df.Z)
      if imag(Z) < 0
        negative_counter += 1
      else
        negative_counter = 0
      end
      if negative_counter == find_at_least_negative
        lowest_freq_idx = i - find_at_least_negative + 1
        break
      end
    end    
    return lowest_freq_idx
  end
  
  #ysz_fitting.typical_plot_exp(EIS_simulation(800, 100, 0.0, use_DRT=false)[1], EIS_df, "!EIS_data_preprocessing")

  
  #lowest frequency cut off
  lowest_freq_idx = get_lowest_freq_idx(EIS_df, find_at_least_negative=10)
  if lowest_freq_idx == -1
    lowest_freq_idx = get_lowest_freq_idx(EIS_df, find_at_least_negative=4)
  end
  if lowest_freq_idx == -1
    println("ERROR: lowest_freq_idx not found!")
    
    lowest_freq_idx=1
    #return throw(Exception)
  end
  
  # intersection with x axis
  x_intersection_freq_idx = -1
  positive_counter = 0
  for i in (lowest_freq_idx + 5):length(EIS_df.f)
    if imag(EIS_df.Z[i]) > 0
      positive_counter += 1
    else  
      positive_counter = 0
    end
    if positive_counter == 8
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
    accepted_inductance_real_axis_threshold = 0.00*real(EIS_df.Z[lowest_freq_idx]) + 1.00*real(EIS_df.Z[x_intersection_freq_idx])
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

function f_interval_preprocessing(EIS_exp::DataFrame, f_interval)
  if f_interval!=Nothing
    if f_interval == "auto"
      #typical_plot_exp(SIM, EIS_exp, "! before")
      return EIS_exp = EIS_data_preprocessing(EIS_exp)
      #typical_plot_exp(SIM, EIS_exp, "! after")
    else
      return EIS_exp = EIS_crop_to_f_interval(EIS_exp, f_interval)
    end
  else
    return EIS_exp
  end
end

function EIS_add_inductance(EIS_df::DataFrame, L)
  EIS_out = deepcopy(EIS_df)
  for (i, f) in enumerate(EIS_df.f)
      EIS_out.Z[i] += im*2*pi*f*L
  end
  return EIS_out
end


