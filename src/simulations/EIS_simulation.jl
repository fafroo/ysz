using CSV
using DataFrames
using PyPlot



const EIS_standard_figure_num = 6

include("../DRT.jl")


######################################

mutable struct EIS_simulation <: abstract_simulation
  TC::Float64
  pO2::Float64
  bias::Float64
  #
  dx_exp::Float64
  f_range::Tuple
  fig_size::Tuple 
  #
  checknodes::Any
  fitness_factor::Float64
  #
  use_DRT::Bool
  DRT_lambda::Float64
  #
  name::String
  ID::Int16
  
  EIS_simulation() = new()
end

function string(SIM::EIS_simulation)
  return "EIS_sim_TC_$(SIM.TC)_pO2_$(SIM.pO2)_bias_$(SIM.bias)"
end

function EIS_simulation(TC, pO2, bias=0.0; dx_exp=-9, f_range=EIS_get_shared_f_range(), fig_size=(9, 6), use_DRT=true, DRT_lambda=0.0)
  output = Array{abstract_simulation}(undef,0)
  for TC_item in TC
    for pO2_item in pO2
      for bias_item in bias
        this = EIS_simulation()
        #
        this.TC = TC_item
        this.pO2 = pO2_item
        this.bias = bias_item
        #
        this.dx_exp = dx_exp
        this.f_range = f_range
        this.fig_size = fig_size
        #
        this.checknodes = get_shared_checknodes(this)
        this.fitness_factor = 1.0
        #
        this.use_DRT = use_DRT
        this.DRT_lambda = DRT_lambda
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
    return f_range = (1.1, 65000, 1.2)
end
  
function get_shared_checknodes(SIM::EIS_simulation)
    return EIS_get_checknodes_geometrical(EIS_get_shared_f_range()...)
end

function setting_legend(SIM::EIS_simulation; latex=true)
  if latex
    return "\$\\theta=$(SIM.TC)\$°C \$\\mathrm{O}_2=$(SIM.pO2)\\% \$ \$\\mathrm{bias}=$(SIM.bias)\$"
  else
    return "TC=$(SIM.TC)°C pO2=$(SIM.pO2)% bias=$(SIM.bias)V"
  end
end

# function typical_plot(SIM::EIS_simulation, EIS_df, additional_string, line_style)
# 
# end

function typical_plot_sim(SIM::EIS_simulation, EIS_df, additional_string="", to_standard_figure=true)
  if to_standard_figure
    figure(EIS_standard_figure_num, figsize=SIM.fig_size)
  end
  
  if length(additional_string)>0 && additional_string[1]=='!'
    my_label = additional_string[2:end]
  else
    my_label = "sim $(setting_legend(SIM))$(additional_string)"
  end

  s1 = subplot(221)
  title("Nyquist plot")
  xlabel("Re\$(Z) \\ [\\Omega]\$")
  ylabel("-Im\$(Z) \\ [\\Omega]\$")
  
  plot(real(EIS_df.Z), -imag(EIS_df.Z), "x-", label = my_label)
  
  if !(my_label == "")
      legend(loc="best")
  end
  grid(true)
  s1.set_aspect(1.0)
  
  s2 = subplot(322)
  title("Bode plot - Re")
  #xlabel("log f \$(\$Hz\$)\$")
  ylabel("Re\$(Z) \\ [\\Omega]\$")
  plot(log10.(EIS_df.f), real(EIS_df.Z), "x-", label = my_label)
    
  s3 = subplot(324)
  title("Bode plot - Im")
  xlabel("log f \$(\$Hz\$)\$")
  ylabel("-Im\$(Z) \\ [\\Omega]\$")
  plot(log10.(EIS_df.f), -imag(EIS_df.Z), "x-", label = my_label)
  
  
  if SIM.use_DRT
    DRT_actual = get_DRT(EIS_df, SIM.DRT_lambda)
    
    s4 = subplot(325)
    plot_DRT_h(DRT_actual, false)
    
    s5 = subplot(326)
    plot_DRT_RC(DRT_actual, false)
  else
    s4 = subplot(325)
    plot([])
    
    s5 = subplot(326)
    plot([])
  end
  
  
  
  plt.subplots_adjust(bottom=0.07, top=0.95)
end

function typical_plot_exp(SIM::EIS_simulation, EIS_df, additional_string="", to_standard_figure=true)
  if to_standard_figure
    figure(EIS_standard_figure_num, figsize=SIM.fig_size)
  end
  
  if length(additional_string)>0 && additional_string[1]=='!'
    my_label = additional_string[2:end]
  else
    my_label = "exp $(setting_legend(SIM))$(additional_string)"
  end
 
  
  s1 = subplot(221)
  title("Nyquist plot")
  xlabel("Re\$(Z) \\ [\\Omega]\$")
  ylabel("-Im\$(Z) \\ [\\Omega]\$")
  
  plot(real(EIS_df.Z), -imag(EIS_df.Z), "x:", label = my_label)
  
  if !(my_label == "")
      legend(loc="best")
  end
  grid(true)
  s1.set_aspect(1.0)
  
  s2 = subplot(322)
  title("Bode plot - Re")
  #xlabel("log f \$(\$Hz\$)\$")
  ylabel("Re\$(Z) \\ [\\Omega]\$")
  plot(log10.(EIS_df.f), real(EIS_df.Z), "x:", label = my_label)
    
  s3 = subplot(324)
  title("Bode plot - Im")
  xlabel("log f \$(\$Hz\$)\$")
  ylabel("-Im\$(Z) \\ [\\Omega]\$")
  plot(log10.(EIS_df.f), -imag(EIS_df.Z), "x:", label = my_label)
  
  if SIM.use_DRT
    DRT_actual = get_DRT(EIS_df, SIM.DRT_lambda)
    
    s4 = subplot(325)
    plot_DRT_h(DRT_actual, false)
    
    s5 = subplot(326)
    plot_DRT_RC(DRT_actual, false)
  else
    s4 = subplot(325)
    plot([])
    
    s5 = subplot(326)
    plot([])
  end
  
  plt.subplots_adjust(bottom=0.07, top=0.95)
end

function fitness_error_report(SIM::EIS_simulation, plot_prms_string, EIS_exp, EIS_sim)
  println("EIS fitness error $(setting_legend(SIM, latex=false)) $(plot_prms_string)  => ", fitnessFunction(SIM, EIS_exp, EIS_sim))
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
      df_out = DataFrame(f = df.f, Z = (df.Re .+ df.Im.*im))
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

function save_file_prms(sim::EIS_simulation, df_out, save_dir, prms, prms_names, scripted_tuple)
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="EIS", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    run(`mkdir -p $out_path`)

    CSV.write(
        string(out_path,out_name),
        DataFrame(f = df_out.f, Re = real(df_out.Z), Im = imag(df_out.Z))
    )
    return
end


function typical_run_simulation(SIM::EIS_simulation, prms_names_in, prms_values_in, pyplot::Int=0)
  EIS_df = ysz_experiments.run_new(
      pyplot=(pyplot == 2 ? true : false), EIS_IS=true, out_df_bool=true, EIS_bias=SIM.bias, f_range=SIM.f_range,
      dx_exp=SIM.dx_exp,
      T=TCtoT(SIM.TC), pO2=pO2tosim(SIM.pO2),
      prms_names_in=prms_names_in,
      prms_values_in=prms_values_in
  )
end

function import_data_to_DataFrame(SIM::EIS_simulation)
  import_EIStoDataFrame(TC=SIM.TC, pO2=SIM.pO2, bias=SIM.bias)
end

function EIS_test_checknodes_range(f_range=EIS_get_shared_f_range())
  Nyquist_plot(apply_checknodes(EIS_simulation(), import_EIStoDataFrame(TC=800,pO2=100,bias=0.0),EIS_get_checknodes_geometrical(f_range...)))
end


function EIS_view_experimental_data(TC_list, pO2_list, bias_list; use_checknodes=false, fig_num=12)    
    
    for TC in TC_list
      for pO2 in pO2_list
        for bias in bias_list
          if use_checknodes
            SIM = EIS_simulation()
            checknodes =  get_shared_checknodes(SIM)
            EIS_exp = apply_checknodes(SIM, import_EIStoDataFrame(TC=TC, pO2=pO2, bias=bias), checknodes)
          else
            EIS_exp = import_EIStoDataFrame(TC=TC, pO2=pO2, bias=bias)
          end
          figure(fig_num)
          typical_plot_exp(EIS_simulation(TC, pO2, bias)..., EIS_exp, "", false)
        end
      end
    end
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


function fitnessFunction(SIM::EIS_simulation, exp_EIS::DataFrame, sim_EIS::DataFrame)
        err = 0.0
        if  exp_EIS.f == sim_EIS.f
                
                for row = 1:size(exp_EIS,1)
                        err +=( (real(exp_EIS.Z[row]) - real(sim_EIS.Z[row]))^2 
                                 + (imag(exp_EIS.Z[row]) - imag(sim_EIS.Z[row]))^2)
                end
        else
                println("ERROR: EIS_fitnesFunction: shape mismatch or different *.f values")
                return Exception()
        end
        # returns average error per checknode
        return sqrt(err)/Float32(size(exp_EIS,1))
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





