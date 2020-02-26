using CSV
using DataFrames
using PyPlot







######################################

mutable struct EIS_simulation <: abstract_simulation
  TC::Int64
  pO2::Int64
  bias::Float64
  #
  dx_exp::Float64
  omega_range::Tuple
  fig_size::Tuple 
  #
  checknodes::Any
  #
  name::String
  
  EIS_simulation() = new()
end

function string(SIM::EIS_simulation)
  return "EIS_sim_TC_$(SIM.TC)_pO2_$(SIM.pO2)_bias_$(SIM.bias)"
end

function EIS_simulation(TC, pO2, bias=0.0; dx_exp=-9, omega_range=EIS_get_shared_omega_range(), fig_size=(9, 6))
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
        this.omega_range = omega_range
        this.fig_size = fig_size
        #
        this.checknodes = get_shared_checknodes(this)
        #
        this.name = string(this)
        
        push!(output, this)
      end
    end
  end
  return output
end

#######################################
# potentially changable stuff #########

function EIS_get_shared_omega_range()
    # experimental range is 0.1 - 65000 Hz
    # omegas = (w0, w1, w_fac)
    return omega_range = (1.1, 51000, 1.2)
end
  
function get_shared_checknodes(SIM::EIS_simulation)
    return EIS_get_checknodes_geometrical(EIS_get_shared_omega_range()...)
end

function experiment_legend(SIM::EIS_simulation; latex=true)
  if latex
    return "\$\\theta=$(SIM.TC)\$°C \$\\mathrm{O}_2=$(SIM.pO2)\\% \$ \$\\mathrm{bias}=$(SIM.bias)\$"
  else
    return "TC=$(SIM.TC)°C pO2=$(SIM.pO2)% bias=$(SIM.bias)V"
  end
end


function typical_plot_sim(SIM::EIS_simulation, EIS_df, additional_string="")
  figure(6, figsize=SIM.fig_size)
  my_label = "exp $(experiment_legend(SIM))$(additional_string)"

  title("Nyquist plot")
  xlabel("Re\$(Z)\$")
  ylabel("-Im\$(Z)\$")
  
  plot(real(EIS_df.Z), -imag(EIS_df.Z), "x-", label = my_label)
  
  if !(my_label == "")
      legend(loc="best")
  end
  grid(true)
  ax = gca()
  ax.set_aspect(1.0)
end

function typical_plot_exp(SIM::EIS_simulation, EIS_df, additional_string="")
  figure(6, figsize=SIM.fig_size)
  my_label = "sim $(experiment_legend(SIM))$(additional_string)"

  title("Nyquist plot")
  xlabel("Re\$(Z)\$")
  ylabel("-Im\$(Z)\$")
  
  plot(real(EIS_df.Z), -imag(EIS_df.Z), "x:", label = my_label)
  
  if !(my_label == "")
      legend(loc="best")
  end
  grid(true)
  ax = gca()
  ax.set_aspect(1.0)  
end

function fitness_error_report(SIM::EIS_simulation, plot_prms_string, EIS_exp, EIS_sim)
  println("EIS fitness error $(experiment_legend(SIM, latex=false)) $(plot_prms_string)  => ", fitnessFunction(SIM, EIS_exp, EIS_sim))
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
      pyplot=(pyplot == 2 ? true : false), EIS_IS=true, out_df_bool=true, EIS_bias=SIM.bias, omega_range=SIM.omega_range,
      dx_exp=SIM.dx_exp,
      T=TCtoT(SIM.TC), pO2=pO2tosim(SIM.pO2),
      prms_names_in=prms_names_in,
      prms_values_in=prms_values_in,
  )
end

function import_data_to_DataFrame(SIM::EIS_simulation)
  import_EIStoDataFrame(TC=SIM.TC, pO2=SIM.pO2, bias=SIM.bias)
end

function EIS_test_checknodes_range(omega_range=EIS_get_shared_omega_range())
  Nyquist_plot(EIS_apply_checknodes(import_EIStoDataFrame(TC=800,pO2=100,bias=0.0),EIS_get_checknodes_geometrical(omega_range...)))
end


function EIS_view_experimental_data(TC_list, pO2_list, bias_list; use_checknodes=false, fig_num=10)    
    figure(fig_num)
    for TC in TC_list
      for pO2 in pO2_list
        for bias in bias_list
          if use_checknodes
            checknodes =  EIS_get_shared_checknodes()
            EIS_exp = EIS_apply_checknodes(import_EIStoDataFrame(TC=TC, pO2=pO2, bias=bias), checknodes)
          else
            EIS_exp = import_EIStoDataFrame(TC=TC, pO2=pO2, bias=bias)
          end
          Nyquist_plot(EIS_exp, "exp $(experiment_legend(EIS_simulation(TC, pO2, bias=bias)))")
        end
      end
    end
end

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


