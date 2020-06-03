using CSV
using DataFrames
using PyPlot



const  CV_standard_figure_num = 5



######################################


mutable struct CV_simulation <: abstract_simulation
  TC::Float32
  pO2::Float32
  data_set::String
  #
  upp_bound::Float32
  low_bound::Float32
  dx_exp::Float64
  sample::Int64
  voltrate::Float32
  #
  fig_size::Tuple
  #
  checknodes::Any
  fitness_factor::Float64
  #
  name::String
  ID::Int16
  
  CV_simulation() = new()
end

function string(SIM::CV_simulation)
  return "CV_sim_TC_$(SIM.TC)_pO2_$(SIM.pO2)"
end

function CV_simulation(TC, pO2; data_set="MONO_110", upp_bound=1.0, low_bound=-1.0, dx_exp=-9, sample=8, voltrate=0.01, fig_size=(9, 6), fitness_factor=10.0)
  output = Array{abstract_simulation}(undef,0)
  for TC_item in TC
    for pO2_item in pO2
      this = CV_simulation()
      
      this.TC = TC_item
      this.pO2 = pO2_item 
      this.data_set = data_set
      #
      this.upp_bound = upp_bound
      this.low_bound = low_bound
      this.dx_exp = dx_exp
      this.sample = sample
      this.voltrate = voltrate
      #
      this.fig_size = fig_size
      #
      this.checknodes = get_shared_checknodes(this)
      this.fitness_factor = fitness_factor
      #
      this.name = "CV"
      this.ID = 1
      
      push!(output, this)
    end
  end
  return output
end

function CV_simulation(TC, pO2, bias; data_set="MONO_110", upp_bound=1.0, low_bound=-1.0, dx_exp=-9, sample=8, voltrate=0.01, fig_size=(9, 6), fitness_factor=10.0)
    CV_simulation(TC, pO2; data_set=data_set, dx_exp=dx_exp, sample=sample, voltrate=voltrate, fig_size=fig_size)
end


#######################################

function get_shared_checknodes(SIM::CV_simulation)
    return CV_get_checknodes(0.005,0.99,-0.99,-0.005,0.04)
end

function get_nice_checknodes(sim::CV_simulation)
  # experimental range is 0.00104 | 0.99990 | -0.9979 | -0.003.21
  return CV_get_checknodes(0.015,0.99,-0.99,-0.004,0.01)
end

function setting_legend(SIM::CV_simulation; latex=true)
  if latex
    return "\$\\theta=$(SIM.TC)\$°C \$\\mathrm{O}_2=$(SIM.pO2)\\%\$ $(SIM.data_set)"
  else
    return "TC=$(SIM.TC)°C pO2=$(SIM.pO2)% $(SIM.data_set)"
  end
end


function apply_checknodes(sim::CV_simulation, CV_in, checknodes)
    DataFrame( U = checknodes[:,1], I = CV_get_I_values(CV_in, checknodes))
end



function load_file_prms(sim::CV_simulation; save_dir, prms, prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), scripted_tuple=(0, 0, 0, 0, 0, 0), throw_exception=true, verbose=true)
    
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="CV", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    
    df_out = DataFrame()
    try throw_exception
      df_out = CSV.read(
          string(out_path,out_name),
      )
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

function save_file_prms(sim::CV_simulation, df_out, save_dir, prms, prms_names, scripted_tuple)
    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="CV", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    mkpath(out_path)

    CSV.write(
        string(out_path,out_name),
        df_out
    )
    return
end

function typical_plot_sim(SIM::CV_simulation, CV_df, additional_string="", to_standard_figure=true; plot_legend=true)
  if to_standard_figure
    figure(CV_standard_figure_num, figsize=SIM.fig_size)
  end
  my_label = "sim $(setting_legend(SIM))$(additional_string)"
  
  PyPlot.title("CV curve")
  PyPlot.xlabel(L"\eta \ [V]")
  PyPlot.ylabel(L"I \ [A]")
  
  PyPlot.plot(CV_df.U, CV_df.I, label=my_label)
  PyPlot.grid(true)
  if !(my_label == "") && plot_legend
      legend(loc="best")
  end
end

function typical_plot_exp(SIM::CV_simulation, CV_df, additional_string="", to_standard_figure=true)
  if to_standard_figure
    figure(CV_standard_figure_num, figsize=SIM.fig_size)
  end
  my_label = "exp $(setting_legend(SIM))$(additional_string)"
  
  PyPlot.title("CV curve")
  PyPlot.xlabel(L"\eta \ [V]")
  PyPlot.ylabel(L"I \ [A]")
  
  PyPlot.plot(CV_df.U, CV_df.I, "--", label=my_label)
  PyPlot.grid(true)

  if !(my_label == "")
      legend(loc="best")
  end
end

function typical_run_simulation(SIM::CV_simulation, prms_names_in, prms_values_in, pyplot::Int=0) 
  ysz_experiments.run_new(
      out_df_bool=true, voltammetry=true, pyplot=(pyplot == 2 ? true : false), 
      dx_exp=SIM.dx_exp, sample=SIM.sample, upp_bound=SIM.upp_bound, low_bound=SIM.low_bound, voltrate=SIM.voltrate,
      T=TCtoT(SIM.TC), pO2=pO2tosim(SIM.pO2),
      prms_names_in=prms_names_in,
      prms_values_in=prms_values_in,
  )
end

function import_data_to_DataFrame(SIM::CV_simulation)
  import_CVtoDataFrame(TC=SIM.TC, pO2=SIM.pO2, data_set=SIM.data_set)
end

function CV_view_experimental_data(TC_list, pO2_list; use_checknodes=false, fig_num=11)    
    figure(fig_num)
    for TC in TC_list
      for pO2 in pO2_list
        if use_checknodes
          checknodes =  CV_get_shared_checknodes()
          CV_exp = CV_apply_checknodes(import_CVtoDataFrame(TC=TC, pO2=pO2), checknodes)
        else
          CV_exp = import_CVtoDataFrame(TC=TC, pO2=pO2)
        end
        typical_plot_exp(CV_simulation(TC, pO2)..., CV_exp, "", false)
      end
    end
end

function fitness_error_report(SIM::CV_simulation, plot_prms_string, CV_exp, CV_sim)
  println("CV fitness error $(setting_legend(SIM, latex=false)) $(plot_prms_string)  => ", fitnessFunction(SIM, CV_exp, CV_sim))
end

function CV_get_I_values(CVraw, checknodes)
    Ix_values = []
    #println(CVraw.U)
    #println(CVraw.I)

    start_i = 1
    for row = 1:size(checknodes,1)
        direction = 1
        Ux = checknodes[row,1]
        xk = -1
        i = start_i
        while i < length(CVraw.U)
            #println(Ux, " // ", i, " // [", CVraw.U[i], ", ", CVraw.U[i+1], "]")
            if direction*(CVraw.U[i+1] - CVraw.U[i]) < 0
                direction *= -1
            end
            if is_between(Ux, CVraw.U[i], CVraw.U[i+1]) & (direction == checknodes[row,2]) 
                xk = i
                start_i = i
                break
            end
            i += 1
        end
        
        if xk == -1
            println("Error: CV_get_I_values: Ux ",Ux," is not in CVraw.U: ",CVraw.U)
            return Exception()
        end
        Ik = CVraw.I[xk]
        Il = CVraw.I[xk+1]
        Ix = Ik + ((Il - Ik)/(CVraw.U[xk+1] - CVraw.U[xk])) * (Ux - CVraw.U[xk])
        
        Ix_values = [Ix_values; Ix]
        end
    
    return Ix_values
end

function linComb(SIM::CV_simulation, CV1,CV2,lambda)
    if CV1.U == CV2.U
        return DataFrame(U = CV1.U, I = CV1.I.*(1-lambda) .+ CV2.I.*(lambda))
    else
        println("ERROR: linComb_CVs: shape or value mismatch:\n",CV1.U,"\n",CV2.U)
        return Exception()
    end
end

function biliComb(SIM::CV_simulation, Q12,Q21,Q22,x,y)
    if Q11.U == Q12.U && Q11.U == Q21.U && Q11.U ==Q22.U
        return DataFrame(U = Q11.U, 
            I = Q11.I.*(1-x).*(1-y) .+ Q21.I*x.*(1-y) .+ Q12.I.*(1-x).*y .+ Q22.I.*x.*y
            )
    else
        println("ERROR: biliComb_CVs: shape mismatch or different *.U values")
        return Exception()
    end
end

function fitnessFunction(SIM::CV_simulation, exp_CV::DataFrame, sim_CV::DataFrame)
        if (count=size(exp_CV,1)) == size(sim_CV,1)
                err = 0.0
                for row = 1:count
                        err += (exp_CV.I[row] - sim_CV.I[row])^2
                end
        else
                println("ERROR: fitnesFunction: shape mismatch")
                return Exception()
        end
        # returns average error per checknode
        return sqrt(err)/Float32(count)
end

function CV_get_checknodes(start_n,upper_n,lower_n,end_n,step_n)
    # nodes of voltage U with direction of sweep where I will be fitted
    # [U direction; ... ]
    # the list must be sorted !!! w.r.t. U and direction 1 -> -1 -> 1
    # like in the real CV measurement
    ll = start_n : step_n : upper_n
    res = [ll [1 for i in 1:size(ll,1)]]
    ll = upper_n : -step_n : lower_n
    res = vcat(res, [ll [-1 for i in 1:size(ll,1)]])
    ll = lower_n : step_n : end_n
    res = vcat(res, [ll [1 for i in 1:size(ll,1)]])
    return res
end

function initialize_trend_tuples(SIM::CV_simulation, CV_ref::DataFrame)
  trend_tuples = DataFrame(prm_value=[])
  for i in 1:size(CV_ref,1)
    trend_tuples[!, Symbol(string(i))] = Array{Float32}(undef, 0)
  end
  return trend_tuples
end

function get_trend_tuple(SIM::CV_simulation, CV_ref::DataFrame, CV_test::DataFrame)
  # plot signed deviation from referent simulation_curve
  trend_tuple = deepcopy(CV_test.I - CV_ref.I)
  return trend_tuple
end

# function plot_trend_tuples(SIM::CV_simulation, trend_tuples)
#   
# end

