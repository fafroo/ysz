using CSV
using DataFrames
using PyPlot

#include("../src/incommon_simulation.jl")


  




######################################


mutable struct CV_simulation <: abstract_simulation
  TC::Int64
  pO2::Int64
  #
  dx_exp::Float64
  sample::Int64
  fig_size::Tuple
  #
  name::String
  
  CV_simulation() = new()
end

function string(SIM::CV_simulation)
  return "CV_sim_TC_$(SIM.TC)_pO2_$(SIM.pO2)"
end

function CV_simulation(TC, pO2; dx_exp=-9, sample=8, fig_size=(9, 6))
  output = Array{abstract_simulation}(undef,0)
  for TC_item in TC
    for pO2_item in pO2
      this = CV_simulation()
      
      this.TC = TC_item
      this.pO2 = pO2_item 
      #
      this.dx_exp = dx_exp
      this.sample = sample
      this.fig_size = fig_size
      #
      name = string(this)
      
      push!(output, this)
    end
  end
  return output
end


#######################################

function get_shared_checknodes(sim::CV_simulation)
    return CV_get_checknodes(0.005,0.99,-0.99,-0.005,0.04)
end

function get_nice_checknodes(sim::CV_simulation)
  # experimental range is 0.00104 | 0.99990 | -0.9979 | -0.003.21
  return CV_get_checknodes(0.015,0.99,-0.99,-0.004,0.01)
end

function experiment_legend(SIM::CV_simulation; latex=true)
  if latex
    return "\$\\theta=$(SIM.TC)\$°C \$\\mathrm{O}_2=$(SIM.pO2)\\%\$"
  else
    return "TC=$(SIM.TC)°C pO2=$(SIM.pO2)%"
  end
end


function apply_checknodes(sim::CV_simulation, CV_in, checknodes)
    DataFrame( U = checknodes[:,1], I = CV_get_I_values(CV_in, checknodes))
end



function load_file_prms(sim::CV_simulation;save_dir, prms, prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), scripted_tuple=(0, 0, 0, 0, 0, 0), throw_exception=true, verbose=true)
    
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
    run(`mkdir -p $out_path`)

    CSV.write(
        string(out_path,out_name),
        df_out
    )
    return
end

function typical_plot_sim(SIM::CV_simulation, CV_df, additional_string="")
  figure(5, figsize=SIM.fig_size)
  my_label = "sim $(experiment_legend(SIM))$(additional_string)"
  
  PyPlot.title("CV curve")
  PyPlot.xlabel(L"\eta \ (V)")
  PyPlot.ylabel(L"I \ (A)")
  
  PyPlot.plot(CV_df.U, CV_df.I, label=my_label)

  if !(my_label == "")
      legend(loc="best")
  end
end

function typical_plot_exp(SIM::CV_simulation, CV_df, additional_string="")
  figure(5, figsize=SIM.fig_size)
  my_label = "exp $(experiment_legend(SIM))$(additional_string)"
  
  PyPlot.title("CV curve")
  PyPlot.xlabel(L"\eta \ (V)")
  PyPlot.ylabel(L"I \ (A)")
  
  PyPlot.plot(CV_df.U, CV_df.I, "--", label=my_label)

  if !(my_label == "")
      legend(loc="best")
  end
end

function typical_run_simulation(SIM::CV_simulation, prms_names_in, prms_values_in, pyplot::Int=0) 
  ysz_experiments.run_new(
      out_df_bool=true, voltammetry=true, dx_exp=SIM.dx_exp, sample=SIM.sample, pyplot=(pyplot == 2 ? true : false),
      T=TCtoT(SIM.TC), pO2=pO2tosim(SIM.pO2),
      prms_names_in=prms_names_in,
      prms_values_in=prms_values_in,
  )
end

function import_data_to_DataFrame(SIM::CV_simulation)
  import_CVtoDataFrame(TC=SIM.TC, pO2=SIM.pO2)
end

function CV_view_experimental_data(TC_list, pO2_list; use_checknodes=false, fig_num=20)    
    figure(fig_num)
    for TC in TC_list
      for pO2 in pO2_list
        if use_checknodes
          checknodes =  CV_get_shared_checknodes()
          CV_exp = CV_apply_checknodes(import_CVtoDataFrame(TC=TC, pO2=pO2), checknodes)
        else
          CV_exp = import_CVtoDataFrame(TC=TC, pO2=pO2)
        end
        CV_plot(CV_exp, "exp $(experiment_legend(CV_simulation(TC, pO2)))")
      end
    end
end

function fitting_report(SIM::CV_simulation, plot_prms_string, CV_exp, CV_sim)
  println("CV_fitting error $(experiment_legend(SIM, latex=false)) $(plot_prms_string)  => ", fitnessFunction(SIM, CV_exp, CV_sim))
end
