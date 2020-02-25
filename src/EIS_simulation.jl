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
        name = string(this)
        
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
  
function get_shared_checknodes(sim::EIS_simulation)
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

function fitting_report(SIM::EIS_simulation, plot_prms_string, EIS_exp, EIS_sim)
  println("EIS_fitting error $(experiment_legend(SIM, latex=false)) $(plot_prms_string)  => ", fitnessFunction(SIM, EIS_exp, EIS_sim))
end


#######################################
# relatively static part ##############

function apply_checknodes(SIM::EIS_simulation, EIS_in,checknodes)
    DataFrame( f = checknodes, Z = EIS_get_Z_values(EIS_in, checknodes))
end

function load_file_prms(sim::EIS_simulation;save_dir, prms, prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), scripted_tuple=(0, 0, 0, 0, 0, 0), throw_exception=true, verbose=true)

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




