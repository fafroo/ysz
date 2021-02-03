using CSV
using DataFrames
using PyPlot
using Printf

function import_CVtoDataFrame_path(f_name)
    df = DataFrame(U = Float32[], I = Float32[], t = Float32[])
    line_is_valid=false
    open(f_name) do file
        for ln in eachline(file)
            if line_is_valid
                push!(df,[parse(Float32,el) for el in split(ln)])
            else
                if ln=="End Comments"
                    line_is_valid=true
                end
            end
        end
    end
    df
end

#  "22-10-2020-point_2_PROKOP.txt"

function import_EIStoDataFrame_path(f_name)
    #df = DataFrame(f = Float32[], Re = Float32[], Im = Float32[])
    if length(f_name) >= 10 && f_name[end-9:end]=="PROKOP.txt"
      # position of frequency column
      freq_col_idx = parse(Int32,f_name[end-10])
      
      df = CSV.File(f_name, decimal=',') |> DataFrame            
      return DataFrame(
        f = reverse(deepcopy(df[!, freq_col_idx])), 
        Z = reverse(deepcopy(df[!, freq_col_idx+1] .- im*df[!, freq_col_idx+2]))
      )
    end
    df = DataFrame(f = Float32[], Z = Complex[])
    line_is_valid=false
    open(f_name) do file
        for ln in eachline(file)
            if line_is_valid
                splitted_row = [parse(Float32,el) for el in split(ln)]
                #push!(df, (splitted_row[1], splitted_row[5], splitted_row[6] ))
                push!(df, (splitted_row[1], splitted_row[5] +  splitted_row[6]*im ))
            else
                if ln=="End Comments"
                    line_is_valid=true
                end
            end
        end
    end
    df.f = reverse(df.f)
    df.Z = reverse(df.Z)    
    return df
end

function get_U_I_from_IVpoint(f_name)
    line_is_valid=false
    last_line = []
    open(f_name) do file
        for ln in eachline(file)
            if line_is_valid
                last_line = [parse(Float32,el) for el in split(ln)]
            else
                if ln=="End Comments"
                    line_is_valid=true
                end
            end
        end
    end
    return (last_line[1], last_line[2])
end

function import_IVtoDataFrame_folder(;TC, pO2, bias_array, folder)
  df = DataFrame(U = Float32[], I = Float32[], t = Float32[])
  for bias in bias_array
    if bias == 0
      push!(df, (0,0,0)) 
    else
      push!(df,(get_U_I_from_IVpoint(folder*"/pol$(Int(bias*1000)).cor")..., 0) ) 
    end
  end
  return df
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

function import_CVtoDataFrame(;TC,pO2, data_set="MONO_110")
  pO2=Int64(pO2)
  pO2_Dan = @sprintf("%02g",pO2)
  #
  TC=Int64(TC)
  #
  if data_set=="OLD_MONO_100"
    fNAME=string("../snehurka/experimental_data_PSS/YSZ_09-2019_oxygen100/100 750to850 0to100%O2/",TC,"C/100 ",TC,"C ",pO2_Dan,"% do 1V/CV.cor")
  elseif data_set == "POLY"
    fNAME=string("../snehurka/experimental_data_PSS/jako asi 6/$(TC) $(pO2_Dan) 6/cv1.cor")
  elseif data_set == "MONO_110"
    fNAME=string("../snehurka/experimental_data_PSS/YSZ 110/110 $(TC) $(pO2_Dan)/cv1.cor")
  #
  #
  elseif data_set=="POLY_I-V"
    # this needs to be added to separate simulation !!! ... IV_simulation
    return import_IVtoDataFrame_folder(TC=TC, pO2=pO2, bias_array=vcat(collect(0 : 0.1 : 1), collect(0.9 : -0.1 : -0.9), collect(-1 : 0.1 : 0)), 
          folder="../snehurka/experimental_data_PSS/jako asi 6/$(TC) $(pO2_Dan) 6/")
  elseif data_set=="MONO_110_I-V"
    return import_IVtoDataFrame_folder(TC=TC, pO2=pO2, bias_array=vcat(collect(0 : 0.1 : 1), collect(0.9 : -0.1 : -0.9), collect(-1 : 0.1 : 0)), 
          folder="../snehurka/experimental_data_PSS/YSZ 110/110 $(TC) $(pO2_Dan)/")
  else
    fNAME=string("../snehurka/experimental_data_PSS/individual_files/$(data_set)")
  end
  return import_CVtoDataFrame_path(fNAME)
end





function import_EIStoDataFrame(;TC, pO2, bias, data_set="MONO_110", extra_tokens=Dict())
  ET = extra_tokens
  #
  pO2=Int64(pO2)
  pO2_Dan = @sprintf("%02g",pO2)
  pO2_Michal = @sprintf("%03g",pO2)
  #
  TC=Int64(TC)
  #
  bias=Float64(bias)
    bias_mV = Int32(bias*1000)
  if abs(bias_mV) < 10
    bias_mV = "0_1"
  end
  #
  if length(data_set) > 3 && (data_set[end-3 : end-1] == "OCV")
    bias_mv = "0_"*data_set[end]
  end
  #
  if length(data_set) >= 4 && data_set[1:4]=="POLY"
    fNAME=string("../snehurka/experimental_data_PSS/jako asi 6/$(TC) $(pO2_Dan) 6/eis_$(bias_mV).z") 
  elseif data_set=="Zahner"
    fNAME=string("../snehurka/experimental_data_PSS/individual_files/TEST DRT - Zahner - dummy cell.z")
  elseif data_set=="Solartron"
    fNAME=string("../snehurka/experimental_data_PSS/individual_files/TEST DRT - Solartron - dummy cell.z")
  elseif data_set=="HebbWagner"
    # TC \in (600 : 20 : 720) ... bias = 0.3 ... pO2 = nizke, temer nulove
    fNAME=string("../snehurka/experimental_data_PSS/HebbWagner/$(TC) C/$(TC)_EIS $(bias)V v ref 50mV amplituda.z")
  
  ##########################################################
  # ! ! ! ! ! pO2 used as a repetition number ! ! ! ! ! ! 
  elseif data_set=="HebbWagner_110_more_biases"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("../snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/02 YSZ 110 Ag-Ag/700 a 800 C vice biasu/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
  
  elseif data_set=="HebbWagner_100"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("../snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/04 YSZ 100 Ag-Ag/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
    
  elseif data_set=="HebbWagner_110"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("../snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/02 YSZ 110 Ag-Ag/vsechny teploty a jeden bias/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
    
  elseif data_set=="HebbWagner_111"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("../snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/03 YSZ 111 Ag-Ag/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(pO2).z")
      
      
      
  elseif data_set=="HW_TEST"
    # TC \in (600 : 20 : 720) ...
    fNAME=string("../snehurka/experimental_data_PSS/Hebb-Wagner_monokrystaly/03 YSZ 111 Ag-Ag/$(TC) C/EIS $(bias)V v ref 50mV amplituda_Rp0$(ET["rep"]).z")
      
      
      
      
      
      
  ##########################################################
  
  
  elseif length(data_set) >= 8 && data_set[1:8]=="MONO_110"
    fNAME=string("../snehurka/experimental_data_PSS/YSZ 110/110 $(TC) $(pO2_Dan)/eis_$(bias_mV).z")
  elseif data_set=="OLD_MONO_100"
    fNAME=string("../snehurka/experimental_data_PSS/YSZ_09-2019_oxygen100/100 750to850 0to100%O2/$(TC)C/100 $(TC)C $(pO2_Dan)% do 1V/is $(bias)DC 50AC.z")
  elseif data_set=="MICHAL"
    if bias == 0
      ocp_token = "_ocp"
    else
      ocp_token = ""
    end
    fNAME=string("../snehurka/experimental_data_PSS/K4/$(TC) C/$(pO2_Michal) O2/EIS_$(bias)DC_50AC$(ocp_token).z")
  
  elseif data_set=="K06"
    if bias == 0
      ocp_token = "_ocp"
    else
      ocp_token = ""
    end
    fNAME=string("../snehurka/experimental_data_PSS/K06/$(TC) C/$(pO2_Michal) O2/EIS_$(bias)DC_50AC$(ocp_token).z")
  
  elseif data_set=="K07"
    if bias == 0
      ocp_token = "_ocp"
    else
      ocp_token = ""
    end
    fNAME=string("../snehurka/experimental_data_PSS/K07/$(TC) C/$(pO2_Michal) O2/EIS_$(bias)DC_50AC$(ocp_token).z")  
  else
    fNAME=string("../snehurka/experimental_data_PSS/individual_files/$(data_set)")
  end
  
  return import_EIStoDataFrame_path(fNAME)
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

function CV_plot(CV_df, my_label="")
#    plot(df.U, df.Ib ,"blue"    ,label="bulk")
    #plot(phi_range[cv_range].-phi0, Ibb_range[cv_range] ,label="bulk_grad")
#    plot(df.U, df.Is, "green"   ,label="surf")
#    plot(df.U, df.Ir, "red"     ,label="reac")    
    PyPlot.title("CV curve")
    PyPlot.xlabel(L"\eta \ (V)")
    PyPlot.ylabel(L"I \ (A)")
    
    PyPlot.plot(CV_df.U, CV_df.I, label=my_label)

    if !(my_label == "")
        legend(loc="best")
    end
end

function Nyquist_plot(EIS_df, my_label="")
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

function main()    
    CVraw = import_CVtoDataFrame(fNAME)
end
