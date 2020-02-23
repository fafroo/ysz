using CSV
using DataFrames
using PyPlot

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

function import_EIStoDataFrame_path(f_name)
    #df = DataFrame(f = Float32[], Re = Float32[], Im = Float32[])
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

function import_CVtoDataFrame(;TC,pO2)
  if pO2==0
    pO2="00"
  end
  fNAME=string("../snehurka/experimental_data_PSS/YSZ_09-2019_oxygen100/100 750to850 0to100%O2/",TC,"C/100 ",TC,"C ",pO2,"% do 1V/CV.cor")
  return import_CVtoDataFrame_path(fNAME)
end

function import_EIStoDataFrame(;TC, pO2, bias)
  if pO2==0
    pO2="00"
  end
  fNAME=string("../snehurka/experimental_data_PSS/YSZ_09-2019_oxygen100/100 750to850 0to100%O2/",TC,"C/100 ",TC,"C ",pO2,"% do 1V/is ",Float64(bias),"DC 50AC.z")
  return import_EIStoDataFrame_path(fNAME)
end

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

function Nyquist_plot(EIS_data, my_label="")
    title("Nyquist plot")
    xlabel("Re\$(Z)\$")
    ylabel("-Im\$(Z)\$")
   
   plot(real(EIS_data.Z), -imag(EIS_data.Z), "x-", label = my_label)
    
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
