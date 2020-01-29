using CSV
using DataFrames
using PyPlot

function import_CVtoDataFrame_path(f_name)
    df = DataFrame(U = Float32[], I = Float32[], t = Float32[])
    line_is_valid=false
    open(string(f_name,".cor")) do file
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
    return "TODO"
end


### TO DELETE ... ###############################
function is_between(x, a, b)
    if (b - a) > 0
        if (a <= x) & (x <= b)
            return true
        end
    else
        if (b <= x) & (x <= a)
            return true
        end
    end
    return false
end

function prepare_checkvalues(checknodes,CVraw,opt_pyplot=false)
    check_array = []
    Ux_values = []
    #println(CVraw.U)
    #println(CVraw.I)

    start_i = 1
    for row = 1:size(checknodes,1)
        direction = 1
        Ux = checknodes[row,1]
        xk = -1
        i = start_i
        while i < length(CVraw.U)-1
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
            print("Error: prepare_checkvalues: Ux ",Ux," is not in CVraw.U")
            return 0
        end
        Ik = CVraw.I[xk]
        Il = CVraw.I[xk+1]
        Ix = Ik + ((Il - Ik)/(CVraw.U[xk+1] - CVraw.U[xk])) * (Ux - CVraw.U[xk])
        
        append!(check_array,[[Ux,Ix,checknodes[row,2]]])
    end
    
    return check_array
end
################################################

function import_CVtoDataFrame(;T,pO2)
    fNAME=string("/home/masicko/Doktor_Linux/experimental_data/YSZ_09-2019_oxygen100/100 750to850 0to100%O2/",T,"C/100 ",T,"C ",pO2,"% do 1V/CV")
    
    
    return import_CVtoDataFrame_path(fNAME)
end


function main()    
    checknodes = [
        0.1 1;
        0.5 1;
        0.9 1;
        0.9 -1;
        0.5 -1
    ]
    CVraw = import_CVtoDataFrame(fNAME)
end
