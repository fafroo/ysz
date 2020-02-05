using CSV
using DataFrames
using PyPlot

include("./import_experimental_data.jl")

#########
# TODO
# - bug in CV_get_I_values >> nekrici, kdyz CV zacina az nad start_i
#########

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

function EIS_linComb(Q1,Q2,lambda)
    if Q1[!, 1] == Q2[!, 1]
        return DataFrame(f = Q1.f, Z = Q1.Z.*(1-lambda) .+ Q2.Z.*(lambda))
    else
        println("ERROR: linComb_EIS: shape or value mismatch:\n",Q1.Z,"\n",Q2.Z)
        return Exception()
    end
end

function EIS_biliComb(Q11,Q12,Q21,Q22,x,y)
    if Q11[!, 1] == Q12[!, 1] && Q11[!, 1] == Q21[!, 1] && Q11[!, 1] ==Q22[!, 1]
        return DataFrame(f = Q11[!, 1], 
             Z = Q11[!, 2].*(1-x).*(1-y) .+ Q21[!, 2].*x.*(1-y) .+ Q12[!, 2].*(1-x).*y .+ Q22[!, 2].*x.*y
            )
    else
        println("ERROR: biliComb_EIS: shape mismatch or different *.f values")
        return Exception()
    end
end

function EIS_fitnessFunction(exp_EIS::DataFrame, sim_EIS::DataFrame)
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

function main()
        checknodes = EIS_get_checknodes_geometrical(1, 8000, 2)
        
        EIS_raw2 = import_EIStoDataFrame(T="800",pO2="00",bias="0.0" )
        
        EIS_raw= import_EIStoDataFrame(T="800", pO2="00",bias="-0.5")
        #plot(real(EIS_raw2.Z), -imag(EIS_raw2.Z), "x-")      
        
        EIS_raw3= import_EIStoDataFrame(T="800", pO2="00",bias="-1.0")
        #plot(real(EIS_raw3.Z), -imag(EIS_raw3.Z), "x-")      
 
        exp_EIS_nodes = DataFrame(
                f=checknodes[:], 
                Z=EIS_get_Z_values(EIS_raw, checknodes), 
        )
        a = figure(1)
        Nyquist_plot(EIS_raw, "raw")
        Nyquist_plot(exp_EIS_nodes, "nodes")
        Nyquist_plot(exp_EIS_nodes, "nodes")

        
        b = figure(2)
        Nyquist_plot(EIS_raw, "baf")
        Nyquist_plot(exp_EIS_nodes, "krach")
        
        return
        
        exp_EIS_nodes = DataFrame(
                f = checknodes[:,1], 
                Z = CV_get_I_values(EISraw, checknodes), 
        )
        sim_CV_nodes = exp_CV_nodes[:, :]
        sim_CV_nodes.I.+=0.1
        
        println(exp_CV_nodes)
        println(sim_CV_nodes)
        
        fitnessFunction(exp_CV_nodes,sim_CV_nodes)
        #sim_CV_nodes = DataFrame(
        #	U=checknodes[:,1],
        #	I=CV_get_I_values(CVsim, checknodes), 
        #)
end
