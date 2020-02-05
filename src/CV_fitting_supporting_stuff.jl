using CSV
using DataFrames
using PyPlot

include("./import_experimental_data.jl")

#########
# TODO
# - bug in CV_get_I_values >> nekrici, kdyz CV zacina az nad start_i
#########



## general functions
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

function linComb_CVs(CV1,CV2,lambda)
    if CV1.U == CV2.U
        return DataFrame(U = CV1.U, I = CV1.I.*(1-lambda) .+ CV2.I.*(lambda))
    else
        println("ERROR: linComb_CVs: shape or value mismatch:\n",CV1.U,"\n",CV2.U)
        return Exception()
    end
end

function biliComb_CVs(Q11,Q12,Q21,Q22,x,y)
    if Q11.U == Q12.U && Q11.U == Q21.U && Q11.U ==Q22.U
        return DataFrame(U = Q11.U, 
            I = Q11.I.*(1-x).*(1-y) .+ Q21.I*x.*(1-y) .+ Q12.I.*(1-x).*y .+ Q22.I.*x.*y
            )
    else
        println("ERROR: biliComb_CVs: shape mismatch or different *.U values")
        return Exception()
    end
end

function fitnessFunction(exp_CV::DataFrame, sim_CV::DataFrame)
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

function get_checknodes_short()
    # nodes of voltage U with direction of sweep where I will be fitted
    # [U direction; ... ]
    # the list must be sorted !!! w.r.t. U and direction 1 -> -1 -> 1
    # like in the real CV measurement
    checknodes = [
        0.05 1;
        0.25 1;
        0.5 1;
        0.75 1;
        0.95 1;
        0.95 -1;
        0.75 -1
        0.5 -1;
        0.25 -1;
        0.0 -1;
        -0.05 -1;
        -0.25 -1;
        -0.5 -1;
        -0.75 -1;
        -0.95 -1;
        -0.95 1;
        -0.75 1
        -0.5 1;
        -0.25 1;
        -0.05 1;
    ]
end

function get_checknodes(start_n,upper_n,lower_n,end_n,step_n)
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

function main()
        CVraw = import_CVtoDataFrame(T="700",pO2="00")
        checknodes = get_checknodes_short()
        exp_CV_nodes = DataFrame(
                U=checknodes[:,1], 
                I=CV_get_I_values(CVraw, checknodes), 
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
