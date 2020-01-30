module ysz_fitting

using PyPlot
using DataFrames

include("../examples/ysz_voltammetry_like.jl")
include("../src/fitting_stuff.jl")
include("../src/import_experimental_data.jl")




function simple_run()
    old_prms = [21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1]
    ysz_voltammetry_like.run_new(
            out_df_bool=true, voltammetry=true, sample=40, pyplot_finall=true,
            prms_in=[19.50, 19.00, 0.0905748, -0.708014, 0.6074566741435283, -1.5]
            )
    return
end


function finall_plot(df, label)
#    plot(df.U, df.Ib ,"blue"    ,label="bulk")
    #plot(phi_range[cv_range].-phi0, Ibb_range[cv_range] ,label="bulk_grad")
#    plot(df.U, df.Is, "green"   ,label="surf")
#    plot(df.U, df.Ir, "red"     ,label="reac")    
    PyPlot.plot(df.U, df.I, label=label)

    PyPlot.title("CV curve")
    PyPlot.xlabel(L"\eta \ (V)")
    PyPlot.ylabel(L"I \ (A)")
    PyPlot.legend(loc="best")
    PyPlot.grid()
end

function plot_all(df_list, prms_list)
    PyPlot.figure(figsize=(6,4))
    #PyPlot.figure(figsize=(10,8))
    for i in 1:size(df_list,1)
        #println(i)
        finall_plot(
            df_list[i],
            string("prms = ",prms_list[i])
        )
    end
    
    CV_orig = import_CVtoDataFrame(;T=800,pO2=100)
    #checknodes = get_checknodes(0.01,0.95,-0.95,-0.01,0.05)
    checknodes = get_checknodes(0.01,0.95,-0.95,-0.01,0.1)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_orig, checknodes))
    
    finall_plot(CV_exp, "exp")
end

function plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
    PyPlot.figure(figsize=(6,4))
    
    
    for i1 in 1:size(rprm1,1)
        for i2 in 1:size(rprm2,1)
            if Bool(ok_matrix[i1,i2])
                #println(" ploting i1 / i2 : ",i1, " / ", i2)
                finall_plot(
                    CV_matrix[i1,i2],
                    string("prms = ", rprm1[i1], " / ", rprm2[i2])
                )
            end
        end
    end

    CV_orig = import_CVtoDataFrame(;T=T,pO2=pO2)
    #checknodes = get_checknodes(0.01,0.95,-0.95,-0.01,0.05)
    checknodes = get_checknodes(0.01,0.95,-0.95,-0.01,0.1)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_orig, checknodes))
    
    finall_plot(CV_exp, "exp")
end


#####################
## 2D optimization ##
#####################

mutable struct FF_2D
    x_matrix
    y_matrix
    err_matrix
end



function surfplot_FF_2D(FF::FF_2D; my_title="Objective function")
    surf(FF.x_matrix, FF.y_matrix, FF.err_matrix, cstride=1, #cmap=ColorMap("gray"), 
        alpha=0.8);
    xlabel("prm1")
    ylabel("prm2")
    title(my_title)
    zlabel("error")
end


function get_FF_interp_2D(CV_list, prm1_list, prm2_list; nx=20, ny=20, T=800, pO2=100)
    # ordered for prms A, B as [(A1,B1), (A1,B2), (A2,B1), (A2, B2)]
    # prms_list ordered as [A1 B1 A2 B2]
    CV_raw = import_CVtoDataFrame(;T=T,pO2=pO2)
    checknodes = get_checknodes(0.05,0.95,-0.95,-0.05,0.01)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_raw, checknodes))
    
    CV_with_checknodes_list = []
    
    for i in 1:size(CV_list,1)
        push!(CV_with_checknodes_list, DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[i], checknodes)))
    end
    
    # bilinear interpolation
    rx = range(0.0, stop=1.0, length=nx)
    ry = range(0.0, stop=1.0, length=ny)
    
    x_matrix = zeros(nx,ny)
    y_matrix = zeros(nx,ny)
    err_matrix = zeros(nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            x_matrix[ix, iy] = rx[ix]*(prm1_list[2]-prm1_list[1]).+prm1_list[1]
            y_matrix[ix, iy] = ry[iy]*(prm2_list[2]-prm2_list[1]).+prm2_list[1]
            err_matrix[ix, iy] = fitnessFunction(
                    CV_exp,
                    DataFrame(
                        U = checknodes[:,1], 
                        I = biliComb_CVs(
                            CV_with_checknodes_list[1],
                            CV_with_checknodes_list[2],
                            CV_with_checknodes_list[3],
                            CV_with_checknodes_list[4],
                            rx[ix],ry[iy]).I
                    )
                )
        end
    end
    
    return FF_2D(x_matrix, y_matrix, err_matrix)
end


function all_four_nodes_ok(ok_matrix, x, y)
    for i in 0:1
        for j in 0:1
            if ok_matrix[x+i, y + j] == 0
                return false
            end
        end
    end
    return true
end

function print_ok_matrix(ok_matrix)
println("... printing ok_matrix ...")
    for i2 in 1:size(ok_matrix,2)
        for i1 in 1:size(ok_matrix,1)
        
            print("|", ok_matrix[i1, i2])
        end
        print("|\n")
    end
end

function find_mins(FF, dx; debug_print_bool=false)
    hard_min = Inf
    count_of_mins = 0
    values = []
    positions = []
    curv = []
    
    if size(FF.err_matrix,1) < 3 || size(FF.err_matrix,2) < 3
        println("ERROR: is_min_there: some dimension of FF.err_matrix is smaller than 3")
        return throw(Exception)
    end
    
    for i1 in 2:size(FF.err_matrix,1)-1
        for i2 in 2:size(FF.err_matrix,2)-1
            pot_min = FF.err_matrix[i1, i2]
            if pot_min < hard_min
                hard_min = pot_min
            end
            
            if debug_print_bool
                println(i1, " ", i2, " > ",pot_min)
                
                println(pot_min < FF.err_matrix[i1-1, i2], " ",
                    pot_min < FF.err_matrix[i1+1, i2], " ",
                    pot_min < FF.err_matrix[i1, i2 - 1], " ",
                    pot_min < FF.err_matrix[i1, i2 + 1])
            end
                
            if (pot_min < FF.err_matrix[i1 - 1, i2 - 1] &&
                pot_min < FF.err_matrix[i1 - 1, i2] &&
                pot_min < FF.err_matrix[i1 - 1, i2 + 1] &&
                pot_min < FF.err_matrix[i1, i2 - 1] &&
                pot_min < FF.err_matrix[i1, i2 + 1] &&
                pot_min < FF.err_matrix[i1 + 1, i2 - 1] &&
                pot_min < FF.err_matrix[i1 + 1, i2] &&
                pot_min < FF.err_matrix[i1 + 1, i2 + 1]
                )
                
                c1 = (FF.err_matrix[i1 - 1, i2] - 2*pot_min + FF.err_matrix[i1 + 1, i2])/(dx[1]*dx[1])
                c2 = (FF.err_matrix[i1, i2 - 1] - 2*pot_min + FF.err_matrix[i1, i2 + 1])/(dx[2]*dx[2])
                c12= (FF.err_matrix[i1 - 1, i2 - 1] + FF.err_matrix[i1 + 1, i2 + 1] 
                    - FF.err_matrix[i1 + 1, i2 - 1] - FF.err_matrix[i1 - 1, i2 + 1])/(4*dx[1]*dx[2])
                count_of_mins += 1
                append!(values,pot_min)
                push!(positions, [FF.x_matrix[i1, i2], FF.y_matrix[i1, i2]])
                push!(curv,[c1, c2,c12])
            end
        end
    end
    
    return count_of_mins, values, positions, curv, hard_min
end

function pickup_min(min_count, min_values, min_positions)
    value = Inf;
    position = [0, 0]
    for i in 1:min_count
        if min_values[i] < value
            value = min_values[i]
            position = min_positions[i]
        end
    end
    if value == Inf
        println("ERROR: pickup_min: value = -1")
        throw(Exception)
    else
        return value, position
    end
end

function scan_2D_recursive(;pyplot=false, 
                            T=800, pO2=100,
                            #rprm1=range(15, stop=20, length=3),
                            #rprm2=range(15, stop=20, length=3),
                            #rprm1=range(18.4375, stop=18.59375, length=3),
                            #rprm2=range(20.59375, stop=21, length=3),
                            rprm1=range(13.0, stop=25., length=8),
                            #rprm2=range(15.0, stop=23, length=7),
                            rprm2=range(-1.5, stop=1.5, length=8),
                            nx=200, ny=200,
                            wp=[1, 1, 0, 0, 0, 0], #TODO !!!
                            depth_remaining=5,
                            approx_min=Inf,
                            approx_min_position=[0, 0],
                            recursive=false,
                            recursive_string="",
                            )

    println("nx = ",nx, " ..... ny = ",ny)
    #PyPlot.close()
    #PyPlot.close()
    
    ######################################################
    #rprm1 = collect(-5.0 : 1.0 : 5.0)
    #rprm1 = collect(-3. : 0.5 : 3.0)
    #rprm1=range(15, stop=20, length=3)
    
    #rprm2 = collect(-3. : 0.5 : 3.0)
    #rprm2=range(15, stop=20, length=3)
    ######################################################
    
    wpn = ["A0","R0","DGA","DGR","beta","A"]
    #A0 = 21.71975544711280
    #R0 = 19.53
    A0 = 19.50
    R0 = 19.85
    DGA = 0.0905748
    DGR = -0.708014
    beta = 0.6074566741435283
    A = -0.356

    println(recursive_string,"computing 2D scan... ")
    lrprm1 = size(rprm1,1)
    lrprm2 = size(rprm2,1)
    
    global_min = Inf
    global_min_position = [0,0]
    
    CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
    ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
    println(recursive_string,"Computing CVs")
    @time(
        for i1 in 1:lrprm1
            for i2 in 1:lrprm2
                try
                    print(recursive_string,rprm1[i1], " / ", rprm2[i2])
                    
                    ##########################################
                    R0 = rprm1[i1]
                    A = rprm2[i2]
                    ##########################################
                    
                    CV_matrix[i1,i2] = ysz_voltammetry_like.run_new(
                        out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
                        prms_in=[A0, R0, DGA, DGR, beta, A]
                        )
                    ok_matrix[i1, i2] = 1
                    println("  .. ok :)")
                catch e
                    if e isa InterruptException
                        rethrow(e)
                    else
                        ok_matrix[i1, i2] = 0
                        println("  :(   ")
                        continue
                    end
                end
            end
        end
    )
    
    #PyPlot.figure(figsize=(6,4))
    #xlabel("parameter")
    #zlabel("error")
    #title("Fitness Function")

    if pyplot && !(recursive)
        figure(figsize=(3,3))
    end
    #print_ok_matrix(ok_matrix)
    
    println(recursive_string,"Computing interpolation ... ")
    @time(
        for i1 in 1:lrprm1-1
            for i2 in 1:lrprm2-1
                if all_four_nodes_ok(ok_matrix,i1,i2)
                    print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
                    
                    FF=get_FF_interp_2D(
                        [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
                        rprm1[i1:i1+1],
                        rprm2[i2:i2+1];
                        nx=nx, ny=ny,
                        T=T, pO2=pO2
                    )
                    
                    min_count, min_values, min_positions, min_curv, hard_min = find_mins(
                        FF, 
                        [(rprm1[i1+1] - rprm1[i1])/nx,
                         (rprm1[i2+1] - rprm1[i2])/ny]
                    )
                    print("... h_min = ",hard_min)
                    
                    
                    if min_count > 0
                        println(" ... min_c ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
                        #println(" ... min ", min_count)#
                        
                        approx_min, approx_min_position = pickup_min(min_count, min_values, min_positions)
                        
                        if depth_remaining > 0
                            println(recursive_string,">> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
                            new_min, new_min_position = scan_2D_recursive(;pyplot=true, 
                                T=T, pO2=pO2,
                                rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
                                rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
                                approx_min=approx_min,
                                approx_min_position=approx_min_position,
                                recursive=true,
                                depth_remaining=depth_remaining - 1,
                                recursive_string=string(recursive_string, "<", depth_remaining - 1, ">"),
                                nx = Int32(round(nx/1.2)), ny = Int32(round(ny/1.2))
                            )
                            if new_min < global_min
                                global_min = new_min
                                global_min_position =new_min_position
                            end
                        else
                            if approx_min < global_min
                                global_min = approx_min
                                global_min_position = approx_min_position
                            end
                            if pyplot
                                surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
                            end
                        end
                    else
                        if pyplot
                            surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
                        end
                        println()
                    end
                end 
            end
        end
    )
    
    
    
    if pyplot && !(recursive)
        plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
    end
    return global_min, global_min_position
end



# function scan_2D(;pyplot=false, 
#                 T=800, pO2=100,
# #                rprm1=range(15, stop=20, length=3),
# #                rprm2=range(15, stop=20, length=3),
#                 #rprm1=range(18.4375, stop=18.59375, length=3),
#                 #rprm2=range(20.59375, stop=21, length=3),
#                 rprm1=range(18.4375, stop=18.59375, length=3),
#                 rprm2=range(20.59375, stop=21, length=3),
#                 wp=[1, 1, 0, 0, 0, 0], #TODO !!!
#                 depth_remaining=2,
#                 recursive=false,
#                 nx=100, ny=100
#                 )
# 
#     
#     #PyPlot.close()
#     #PyPlot.close()
#     
#     ######################################################
#     #rprm1 = collect(-5.0 : 1.0 : 5.0)
#     #rprm1 = collect(-3. : 0.5 : 3.0)
#     #rprm1=range(15, stop=20, length=3)
#     
#     #rprm2 = collect(-3. : 0.5 : 3.0)
#     #rprm2=range(15, stop=20, length=3)
#     ######################################################
#     
#     wpn = ["A0","R0","DGA","DGR","beta","A"]
#     A0 = 21.71975544711280
#     R0 = 19.53
#     DGA = 0.0905748
#     DGR = -0.708014
#     beta = 0.6074566741435283
#     A = -0.356
# 
#     println("computing 2D scan... ")
#     lrprm1 = size(rprm1,1)
#     lrprm2 = size(rprm2,1)
#     
#     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
#     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
#     println("Computing CVs")
#     @time(
#         for i1 in 1:lrprm1
#             for i2 in 1:lrprm2
#                 try
#                     print(rprm1[i1], " / ", rprm2[i2])
#                     
#                     ##########################################
#                     A0 = rprm1[i1]
#                     R0 = rprm2[i2]
#                     ##########################################
#                     
#                     CV_matrix[i1,i2] = ysz_voltammetry_like.run_new(
#                         out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
#                         prms_in=[A0, R0, DGA, DGR, beta, A]
#                         )
#                     ok_matrix[i1, i2] = 1
#                     println("  .. ok :)")
#                 catch e
#                     if e isa InterruptException
#                         rethrow(e)
#                     else
#                         ok_matrix[i1, i2] = 0
#                         println("  :(   ")
#                         continue
#                     end
#                 end
#             end
#         end
#     )
#     
#     #PyPlot.figure(figsize=(6,4))
#     #xlabel("parameter")
#     #zlabel("error")
#     #title("Fitness Function")
# 
#     if pyplot && !(recursive)
#         figure(figsize=(3,3))
#     end
#     print_ok_matrix(ok_matrix)
#     
#     println("Computing interpolation ... ")
#     @time(
#         for i1 in 1:lrprm1-1
#             for i2 in 1:lrprm2-1
#                 if all_four_nodes_ok(ok_matrix,i1,i2)
#                     print(" intrp i1 / i2 : ",i1, " / ", i2)
#                     
#                     FF=get_FF_interp_2D(
#                         [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
#                         rprm1[i1:i1+1],
#                         rprm2[i2:i2+1];
#                         nx=nx, ny=ny,
#                         T=T, pO2=pO2
#                     )
#                     
#                     min_count, min_values, min_positions, min_curv = find_mins(
#                         FF.err_matrix, 
#                         [(rprm1[i1+1] - rprm1[i1])/nx,
#                          (rprm1[i2+1] - rprm1[i2])/ny]
#                     )
#                     
# 
#                     
#                     if pyplot
#                         println("jsem tu")
#                         surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
#                     end
#                     
#                     if min_count > 0
#                         print(" ... min ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
#                         
#                         if depth_remaining > 0
#                             println("\n >>>>>>>> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
#                             scan_2D(;pyplot=true, 
#                                 T=T, pO2=pO2,
#                                 rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
#                                 rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
#                                 depth_remaining=depth_remaining - 1,
#                                 recursive=true
#                                 )
#                         end
#                     end
#                     
#                     print("\n")
#                 end 
#             end
#         end
#     )
#     
#     
#     
#     if pyplot && !(recursive)
#         plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
#     end
#     return
# end



function del_pyplot()
    for i in 1:10
        PyPlot.close()
    end
end










#####################
## 1D optimization ##
#####################
function get_FF_interp_1D(CV_list, prms_list; n=20)
    CV_raw = import_CVtoDataFrame(;T=800,pO2=100)
    checknodes = get_checknodes(0.05,0.95,-0.95,-0.05,0.01)
    CV_exp = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_raw, checknodes))

    CV1 = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[1],checknodes))
    CV2 = DataFrame(U = checknodes[:,1], I = CV_get_I_values(CV_list[end],checknodes))
    node_list = range(0.0, stop=1.0, length=n)
    err_list_interp = [fitnessFunction(
            CV_exp,
            DataFrame(
                U = checknodes[:,1], 
                I = (linComb_CVs(CV1, CV2, i)).I
            )
        )
        for i in node_list
    ]
    prms_list_interp = node_list.*(prms_list[end]-prms_list[1]).+prms_list[1]
    FF_out = DataFrame(prms = prms_list_interp, err = err_list_interp)
end

function comparison()
    CV_list = []
    prms_list = []
    PyPlot.close()
    
    #for A in collect(-0.36: 0.1 : -0.32)
    #for A in [-0.345, -0.33]
    for A in [-0.356, -0.354]
    #for A in [-0.1, -0.6]
    #for A in collect(-0.1 : -0.01 : -0.6)
    #for A in [-0.3893]
    #for R0 in [21.4, 21.8]
    #for R0 in collect(21.4: 0.02 : 21.8)
    for R0 in [20.53]
        println(A)
        ww=DataFrame()
        push!(CV_list,ysz_voltammetry_like.run_new(
            out_df_bool=true, voltammetry=true, sample=10,
            prms_in=[21.71975544711280, R0, 0.0905748, -0.708014, 0.6074566741435283, A]
            )
        )
        push!(prms_list,A)
    end
    end
    println("jdeme dal")
    plot_all(CV_list, prms_list)
    plot_error(CV_list, prms_list)
    return
end


function scan_1D()
    PyPlot.close()
    PyPlot.close()
    #rA = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.2]    
    #rA = collect(-0.4 : 0.02 : -0.3)
    
    #rprms = [19, 20, 21, 22]
    rprms = collect(-0.5 : 0.1 : 0.0)
    
    
    A0 = 21.71975544711280
    R0 = 20.606423236896422
    DGA = 0.0905748
    DGR = -0.708014
    beta = 0.6074566741435283
    A = -0.356
    
    println("computing ... ")
    CV_list = []
    ok_prms_list = []
    for i in 1:size(rprms,1)
        try
        print(rprms[i])
        
        ########
        A = rprms[i]
        ########
        
        push!(CV_list,ysz_voltammetry_like.run_new(
            out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
            prms_in=[A0, R0, DGA, DGR, beta, A]
            )
        )
        push!(ok_prms_list,rprms[i])
        println("  .. ok")
        catch e
            if e isa InterruptException
                rethrow(e)
            else
                println()
                continue
            end
        end
    end

    PyPlot.figure(figsize=(6,4))
    xlabel("parameter")
    ylabel("error")
    title("Fitness Function")
    FF = DataFrame()
    for i in 1:(size(ok_prms_list,1)-1)
        FF = get_FF_interp_1D([CV_list[i], CV_list[i+1]],[ok_prms_list[i], ok_prms_list[i+1]])
        plot(FF.prms, FF.err, linewidth=3.0, label=string("FF ",i))
        plot(FF.prms[1],FF.err[1], "kx", markersize=12.0, ) 
        #plot(FF.prms, FF.err)
    end
    plot(FF.prms[end],FF.err[end], "kx", markersize=12.0)
    legend(loc="best")
    grid()
    
    plot_all(CV_list, ok_prms_list)
    
    return
end


end

