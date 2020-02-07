module ysz_fitting
#######################
# TODO 
# [  ] general global search using projection to each variable
# [  ] compute EIS exacly on checknodes and therefore remove plenty of "EIS_apply_checknodes"
# [  ] put appropriate and finished stuff into "CV_fitting_supporting_stuff"
# [o] implement LM algorithm
# [  ] sjednotit, co znamena pO2, jeslti jsou to procenta a jak se prenasi do simulace!  a taky T .. Celsia a Kelvina !!!
# [x] vymyslet novy relevantni vektor parametru
# [  ] aplikovat masku pro fitting pro scan_2D_recursive
#######################



using PyPlot
using DataFrames
using LeastSquaresOptim

include("../examples/ysz_experiments.jl")
include("../src/CV_fitting_supporting_stuff.jl")
include("../src/EIS_fitting_supporting_stuff.jl")
include("../src/import_experimental_data.jl")


###########
###########




function EIS_get_universal_checknodes()
    return EIS_get_checknodes_geometrical(1, 7500, 1.2)
end

function get_shared_prms() 
    #     old_prms =  [21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1]
    # first dump fit = (A0, R0, DGA, DGR, beta, A) = (18.8, 19.2,     -0.14,      -0.5,   0.6074566741435283,   0.956)
    return (A0, R0, DGA, DGR, beta, A) = (18.8, 19.2,     -0.14,      -0.5,   0.6074566741435283,   0.956)
end

function get_shared_add_prms()
    # cap_fit + resistance fit = (DD, nu, nus, ms_par) = (3.05e-13,    0.85,      0.21,   0.05)
    return (DD, nu, nus, ms_par) = (3.05e-13,    0.85,      0.21,   0.05)
end








##########
##########





function EIS_apply_checknodes(EIS_in,checknodes)
    DataFrame( f = checknodes, Z = EIS_get_Z_values(EIS_in, checknodes))
end


function CV_simple_run()
    old_prms = [21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1]
   @time  ysz_experiments.run_new(
        out_df_bool=true, voltammetry=true, sample=40, pyplot_finall=true,
        prms_in=[19.50, 19.00, 0.0905748, -0.708014, 0.6074566741435283, -1.5]
    )
    return
end

function EIS_simple_run()    
    
    pO2_in_sim = 1.0
    EIS_bias=0.0
    
    (A0, R0, DGA, DGR, beta, A) = get_shared_prms()
    A0 = 18
    DGR=0.1
    DGA=0.1
    
    (DD, nu, nus, ms_par) = get_shared_add_prms()
    DD =  1.0e-12
    
    EIS_df = ysz_experiments.run_new(
        pyplot=false, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias,
        dx_exp=-11,
        # TODO !!! T a pO2
        prms_in=[A0, R0, DGA, DGR, beta, A],
        add_prms_in=(DD, nu, nus, ms_par)
    )
    @show EIS_df
    checknodes =  EIS_get_universal_checknodes()
    Nyquist_plot(EIS_apply_checknodes(EIS_df,checknodes), "s_r")
end

# # # function get_dfs_to_interpolate(
# # #                             EIS_bool=true, pyplot=false,
# # #                            T=800, pO2=100, EIS_bias=0.0,
# # #                             #rprm1=range(15, stop=20, length=3),
# # #                             #rprm2=range(15, stop=20, length=3),
# # #                             #rprm1=range(18.4375, stop=18.59375, length=3),
# # #                             #rprm2=range(20.59375, stop=21, length=3),
# # #                             rprm1=range(19.3, stop=21.0, length=2),
# # #                             #rprm2=range(15.0, stop=23, length=7),
# # #                             rprm2=range(-3.0, stop=-1.6, length=2),
# # #                             nx=100, ny=100,
# # #                             wp=[1, 1, 0, 0, 0, 0], #TODO !!!
# # #                             depth_remaining=1,
# # #                             approx_min=Inf,
# # #                             approx_min_position=[0, 0],
# # #                             recursive=false,
# # #                             recursive_string="",
# # #                             )
# # # 
# # #     println("nx = ",nx, " ..... ny = ",ny)
# # #     #PyPlot.close()
# # #     #PyPlot.close()
# # #     
# # #     ######################################################
# # #     #rprm1 = collect(-5.0 : 1.0 : 5.0)
# # #     #rprm1 = collect(-3. : 0.5 : 3.0)
# # #     #rprm1=range(15, stop=20, length=3)
# # #     
# # #     #rprm2 = collect(-3. : 0.5 : 3.0)
# # #     #rprm2=range(15, stop=20, length=3)
# # #     ######################################################
# # #     
# # #     wpn = ["A0","R0","DGA","DGR","beta","A"]
# # #     #A0 = 21.71975544711280
# # #     #R0 = 19.53
# # #     A0 = 19.50
# # #     R0 = 19.85
# # #     DGA = 0.0905748
# # #     DGR = -0.708014
# # #     beta = 0.6074566741435283
# # #     A = -0.356
# # #     
# # #     (A0, R0, DGA, DGR, beta, A) = (21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1)
# # # 
# # #     println(recursive_string,"computing 2D scan... ")
# # #     lrprm1 = size(rprm1,1)
# # #     lrprm2 = size(rprm2,1)
# # #     
# # #     global_min = Inf
# # #     global_min_position = [0,0]
# # #     
# # #     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
# # #     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
# # #     if EIS_bool
# # #         println(recursive_string,"Computing EISs")
# # #     else
# # #         println(recursive_string,"Computing CVs")
# # #     end
# # #     @time(
# # #         for i1 in 1:lrprm1
# # #             for i2 in 1:lrprm2
# # #                 try
# # #                     print(recursive_string,rprm1[i1], " / ", rprm2[i2])
# # #                     
# # #                     ##########################################
# # #                     R0 = rprm1[i1]
# # #                     A = rprm2[i2]
# # #                     ##########################################
# # #                     
# # #                    if EIS_bool 
# # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # #                             pyplot=false, EIS_IS=true, out_df_bool=true, 
# # #                             tref=0,
# # #                             prms_in=[A0, R0, DGA, DGR, beta, A],  
# # #                             nu_in=0.9, pO2_in=1.0, DD_in=9.5658146540360312e-10
# # #                         )
# # #                     else
# # #                         
# # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # #                             out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
# # #                             prms_in=[A0, R0, DGA, DGR, beta, A]
# # #                         )
# # #                     end
# # #                     ok_matrix[i1, i2] = 1
# # #                     println("  .. ok :)")
# # #                 catch e
# # #                     if e isa InterruptException
# # #                         rethrow(e)
# # #                     else
# # #                         ok_matrix[i1, i2] = 0
# # #                         println("  :(   ")
# # #                         continue
# # #                     end
# # #                 end
# # #             end
# # #         end
# # #     )
# # #     
# # #     #PyPlot.figure(figsize=(6,4))
# # #     #xlabel("parameter")
# # #     #zlabel("error")
# # #     #title("Fitness Function")
# # # 
# # #     if pyplot && !(recursive)
# # #         figure(figsize=(3,3))
# # #     end
# # #     #print_ok_matrix(ok_matrix)
# # #     
# # #     println(recursive_string,"Preparing for interpolation ... ")
# # #         for i1 in 1:lrprm1-1
# # #             for i2 in 1:lrprm2-1
# # #                 if all_four_nodes_ok(ok_matrix,i1,i2)
# # #                     print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
# # #                     
# # #                     if EIS_bool
# # #                             EIS_list = [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]]
# # #                             EIS_with_checknodes_list = []
# # #                             checknodes =  EIS_get_universal_checknodes()
# # #                             for i in 1:size(EIS_list,1)
# # #                                 push!(EIS_with_checknodes_list, DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_list[i], checknodes)))
# # #                             end
# # #                             return EIS_with_checknodes_list
# # #                     else
# # #                         # TODO !!!
# # #                         FF=CV_get_FF_interp_2D(
# # #                             [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
# # #                             rprm1[i1:i1+1],
# # #                             rprm2[i2:i2+1];
# # #                             nx=nx, ny=ny,
# # #                             T=T, pO2=pO2
# # #                         )
# # #                     end
# # #                 end 
# # #             end
# # #         end
# # # end

function EIS_view_exp(T=800, pO2=100, EIS_bias=0.0,)
    EIS_raw = import_EIStoDataFrame(;T=T,pO2=pO2,bias=EIS_bias)
    checknodes =  EIS_get_universal_checknodes()
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    Nyquist_plot(EIS_exp,"exp")
end

function EIS_view_interpolation_at(Q_list, x, y)
    EIS_intrp = DataFrame(
        f = Q_list[1].f,
        Z = EIS_biliComb(
            Q_list[1],
            Q_list[2],
            Q_list[3],
            Q_list[4],
            x,y).Z
    )
    Nyquist_plot(EIS_intrp, string("x/y = ",x,"/",y))
end

function finall_plot(df, label="")
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



function plot_EIS_matrix(EIS_matrix, ok_matrix, rprm1, rprm2, T, pO2, EIS_bias=0.0)
    PyPlot.figure(figsize=(6,4))
    
    checknodes =  EIS_get_universal_checknodes()  
    
    for i1 in 1:size(rprm1,1)
        for i2 in 1:size(rprm2,1)
            if Bool(ok_matrix[i1,i2])
                #println(" ploting i1 / i2 : ",i1, " / ", i2)
                Nyquist_plot(
                    EIS_apply_checknodes(EIS_matrix[i1, i2], checknodes),
                    string("prms = ", rprm1[i1], " / ", rprm2[i2])
                )
            end
        end
    end

    EIS_raw = import_EIStoDataFrame(;T=T,pO2=pO2,bias=EIS_bias)
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    Nyquist_plot(EIS_exp, "exp")
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


function CV_get_FF_interp_2D(CV_list, prm1_list, prm2_list; nx=20, ny=20, T=800, pO2=100)
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



function EIS_get_FF_interp_2D(EIS_list, prm1_list, prm2_list; nx=20, ny=20, T=800, pO2=100, EIS_bias=0.0)
    # ordered for prms A, B as [(A1,B1), (A1,B2), (A2,B1), (A2, B2)]
    # prms_list ordered as [A1 B1 A2 B2]
    EIS_raw = import_EIStoDataFrame(;T=T,pO2=pO2,bias=EIS_bias)
    checknodes =  EIS_get_universal_checknodes()
    EIS_exp = DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_raw, checknodes))
    
    EIS_with_checknodes_list = []
    
    for i in 1:size(EIS_list,1)
        push!(EIS_with_checknodes_list, DataFrame(f = checknodes[:], Z = EIS_get_Z_values(EIS_list[i], checknodes)))
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
            err_matrix[ix, iy] = EIS_fitnessFunction(
                    EIS_exp,
                    DataFrame(
                        f = checknodes[:], 
                        Z = EIS_biliComb(
                            EIS_with_checknodes_list[1],
                            EIS_with_checknodes_list[2],
                            EIS_with_checknodes_list[3],
                            EIS_with_checknodes_list[4],
                            rx[ix],ry[iy]).Z
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
                            EIS_bool=true,
                            T=800, pO2=100, EIS_bias=0.0,
                            
                            wp=[1, 1, 0, 0, 0, 0], #TODO !!!
                            rprm1=range(17, stop=18, length=2),
                           #rprm2=range(15.0, stop=23, length=7),
                            rprm2=range(18, stop=19, length=2),
                            
                            nx=30, ny=30,

                            depth_remaining=1,
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

    (A0, R0, DGA, DGR, beta, A) = get_shared_prms()
   (DD, nu, nus, ms_par)  =  get_shared_add_prms()
    
   
    println(recursive_string,"computing 2D scan... ")
    lrprm1 = size(rprm1,1)
    lrprm2 = size(rprm2,1)
    
    global_min = Inf
    global_min_position = [0,0]
    
    CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
    ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
    if EIS_bool
        println(recursive_string,"Computing EISs")
    else
        println(recursive_string,"Computing CVs")
    end
    @time(
        for i1 in 1:lrprm1
            for i2 in 1:lrprm2
                try
                    print(recursive_string,rprm1[i1], " / ", rprm2[i2])
                    
                    ##########################################
                    A0 = rprm1[i1]
                    R0 = rprm2[i2]
                    ##########################################
                    
                   if EIS_bool
                        CV_matrix[i1,i2] = ysz_experiments.run_new(
                            pyplot=false, EIS_IS=true, out_df_bool=true, EIS_bias=EIS_bias,
                            dx_exp=-11,
                            # TODO !!! T a pO2
                            prms_in=[A0, R0, DGA, DGR, beta, A],
                            add_prms_in=(DD, nu, nus, ms_par)
                        )
                    else
                        CV_matrix[i1,i2] = ysz_experiments.run_new(
                            out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
                            prms_in=[A0, R0, DGA, DGR, beta, A],
                            add_prms_in=(DD, nu, nus, ms_par)
                        )
                    end
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
                    
                    if EIS_bool
                         FF=EIS_get_FF_interp_2D(
                            [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
                            rprm1[i1:i1+1],
                            rprm2[i2:i2+1];
                            nx=nx, ny=ny,
                            T=T, pO2=pO2, EIS_bias=EIS_bias
                        )                       
                    else
                        FF=CV_get_FF_interp_2D(
                            [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
                            rprm1[i1:i1+1],
                            rprm2[i2:i2+1];
                            nx=nx, ny=ny,
                            T=T, pO2=pO2
                        )
                    end
                    
                    min_count, min_values, min_positions, min_curv, hard_min = find_mins(
                        FF, 
                        [(rprm1[i1+1] - rprm1[i1])/nx,
                         (rprm2[i2+1] - rprm2[i2])/ny]
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
        if EIS_bool 
            plot_EIS_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2, EIS_bias)
        else
            plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
        end
    end
    return global_min, global_min_position
end


# # # # function scan_2D_recursive(;pyplot=false, 
# # # #                             T=800, pO2=100,
# # # #                             #rprm1=range(15, stop=20, length=3),
# # # #                             #rprm2=range(15, stop=20, length=3),
# # # #                             #rprm1=range(18.4375, stop=18.59375, length=3),
# # # #                             #rprm2=range(20.59375, stop=21, length=3),
# # # #                             rprm1=range(10.0, stop=28., length=4),
# # # #                             #rprm2=range(15.0, stop=23, length=7),
# # # #                             rprm2=range(-2.5, stop=1.0, length=4),
# # # #                             nx=100, ny=100,
# # # #                             wp=[1, 1, 0, 0, 0, 0], #TODO !!!
# # # #                             depth_remaining=1,
# # # #                             approx_min=Inf,
# # # #                             approx_min_position=[0, 0],
# # # #                             recursive=false,
# # # #                             recursive_string="",
# # # #                             )
# # # # 
# # # #     println("nx = ",nx, " ..... ny = ",ny)
# # # #     #PyPlot.close()
# # # #     #PyPlot.close()
# # # #     
# # # #     ######################################################
# # # #     #rprm1 = collect(-5.0 : 1.0 : 5.0)
# # # #     #rprm1 = collect(-3. : 0.5 : 3.0)
# # # #     #rprm1=range(15, stop=20, length=3)
# # # #     
# # # #     #rprm2 = collect(-3. : 0.5 : 3.0)
# # # #     #rprm2=range(15, stop=20, length=3)
# # # #     ######################################################
# # # #     
# # # #     wpn = ["A0","R0","DGA","DGR","beta","A"]
# # # #     #A0 = 21.71975544711280
# # # #     #R0 = 19.53
# # # #     A0 = 19.50
# # # #     R0 = 19.85
# # # #     DGA = 0.0905748
# # # #     DGR = -0.708014
# # # #     beta = 0.6074566741435283
# # # #     A = -0.356
# # # # 
# # # #     println(recursive_string,"computing 2D scan... ")
# # # #     lrprm1 = size(rprm1,1)
# # # #     lrprm2 = size(rprm2,1)
# # # #     
# # # #     global_min = Inf
# # # #     global_min_position = [0,0]
# # # #     
# # # #     CV_matrix = Array{DataFrame}(undef,lrprm1, lrprm2)
# # # #     ok_matrix = Array{Int32}(undef,lrprm1, lrprm2)
# # # #     println(recursive_string,"Computing CVs")
# # # #     @time(
# # # #         for i1 in 1:lrprm1
# # # #             for i2 in 1:lrprm2
# # # #                 try
# # # #                     print(recursive_string,rprm1[i1], " / ", rprm2[i2])
# # # #                     
# # # #                     ##########################################
# # # #                     R0 = rprm1[i1]
# # # #                     A = rprm2[i2]
# # # #                     ##########################################
# # # #                     
# # # #                    if true
# # # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # # #                             out_df_bool=true, voltammetry=true, sample=10, #pyplot=true,
# # # #                             prms_in=[A0, R0, DGA, DGR, beta, A]
# # # #                         )
# # # #                     else
# # # #                         CV_matrix[i1,i2] = ysz_experiments.run_new(
# # # #                             pyplot=false, EIS_IS=true, out_df_bool=false, 
# # # #                             tref=0,
# # # #                             prms_in=[A0, R0, DGA, DGR, beta, A],  
# # # #                             nu_in=0.9, pO2_in=1.0, DD_in=9.5658146540360312e-10
# # # #                         )
# # # #                     end
# # # #                     ok_matrix[i1, i2] = 1
# # # #                     println("  .. ok :)")
# # # #                 catch e
# # # #                     if e isa InterruptException
# # # #                         rethrow(e)
# # # #                     else
# # # #                         ok_matrix[i1, i2] = 0
# # # #                         println("  :(   ")
# # # #                         continue
# # # #                     end
# # # #                 end
# # # #             end
# # # #         end
# # # #     )
# # # #     
# # # #     #PyPlot.figure(figsize=(6,4))
# # # #     #xlabel("parameter")
# # # #     #zlabel("error")
# # # #     #title("Fitness Function")
# # # # 
# # # #     if pyplot && !(recursive)
# # # #         figure(figsize=(3,3))
# # # #     end
# # # #     #print_ok_matrix(ok_matrix)
# # # #     
# # # #     println(recursive_string,"Computing interpolation ... ")
# # # #     @time(
# # # #         for i1 in 1:lrprm1-1
# # # #             for i2 in 1:lrprm2-1
# # # #                 if all_four_nodes_ok(ok_matrix,i1,i2)
# # # #                     print(recursive_string," intrp i1 / i2 : ",i1, " / ", i2, " for prms: ", rprm1[i1], " / ", rprm2[i2])
# # # #                     
# # # #                     FF=get_FF_interp_2D(
# # # #                         [CV_matrix[i1,i2], CV_matrix[i1,i2+1], CV_matrix[i1+1,i2], CV_matrix[i1+1,i2+1]], 
# # # #                         rprm1[i1:i1+1],
# # # #                         rprm2[i2:i2+1];
# # # #                         nx=nx, ny=ny,
# # # #                         T=T, pO2=pO2
# # # #                     )
# # # #                     
# # # #                     min_count, min_values, min_positions, min_curv, hard_min = find_mins(
# # # #                         FF, 
# # # #                         [(rprm1[i1+1] - rprm1[i1])/nx,
# # # #                          (rprm1[i2+1] - rprm1[i2])/ny]
# # # #                     )
# # # #                     print("... h_min = ",hard_min)
# # # #                     
# # # #                     
# # # #                     if min_count > 0
# # # #                         println(" ... min_c ", min_count, "\n    val ", min_values, "\n    pos ", min_positions, "\n    curv ", min_curv)
# # # #                         #println(" ... min ", min_count)#
# # # #                         
# # # #                         approx_min, approx_min_position = pickup_min(min_count, min_values, min_positions)
# # # #                         
# # # #                         if depth_remaining > 0
# # # #                             println(recursive_string,">> Going deeper in intrp i1 / i2 : ",i1, " / ", i2)
# # # #                             new_min, new_min_position = scan_2D_recursive(;pyplot=true, 
# # # #                                 T=T, pO2=pO2,
# # # #                                 rprm1=range(rprm1[i1], stop=rprm1[i1+1], length=3),
# # # #                                 rprm2=range(rprm2[i2], stop=rprm2[i2+1], length=3),
# # # #                                 approx_min=approx_min,
# # # #                                 approx_min_position=approx_min_position,
# # # #                                 recursive=true,
# # # #                                 depth_remaining=depth_remaining - 1,
# # # #                                 recursive_string=string(recursive_string, "<", depth_remaining - 1, ">"),
# # # #                                 nx = Int32(round(nx/1.2)), ny = Int32(round(ny/1.2))
# # # #                             )
# # # #                             if new_min < global_min
# # # #                                 global_min = new_min
# # # #                                 global_min_position =new_min_position
# # # #                             end
# # # #                         else
# # # #                             if approx_min < global_min
# # # #                                 global_min = approx_min
# # # #                                 global_min_position = approx_min_position
# # # #                             end
# # # #                             if pyplot
# # # #                                 surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
# # # #                             end
# # # #                         end
# # # #                     else
# # # #                         if pyplot
# # # #                             surfplot_FF_2D(FF, my_title=string("Depth_remaining = ", depth_remaining))
# # # #                         end
# # # #                         println()
# # # #                     end
# # # #                 end 
# # # #             end
# # # #         end
# # # #     )
# # # #     
# # # #     
# # # #     
# # # #     if pyplot && !(recursive)
# # # #         plot_CV_matrix(CV_matrix, ok_matrix, rprm1, rprm2, T, pO2)
# # # #     end
# # # #     return global_min, global_min_position
# # # # end



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
#                     CV_matrix[i1,i2] = ysz_experiments.run_new(
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
        push!(CV_list,ysz_experiments.run_new(
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
        
        push!(CV_list,ysz_experiments.run_new(
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



####################################
####################################
####################################
######## Levenberg-Maquard ##########
####################################

function prepare_prms(mask, x0, x)
    prms = zeros(0)
    xi = 1
    for i in collect(1 : 1 :length(mask))
        if convert(Bool,mask[i])
            append!(prms, x[xi])
            xi += 1
        else
            append!(prms, x0[i])
        end
    end
    return prms
end

function LM_optimize()
    function to_optimize(x) 
        #err = run_new(print_bool=false, fitting=true, voltametry=true, pyplot=false, voltrate=0.005, sample=8, bound=0.41, 
        #    prms_in=x)
        prms = prepare_prms(mask, x0, x)
        print(" >> mask = ",mask)
        print(" || prms = ",prms)
        
        #err = run_new(print_bool=false, fitting=true, voltametry=true, pyplot=true, voltrate=0.005, sample=50, upp_bound=0.474, low_bound=-0.429, 
        #    prms_in=prms)
            
        err = run_new(print_bool=true, fitting=true, voltametry=true, dlcap=false, pyplot=true, voltrate=0.001, sample=10, upp_bound=0.55, low_bound=-0.548, 
            prms_in=prms, width=0.45e-3)
        
        println(" || err =", err)
        return [err]
    end
    #x0 = zeros(2)
    #optimize(rosenbrock, x0, LevenbergMarquardt())
    
    lower_bounds = [-20, -20, -1.2, -1.2, 0.1, 0.1]
    upper_bounds = [15, 22, 0.2, 0.2, 0.9, 1.0]
    
    #x0 = [3.1149, -9.59917, 0.498534, -0.597552, 0.724949, 1.17139]
    #x0 = [3.40882, -4.1195, 0.0, 0.0, 0.5, 1.17139]
    #x0 = [4.42991, 19.03254, 0.0, 0.0, 0.5, 0.301]
    #x0 = [4.42991, 20.03254, 0.0, 0.0, 0.5, 0.151]
    
    #x0 = [9.00137, 20.1747, 0.5, 0.0, 0.5, 0.206002]
    x0 = [2.3, 19.5, 1, -1, 0.6, 1.2]
    #x0 = [2.32407, 18.4458, 1.0, -1.0, 0.769217, 1.03349]
    #x0 = [2.34223, 18.4039, 1.0, -1.0, 0.72743, 0.95021] # quite good fit <<<<<<<<<<
    #x0 = [2.34223, 20.8039, 1.0, -1.0, 0.52743, 0.25021] # hand fit <<< usable
    #x0 = [2.74223, 19.7039, 1.0, -1.0, 0.58, 0.25021]
    #x0 = [2.44223, 20.5039, 1.0, -1.0, 0.61, 0.25021] # not bad :)) <<<<<<<<<
    
    # err metric <- check_nodes_long()
    x0 = [2.54223, 20.5039, 1.0000000, -1.000000, 0.500000, 0.250210]  # by hand
    x0 = [2.54223, 20.5039, 1.0000000, -1.000000, 0.632507, 0.250210] # fitted beta << GOOD << err =0.006814132871406775
    x0 = [2.54223, 20.5039, 1.0000000, -1.000000, 0.632507, 0.253287] # fitted A  << err =0.006807257238292433
    x0 = [2.54223, 20.5033, 1.0000000, -1.000000, 0.632507, 0.253287] # fitted R0 << err =0.0068072339219288356
    x0 = [2.55263, 20.5033, 1.0000000, -1.000000, 0.632507, 0.253287] # fitted A0 << err =0.006745444782481952
    x0 = [2.58082, 20.4100, 1.0000000, -1.000000, 0.632507, 0.253287] # fitted A0, R0 << err =0.006574359568544145
    x0 = [2.58082, 20.4100, 0.0905748, -0.708014, 0.632507, 0.253287] # fitted DGA, DGR << err =0.006573859072513787
    x0 = [2.70077, 20.5264, 0.0905748, -0.708014, 0.598817, 0.143269] # fittet 110011 << err =0.005685020802399199
    
    # err metric <- check_nodes_whole()
    x0 = [2.74851, 20.5631, 0.0905748, -0.708014, 0.605159, 0.105409] # fitting... 110011 << err =0.0054800474768966585
    x0 = [2.73650, 20.6063, 0.0905748, -0.708014, 0.607443, 0.100000] # fitting... 110011 << err =0.005415825589421705
    x0 = [2.73645, 20.6064, 0.0905748, -0.708014, 0.607457, 0.100000] # fitted 110011 <<  err =0.0054158249335496105
    x0 = [2.736451985137371, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1] # ACURATE FINALL FIT
    # >> x0corr = [21.71975544711280, 20.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1] # ACURATE FINALL FIT
    #x0 = [6.736451985137371, 24.606423236896422, 0.0905748, -0.708014, 0.6074566741435283, 0.1]
    #x0 = [2.736451985137371, 20.606423236896422, 0.0905748, -0.1508014, 0.6074566741435283, 0.1]
    #x0 = [1.736451985137371, 15.106423236896422, 0.0905748, -0.108014, 0.6074566741435283, 0.1] # ploted to slides
    # >> x0corr = [20.71975544711280, 15.106423236896422, 0.0905748, -0.108014, 0.6074566741435283, 0.1] # ploted to slides
    
    mask = [0, 0, 1, 1, 0, 0] # determining, which parametr should be fitted
    
    x0M = zeros(0)
    lowM = zeros(0)
    uppM = zeros(0)
    for i in collect(1 : 1 : length(mask))
        if convert(Bool,mask[i])
            append!(x0M, x0[i])
            append!(lowM, lower_bounds[i])
            append!(uppM, upper_bounds[i])
        end
    end
    
    PyPlot.close()
    to_optimize(x0M)
    #println(optimize(to_optimize, x0M, lower=lowM, upper=uppM, Δ=1000, f_tol=1.0e-14, g_tol=1.0e-14, LevenbergMarquardt()))
    #optimize(to_optimize, x0M, Dogleg())
    return
end


end

