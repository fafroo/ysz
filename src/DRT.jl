using Printf
using PyPlot
#using Plots
using DataFrames
using LeastSquaresOptim

using NNLS
using LinearAlgebra


#include("../src/simulations/EIS_simulation.jl")

using Random

const DRT_standard_figure = 33

mutable struct DRT_control_struct
  lambda::Float32
  tau_min_fac::Float32
  tau_max_fac::Float32
  tau_range_fac::Float32
  peak_merge_tol::Float32
end

mutable struct DRT_struct
  EIS_df::DataFrame
  tau_range::Array{Float64}
  h::Array{Float64}
  R_ohm::Float64
  L::Float64
  peaks_df::DataFrame
  control::DRT_control_struct
end




function get_R_C_from_DRT(tau_range=Nothing, h_tau=Nothing)
  mcounter = 1
  
  function R_peak(tau, tau_c, h, sigma2, alpha)
    return h*exp.(
      -(
        ((log10.(tau) .- tau_c).*(1 .+ alpha*sign.(log10.(tau) .- tau_c))).^2.0
      )/
      (2*exp(sigma2))
    )
  end
  
  function inte(x)
    R_peak(x, -5, 10, 2, 0.4)
  end
  

  
  if h_tau == Nothing
    tau_range = collect(-10 : 0.001 : 10)
    to_fit = inte(tau_range)
    #plot(tau_range, to_fit)
  else
    to_fit = h_tau
  end

  function discrete_integrate(tau_range, h_tau)
    sum = 0
    for i in 1:length(tau_range)
      sum += h_tau[i]
    end
    return sum
  end
  
  function to_optimize(x)
    res = sqrt(norm(to_fit - R_peak(tau_range, x...)))
    #plot(tau_range, R_peak(tau_range, x...))
    #@show res, x
    res
  end
  
  x_0 = [-3, 1., log(1.0), 0.]
  #to_optimize(x_0)
  
  a = optimize(to_optimize, x_0, BFGS(), Optim.Options(iterations=30))
  #a = optimize(to_optimize, x_0, ConjugateGradient(), Optim.Options(iterations=100))
  #a = optimize(to_optimize, x_0, LevenbergMarquardt(), Optim.Options(iterations=100))
  
  
  #plot(log10.(tau_range), R_peak(tau_range, a.minimizer...))
  
  R_dis = discrete_integrate(tau_range, R_peak(tau_range, a.minimizer...))
  R_orig = discrete_integrate(tau_range, h_tau)
  
  #(R, e) = quadgk(x -> R_peak(x, a.minimizer...), 10^(-4), 10^(-1))
  C_dis = 10^(a.minimizer[1])/R_dis
  C_orig = 10^(a.minimizer[1])/R_orig
  #println(" >>>>>>> R_dis = $(R_dis)  C_dis = $(C_dis)  R_orig = $(R_orig)   C_orig = $(C_orig)  ")
  return R_orig, C_orig
end


function evaluate_RC_peaks_from_DRT(DRT::DRT_struct)
  
  h_max = maximum(DRT.h)
  threshold = h_max/50.0
  
  function peak_assesement(temp_peak_df)
    R = 0.0
    tau_c = 0.0
    for i in 1:size(temp_peak_df, 1)
      R += temp_peak_df.h[i]
      # only an intermediate step ... weighted average of peak points
      tau_c += temp_peak_df.h[i] * temp_peak_df.tau[i]
    end
    tau_c = tau_c/R
    C = tau_c/R
    return tau_c, R, C
  end
  
  #peak finding
  peaks = DataFrame(tau_c = [], R = [], C = [])
  
  temp_peak_df = DataFrame(tau = [], h = [])
  in_peak_bool = false
  for (i, tau) in enumerate(DRT.tau_range)
    if DRT.h[i] > threshold
      push!(temp_peak_df, (tau, DRT.h[i]))
    else
      if size(temp_peak_df, 1) > 0
        push!(peaks, peak_assesement(temp_peak_df))
        
        temp_peak_df = DataFrame(tau = [], h = [])
      end
    end
  end
  
  DRT.peaks_df = peaks
  return
end








# function plot_DRT(DRT::DRT_struct, to_standard_figure=true)
#   if to_standard_figure
#     figure(DRT_standard_figure)
#   end
#   title("DRT")
#   plot(log10.(DRT.tau_range), DRT.h, "-x", label="l=$(DRT.lambda)")
#   xlabel("log10(\$\\tau\$ [s])")
#   ylabel("\$h(\\tau)\$ [Ohm]")
#   if to_standard_figure
#     legend()
#   end
#   peaks_df = get_RC_df_from_DRT(DRT.tau_range, DRT.h)
#   println("DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)")
#   for i in 1:size(peaks_df,1)
#     println(">> R$i = $(peaks_df.R[i])   C$i = $(peaks_df.C[i])    ... (tau_c$i = $(peaks_df.tau_c[i]))")
#   end
# end



function plot_DRT_h(DRT::DRT_struct, to_standard_figure=true, print_bool=false, plot_lambda=true)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  
  title("DRT")
  plot(log10.(DRT.tau_range), DRT.h, "-x", label="\$\\lambda\$=$(DRT.control.lambda)")
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("\$h(\\tau)\$ [Ohm]")
  if (to_standard_figure || plot_lambda)
    legend()
  end
  if print_bool
    println("DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)")  
  end
end
  
function plot_DRT_RC(DRT::DRT_struct, to_standard_figure=true, print_bool=true)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
#   
  title("RC diagram")
  xlabel("C [F]")
  ylabel("R [Ohm]")
  plot(peaks_df.C, peaks_df.R, "x")
  grid(true)
  
  if print_bool
    println("DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)") 
    for i in 1:size(peaks_df,1)
      println(">> f$i = $(1/(2*pi*peaks_df.C[i]*peaks_df.R[i]))    R$i = $(peaks_df.R[i])   C$i = $(peaks_df.C[i])   ... (tau_c$i = $(peaks_df.tau_c[i]))")
    end
  end
end

function plot_DRT_Rtau(DRT::DRT_struct, to_standard_figure=true, print_bool=false)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
  
  title("Rtau diagram")
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("R [Ohm]")
  xlim(log10.([DRT.tau_range[1], DRT.tau_range[end]])...)
  plot(log10.(peaks_df.C.*peaks_df.R), peaks_df.R, "o")
  grid(true)

end

function plot_DRT_Rf(DRT::DRT_struct, to_standard_figure=true, print_bool=false)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
  
  title("Rf diagram")
  xlabel("log10(\$f\$ [Hz])")
  ylabel("R [Ohm]")
  xlim(log10.(1 ./(2*pi*[DRT.tau_range[end], DRT.tau_range[1]]))...)
  plot(log10.(1 ./(2*pi*peaks_df.C.*peaks_df.R)), peaks_df.R, "o")
  grid(true)
  
  if print_bool
    for i in 1:size(peaks_df,1)
      println(">> f$i = $(1/(2*pi*peaks_df.C[i]*peaks_df.R[i]))    R$i = $(peaks_df.R[i])    ... (C$i = $(peaks_df.C[i]))")
    end
  end
end


function get_expspace(A, B, n)
  res = []
  for i in 1:n
    append!(res, A*((Float64(B)/A)^((i - 1)/(n - 1))))
  end
  res
end


function get_DRT(EIS_df::DataFrame, control::DRT_control_struct, debug_mode=false)
  #println(control.lambda)
  #@show tau_max_fac
  tau_min = 1.0/(2*pi*EIS_df.f[end]) / control.tau_min_fac
  tau_max = 1.0/(2*pi*EIS_df.f[1]) * control.tau_max_fac
  tau_range = get_expspace(tau_min, tau_max, control.tau_range_fac*size(EIS_df.f, 1)-2)
  #tau_range = [0.001]
  
  N_f = size(EIS_df.f, 1)
  N_tau = size(tau_range, 1)
  
  # h, R_ohm, L
  n_cols = N_tau + 2  
  # real(Z), imag(Z), regularization
  if control.lambda == 0.0
    n_rows = 2*N_f
  else
    n_rows = 2*N_f + n_cols
  end
#   @show tau_range
#   @show EIS_df.f
  
  A = Matrix{Float64}(undef, n_rows, n_cols)
  b = Vector{Float64}(undef, n_rows)
  
  # assemble "RC" part of A and b
  for (i, f) in enumerate(EIS_df.f)
    #RC
    for (j, tau) in enumerate(tau_range)
      A[i, j]       = real(1/(1 + im*(2*pi*f)*tau))
      A[N_f + i, j] = imag(1/(1 + im*(2*pi*f)*tau))
    end
    # R_ohm
    A[i, N_tau + 1] = 1
    A[N_f + i, N_tau + 1] = 0
    # L
    A[i, N_tau + 2] = 0
    A[N_f + i, N_tau + 2] = 2*pi*f
    
    # b
    b[i] = real(EIS_df.Z[i])
    b[N_f + i] = imag(EIS_df.Z[i])
  end
  
  
#   @show A[N_f,:]
#   @show b[N_f]
#   @show 2222222
#   @show A[2*N_f,:]
#   @show b[2*N_f]
  
  if control.lambda != 0.0
    # assemble "regularization" part of A and b
    A[2*N_f + 1 : end, :] .= Diagonal([control.lambda for i in 1:n_cols])
    b[2*N_f + 1 : end] .= 0
  end
  
  work_l = NNLSWorkspace(A, b);
  solution = solve!(work_l)  
  
  if solution[1] > 0.000001
    println(" !!!!!!!!!!!!!!!!!!!!!!!  Low frequency adjustment ---------------------- ")
    R_ohm_estimate = solution[N_tau + 1] + solution[1]
  
    next_row = zeros(N_tau + 2)
    next_row[N_tau + 1] = 1
    A = [A; transpose(next_row)]
    b = [b; R_ohm_estimate]
    
    work_l = NNLSWorkspace(A, b);
    solution = solve!(work_l)     
  end
  
  if size(A,1) == size(A,2)
#     @show det(A)
#     @show A
#     exact_solution = A \ b
#     @show exact_solution
#     @show solution
#     @show exact_solution - solution
#     println()
#     @show A*exact_solution - b
#     @show A*solution - b
  end
  
#   @show solution
#   
#   @show A \ b
#   @show A
#   @show norm(A*solution - b)
  
  # reconstructin EIS from DRT
  # getting rid of the corner effects
  if control.lambda == 0.0 && false
    solution[end-7 : end-2] .= 0
    solution[1 : 3] .=0
  end
  
  EIS_post = A*solution
  
  EIS_new = deepcopy(EIS_df)
  for i in 1:N_f
    EIS_new.Z[i] = EIS_post[i] + im*EIS_post[N_f + i]
  end
  
  DRT_out = DRT_struct(EIS_new, tau_range, solution[1:end - 2], solution[end-1], solution[end], DataFrame(), control)
  evaluate_RC_peaks_from_DRT(DRT_out)
  
  # the following do NOT change DRT function, only changes peaks_df
  if control.peak_merge_tol > 0.0
    work_to_be_done = true
    new_loop = false
    while work_to_be_done || new_loop
      previous_peak_R = 0
      previous_peak_tau = 10.0^(-20)
      new_loop = false
      for (i, tau) in enumerate(DRT_out.peaks_df.tau_c)
        if abs(log10(previous_peak_tau) - log10(tau)) < control.peak_merge_tol
          
          #
          DRT_out.peaks_df.tau_c[i-1] = 10^((log10(DRT_out.peaks_df.tau_c[i-1] * DRT_out.peaks_df.tau_c[i]))/2)
          DRT_out.peaks_df.R[i-1] += DRT_out.peaks_df.R[i]
          DRT_out.peaks_df.C[i-1] = DRT_out.peaks_df.tau_c[i-1]/DRT_out.peaks_df.R[i-1]
          deleterows!(DRT_out.peaks_df, i)
          #
          
          new_loop = true
        else
          previous_peak_tau = tau
          previous_peak_R = DRT_out.peaks_df.R[i]
        end
      end
      work_to_be_done = false
    end
  end
  
  return DRT_out
end


# function get_A_b()
#   m = rand(20:100)
#   n = rand(20:100)
#   A = randn(m, n)
#   b = randn(m)
#   (A, b)
# end
# 
# function get_sup(A, b, lambda)
#   n = size(b, 1)
#   nuly = [1 for i in 1:n]
#   lambdaI = diagm(n, n, [lambda for i in 1:n])
#   A_sup = vcat(A, lambdaI)
#   b_sup = vcat(b, nuly)
#   (A_sup, b_sup)
# end
# 
# function test_lambda()
#   for i in collect(-6 : 1 : 5) 
#     work_l = NNLSWorkspace(get_sup(A, b, 10.0^i)...); plot(solve!(work_l), label="$i")
#     legend()
#   end
# end


# function assemble_A(f_list, tau_range)
#   N_f = size(f_list, 1)
#   N_tau = size(tau_range, 1)
#   
#   # taus, R_ohm, L
#   n_cols = N_tau + 2
#   
#   # real(Z), imag(Z), regularization
#   n_rows = 2*N_f + n_cols
#   
#   A = Matrix{Float64}(undef, n_rows, n_cols)
# 
#   for (i, f) in enumerate(f_list)
#     #RC
#     for (j, tau) in enumerate(tau_range)
#       A[i, j]       = real(1/(1 + im*(2*pi*f)*tau))
#       A[N_f + i, j] = imag(1/(1 + im*(2*pi*f)*tau))
#     end
#     # R_ohm
#     A[i, N_tau + 1] = 1
#     A[N_f + i, N_tau + 1] = 0
#     # L
#     A[i, N_tau + 2] = 0
#     A[N_f + i, N_tau + 2] = 2*pi*f
#   end
# 
#   # assemble "regularization" part of A
#   A[2*N_f + 1 : end, :] .= diagm(n_cols, n_cols, [lambda for i in 1:n_cols])
#   
#   return A
# end
# 
# function assemble_b(f_list, Z_list)
#   N_f = size(f_list, 1)
#   N_tau = size(tau_range, 1)
#   
#   # taus, R_ohm, L
#   n_cols = N_tau + 2
#   
#   # real(Z), imag(Z), regularization
#   n_rows = 2*N_f + n_cols
#   
#   b = Vector{Float64}(undef, n_rows)
#   
#   for (i, f) in enumerate(EIS_df.f)
#     #RC
#     for (j, tau) in enumerate(tau_range)
#       A[i, j]       = real(1/(1 + im*(2*pi*f)*tau))
#       A[N_f + i, j] = imag(1/(1 + im*(2*pi*f)*tau))
#     end
#     # R_ohm
#     A[i, N_tau + 1] = 1
#     A[N_f + i, N_tau + 1] = 0
#     # L
#     A[i, N_tau + 2] = 0
#     A[N_f + i, N_tau + 2] = 2*pi*f
#     
#     # b
#     b[i] = real(EIS_df.Z[i])
#     b[N_f + i] = real(EIS_df.Z[i])
#   end
#   
#   # assemble "regularization" part of b
#   b[2*N_f + 1 : end] .= 0
#   
#   return b
# end

function EIS_get_RC_CPE_elements(R1, C1, R2, C2, alpha, Rohm=0; f_range=Nothing)
  EIS_RC = DataFrame( f = [], Z = [])
  if f_range==Nothing
    checknodes = get_shared_checknodes(EIS_simulation(800,100,0.0)...)
  else
    checknodes = EIS_get_checknodes_geometrical(f_range...)
  end
  for f in checknodes
    push!(
      EIS_RC, 
      (f, Rohm + R1/(1 + im*2*pi*f*R1*C1) + R2/(1 + ((im*2*pi*f*R2*C2)^alpha))  )
    ) 
  end
  return EIS_RC  
end


