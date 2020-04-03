using Printf
using PyPlot
#using Plots
using DataFrames
using LeastSquaresOptim

using NNLS
using LinearAlgebra


#include("../src/simulations/EIS_simulation.jl")

using Random

mutable struct DRT_struct
  EIS_df::DataFrame
  tau_range::Array{Float64}
  h::Array{Float64}
  R_ohm::Float64
  L::Float64
  lambda::Float64
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


function get_RC_df_from_DRT(tau_range=Nothing, h_tau=Nothing)
  h_max = maximum(h_tau)
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
  for (i, tau) in enumerate(tau_range)
    if h_tau[i] > threshold
      push!(temp_peak_df, (tau, h_tau[i]))
    else
      if size(temp_peak_df, 1) > 0
        push!(peaks, peak_assesement(temp_peak_df))
        
        temp_peak_df = DataFrame(tau = [], h = [])
      end
    end
  end

  return peaks
end








# function plot_DRT(DRT::DRT_struct, to_standard_figure=true)
#   if to_standard_figure
#     figure(33)
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



function plot_DRT_h(DRT::DRT_struct, to_standard_figure=true)
  if to_standard_figure
    figure(33)
  end
  
  suptitle("DRT plot")
  plot(log10.(DRT.tau_range), DRT.h, "-x", label="l=$(DRT.lambda)")
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("\$h(\\tau)\$ [Ohm]")
  if to_standard_figure
    legend()
  end
end
  
function plot_DRT_RC(DRT::DRT_struct, to_standard_figure=true, print_bool=true)
  if to_standard_figure
    figure(33)
  end
  peaks_df = get_RC_df_from_DRT(DRT.tau_range, DRT.h)
  
  #title("RC characteristics")
  xlabel("R [Ohm]")
  ylabel("C [F]")
  plot(peaks_df.R, peaks_df.C, "x")
  grid(true)
  
  if print_bool
    println("DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)")
    for i in 1:size(peaks_df,1)
      println(">> R$i = $(peaks_df.R[i])   C$i = $(peaks_df.C[i])    ... (tau_c$i = $(peaks_df.tau_c[i]))")
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


function get_DRT(EIS_df::DataFrame, lambda=0.0)
  #println(lambda)
  tau_min = 1.0/(2*pi*EIS_df.f[end]) / 0.1
  tau_max = 1.0/(2*pi*EIS_df.f[1]) * 1
  tau_range = get_expspace(tau_min, tau_max, 5*size(EIS_df.f, 1))
  #tau_range = [0.001]
  
  N_f = size(EIS_df.f, 1)
  N_tau = size(tau_range, 1)
  
  # h, R_ohm, L
  n_cols = N_tau + 2  
  # real(Z), imag(Z), regularization
  if lambda == 0.0
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
  
  if lambda != 0.0
    # assemble "regularization" part of A and b
    A[2*N_f + 1 : end, :] .= Diagonal([lambda for i in 1:n_cols])
    b[2*N_f + 1 : end] .= 0
  end
  
  work_l = NNLSWorkspace(A, b);
  solution = solve!(work_l)  
  
#   @show x2
#   @show solution
#   
#   @show A \ b
#   @show A
#   @show norm(A*solution - b)
  
  # reconstructin EIS from DRT
  EIS_post = A*solution
  
  EIS_new = deepcopy(EIS_df)
  for i in 1:N_f
    EIS_new.Z[i] = EIS_post[i] + im*EIS_post[N_f + i]
  end
  
  DRT_out = DRT_struct(EIS_new, tau_range, solution[1:end - 2], solution[end-1], solution[end], lambda)
  
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



