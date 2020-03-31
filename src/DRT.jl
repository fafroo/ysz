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
  tau_list::Array{Float64}
  h::Array{Float64}
  R_ohm::Float64
  L::Float64
  lambda::Float64
end


function plot_DRT(DRT::DRT_struct, to_standard_figure=true)
  if to_standard_figure
    figure(33)
  end
  title("DRT")
  plot(log10.(DRT.tau_list), DRT.h, "-x", label="l=$(DRT.lambda)")
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("\$h(\\tau)\$ [Ohm]")
  if to_standard_figure
    legend()
  end
  println("R_ohm = $(DRT.R_ohm)    L = $(DRT.L)")
end


function get_expspace(A, B, n)
  res = []
  for i in 1:n
    append!(res, A*((Float64(B)/A)^((i - 1)/(n - 1))))
  end
  res
end


function get_DRT(EIS_df::DataFrame, lambda=0.0)
  println(lambda)
  tau_min = 1.0/(2*pi*EIS_df.f[end]) / 0.1
  tau_max = 1.0/(2*pi*EIS_df.f[1]) * 1
  tau_list = get_expspace(tau_min, tau_max, 5*size(EIS_df.f, 1))
  #tau_list = [0.001]
  
  N_f = size(EIS_df.f, 1)
  N_tau = size(tau_list, 1)
  
  # h, R_ohm, L
  n_cols = N_tau + 2  
  # real(Z), imag(Z), regularization
  if lambda == 0.0
    n_rows = 2*N_f
  else
    n_rows = 2*N_f + n_cols
  end
#   @show tau_list
#   @show EIS_df.f
  
  A = Matrix{Float64}(undef, n_rows, n_cols)
  b = Vector{Float64}(undef, n_rows)
  
  # assemble "RC" part of A and b
  for (i, f) in enumerate(EIS_df.f)
    #RC
    for (j, tau) in enumerate(tau_list)
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
    A[2*N_f + 1 : end, :] .= diagm(n_cols, n_cols, [lambda for i in 1:n_cols])
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
  
  DRT_out = DRT_struct(EIS_new, tau_list, solution[1:end - 2], solution[end-1], solution[end], lambda)
  
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


# function assemble_A(f_list, tau_list)
#   N_f = size(f_list, 1)
#   N_tau = size(tau_list, 1)
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
#     for (j, tau) in enumerate(tau_list)
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
#   N_tau = size(tau_list, 1)
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
#     for (j, tau) in enumerate(tau_list)
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



