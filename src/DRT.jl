using Printf
using PyPlot
#using Plots
using DataFrames
using LeastSquaresOptim

using NNLS

using Random


function get_A_b()
  m = rand(20:100)
  n = rand(20:100)
  A = randn(m, n)
  b = randn(m)
  (A, b)
end

function get_sup(A, b, lambda)
  nuly = [1 for i in 1:n]
  lambdaI = diagm(n, n, [lambda for i in 1:n])
  A_sup = vcat(A, lambdaI)
  b_sup = vcat(b, nuly)
  (A_sup, b_sup)
end

function test_lambda()
  for i in collect(-6 : 1 : 5) 
    work_l = NNLSWorkspace(get_sup(A, b, 10.0^i)...); plot(solve!(work_l), label="$i")
    legend()
  end
end

function get_DRT(EIS_df::DataFrame)
  tau_min = 1.0/(2*pi*EIS_df.f[end]) - 1
  tau_max = 1.0/(2*pi*EIS_df.f[end]) + 1

end
