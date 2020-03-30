#!/usr/local/pkg/Installs/linux-ubuntu18.04-x86_64-gcc7.4.0/bin-julia/1.0.5/5fgzaay5eavbfm4v/bin/julia
#
#SBATCH --job-name=NYQ
#SBATCH -p express3
#SBATCH --time=0-6

include(string(pwd(),"/../src/ysz_fitting.jl"))

@time ysz_fitting.run_par_study_script_wrap(ARGS[1], ARGS[2], ARGS[3], ARGS[4], ARGS[5], ARGS[6], ARGS[7], ARGS[8], ARGS[9], ARGS[10], ARGS[11])
