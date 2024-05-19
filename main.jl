time_taken = @elapsed begin # calculat time

workdir = @__DIR__
push!(LOAD_PATH, joinpath(workdir,"module"))
push!(LOAD_PATH, joinpath(workdir,"src"))
# include(joinpath(@__DIR__,"module/testpkg.jl")) # 检测包完整性
# include(joinpath(@__DIR__,"src/equations/quark_gap.jl"))
# include(joinpath(@__DIR__,"module/inter.jl"))
# using .readdse
using testpkg
using dse_add_on


using TOML
using ProgressMeter
using LinearAlgebra
using JLD2
using FastGaussQuadrature
using Gaussquad # View in github: https://github.com/kangjiayin/Gaussquad.jl
using DataFrames
using CSV
using Arpack # 大型稀疏矩阵 eigen
using Dates


using Quark_Gap_Equation
dataset = TOML.parsefile("config.toml")
calcType = dataset["calcType"]["calcType"]

# if calcType == "DSE"
#     # 执行DSE相关的计算
# elseif calcType == "BSE"
#     include(joinpath(workdir,"src/BSE.jl"))
# elseif calcType == "BSEmovingframe"
#     include(joinpath(workdir,"src/BSEmovingframe.jl"))
# elseif calcType == "DecayConstant"
#     # 执行衰变常数的计算
# else
#     error("Unknown calculation type")
# end    








end # time calculation
println("taken time: $time_taken 秒")


