time_taken = @elapsed begin
    
    
    
using TOML
workdir = @__DIR__
# include(joinpath(@__DIR__,"module/testpkg.jl")) # 检测包完整性
# include(joinpath(@__DIR__,"src/equations/quark_gap.jl"))
# include(joinpath(@__DIR__,"module/inter.jl"))

const IS_MODULE_LOADED = Ref(false)

if !IS_MODULE_LOADED[]
    include(joinpath(workdir,"module/readdse.jl"))
    IS_MODULE_LOADED[] = true
    println("模块已加载。")
else
    println("模块已经在之前的REPL会话中加载过了。")
end

dataset = TOML.parsefile("config.toml")
calcType = dataset["calcType"]["calcType"]

if calcType == "DSE"
    # 执行DSE相关的计算
elseif calcType == "BSE"
    include(joinpath(workdir,"src/BSE.jl"))
elseif calcType == "BSEmovingframe"
    include(joinpath(workdir,"src/BSEmovingframe.jl"))
elseif calcType == "DecayConstant"
    # 执行衰变常数的计算
else
    error("Unknown calculation type")
end    








end
println("taken time: $time_taken 秒")


