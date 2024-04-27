# module mesonbse
using ProgressMeter
using ProgressBars
using TOML
using LinearAlgebra
using Dierckx
using JLD2
using Gaussquad # View in github: https://github.com/kangjiayin/Gaussquad.jl
using FastGaussQuadrature
# using ChebyshevFu
using DataFrames
using CSV
using Arpack
using Dates

dataset = TOML.parsefile("config.toml")
# delete!(dataset,"owner")
# delete!(dataset,"title")
delete!(dataset,"readsetting")
# hashname = hash(dataset)

# Read quark system

now_time = Dates.format(Dates.now(), "YYYY-mm-dd-HH:MM")
################################
##这里使用程慕阳的程序,请更改原文件##
################################

quarkrepoint = dataset["quarkDSE"]["repoint"]
quarkintstep = dataset["quarkDSE"]["quarkintstep"]
quarkm = dataset["quarkDSE"]["quarkmass"]

# data
logofcutoff = dataset["mesonBSE"]["logofcutoff"]
mt = dataset["data"]["mt"]
τ = dataset["data"]["tau"]
Λ = dataset["data"]["lambda"]
ω = dataset["data"]["omega"]
dd = dataset["data"]["dd"]
Nf = dataset["data"]["Nf"]
rm = dataset["data"]["rm"]
z2 = dataset["data"]["z2"]
z4 = dataset["data"]["z4"]
cutup = 10. ^logofcutoff
cutdown = 10. ^(-4)

kstep = dataset["mesonBSE"]["kstep"]
zstep = dataset["mesonBSE"]["zstep"]
Pstep = dataset["mesonBSE"]["Pstep"]
massRange = dataset["mesonBSE"]["massRange"]
dim = kstep * zstep
# 取点方式
plist = [massRange[1] + (massRange[2]-massRange[1])/(Pstep-1) * t for t in 0:(Pstep-1)]

meshk,weightk = gausslegendremesh(cutdown,cutup,kstep,2);
meshz,weightz = gausschebyshev(zstep,2);

D(t) = 8 * pi^2 * (dd * exp(-t/(ω^2)) / ω^4 + rm * ( (-expm1(-t/(2 * mt)^2))/t )/ log(τ+(1+t/Λ^2)^2))
# D(t::Float64)=8*pi^2*(rm/log(τ+(1+t/Λ^2)^2))*(1/(t+(ω^4/(t+ω^2))))*(1+dd*ω^2/(t+(ω^4/(t+ω^2)))) # QC-2-04-12


function Inport()
    global z2, z4, AA1, BB1
    local A, B
    AA1 = Main.readdse.realA
    BB1 = Main.readdse.realB
    mk = Main.readdse.uv
    return mk
end

maxofk = Inport()
function AA(x)
    if real(x) > maxofk
        x = complex(maxofk, imag(x))
    end
    return AA1(x)
end
function BB(x)
    if real(x) > maxofk
        x = complex(maxofk, imag(x))
    end
    return BB1(x)
end
# AA(x) = AA1(x)
# BB(x) = BB1(x)
branchfunction(x::Float64)=(x*AA(x)^2+BB(x)^2)

# 处理实数输入
function safe_sqrt(x::Real)
    if x < 0
        return sqrt(Complex(x))
    else
        return sqrt(x)
    end
end

# 处理复数输入
safe_sqrt(x::Complex) = sqrt(x)

print("参数导入完毕,开始计算mesonBSA\n")
if dataset["mesonBSE"]["mesonmode"] == 1
    println("将计算psmeson")
    include(joinpath(pwd(),"src/equations/psmeson.jl"))
elseif dataset["mesonBSE"]["mesonmode"] == 2
    println("将计算scalarmeson")
    include(joinpath(pwd(),"src/equations/scalarmeson.jl"))
elseif dataset["mesonBSE"]["mesonmode"] == 3
    println("将计算vectormeson")
    include(joinpath(pwd(),"src/equations/vectormeson.jl"))
elseif dataset["mesonBSE"]["mesonmode"] == 4
    println("将计算avmeson")
    include(joinpath(pwd(),"src/equations/avmeson.jl"))
end # if for mode

# end # module
