# module mesonbse
using ProgressMeter
# using ProgressBars
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


# Read quark system

now_time = Dates.format(Dates.now(), "YYYY-mm-dd-HH:MM")
################################
##这里使用程慕阳的程序,请更改原文件##
################################

# quarkrepoint = dataset["quarkDSE"]["repoint"]
# quarkintstep = dataset["quarkDSE"]["quarkintstep"]
# quarkm = dataset["quarkDSE"]["quarkmass"]

# data
logofcutoff = dataset["mesonBSE"]["logofcutoff"]
mt = dataset["data"]["mt"]
τ = dataset["data"]["tau"]
Λ = dataset["data"]["lambda"]
ω = dataset["data"]["omega"]
dd = dataset["data"]["dd"]
Nf = dataset["data"]["Nf"]
rm = dataset["data"]["rm"]
Parameters_Xi = dataset["data"]["Parameters_Xi"]
z2 = dataset["data"]["z2"]
z4 = dataset["data"]["z4"]
cutup = 10. ^ logofcutoff
cutdown = 10. ^ (-4.)

kstep = dataset["movingframeBSE"]["kstep"]
z1step = dataset["movingframeBSE"]["z1step"]
z2step = dataset["movingframeBSE"]["z2step"]
Pstep = dataset["movingframeBSE"]["Pstep"]
zintstep = dataset["movingframeBSE"]["zintstep"]
massRange = dataset["movingframeBSE"]["massRange"]
externalP2 = dataset["movingframeBSE"]["externalP2"]
dim = kstep * z1step * z2step
# 取点方式
if Pstep == 1
    plist = [massRange[1]]
else
    plist = [massRange[1] + (massRange[2]-massRange[1])/(Pstep-1) * t for t in 0:(Pstep-1)]
end

# meshk,weightk = complex.(gausslegendremesh(cutdown,cutup,kstep,2));
# meshz,weightz = complex.(gausschebyshev(zstep,2));
meshk,weightk = gausslegendremesh(cutdown,cutup,kstep,2)
meshz1,weightz1 = gausschebyshev(z1step,2)
meshz2,weightz2 = gausslegendre(z2step)
# meshz2,weightz2 = gausschebyshev(z2step,2)
D(t) = 8 * pi^2 * (dd * exp(-t/(ω^2)) / ω^4 + rm * ( (-expm1(-t/(2 * mt)^2))/t )/ log(τ + (1+t/Λ^2)^2))
D_infrared = D
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
# branchfunction(x)=(x*AA(x)^2+BB(x)^2)

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

computeIndex(a::Int, b::Int) = (a - 1) * zstep + b  
computeIndex(a::Int, b::Int, c::Int) = (a - 1) * z1step * z2step + (b - 1) * z2step + c
# computeIndex(i,s)

print("参数导入完毕,开始计算mesonBSA\n")
if dataset["calcType"]["mesonmode"] == 1
    println("将计算psmeson")
    include(joinpath(workdir,"src/equations/psmeson_mf.jl"))
elseif dataset["calcType"]["mesonmode"] == 2
    println("将计算scalarmeson")
    include(joinpath(workdir,"src/equations/scalarmeson_mf.jl"))
elseif dataset["calcType"]["mesonmode"] == 3
    println("将计算vectormeson")
    include(joinpath(workdir,"src/equations/vectormeson_mf.jl"))
elseif dataset["calcType"]["mesonmode"] == 4
    println("将计算avmeson")
    include(joinpath(workdir,"src/equations/avmeson_mf.jl"))
end # if for mode

# end # module
