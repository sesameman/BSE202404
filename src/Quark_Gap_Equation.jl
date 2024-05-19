module Quark_Gap_Equation
export z2

# module ComplexDSE

# greet() = print("Hello World!")

# end # module complexDSE
# ---------------------------------- Module ---------------------------------- #
# using BenchmarkTools
using Printf
using LinearAlgebra
using Gaussquad
using FastGaussQuadrature
import Base.:+
# ---------------------------------------------------------------------------- #



elapsed_time = @elapsed begin
# --------------------------------- Constant --------------------------------- #
const mt = 0.5;
const τ = ℯ^2-1;
const Λ = 0.234;
const ω = 0.5;
const dd = (0.87)^3/ω;
const Nf = 4;
const rm = 12/(33 - 2*Nf);
renorm = 19.
# ---------------------------------------------------------------------------- #

# -------------------------------- Parameters -------------------------------- #
z2 = 0.991882
z4 = 0.8
currentQuarkMass = 0.004;
# quadrature step
parabolastep = 100 # 抛物线下半部分取点
verticlestep = 30  # uv cut upper half 
cauchystep = parabolastep + verticlestep
intQstep = 32
intZstep = 12
uv = 10000.
ir = 0.0001
x0_complex = 0.4
bigM = x0_complex ^ 2
# ? 这里用列表可以整个插入complex函数里
iniAnum = [1.6]
iniBnum = [0.8]
# ? 迭代参数
iteration = 0
max_iterations = 100
tolerance_ab = 10^(-5)
tolerance_z2 = 10^(-2)
tolerance_fab = 0.0001

# ---------------------------------------------------------------------------- #

# ----------------------------------- Gluon ---------------------------------- #
F(x) = (1 - exp(-x / (2*mt)^2)) / x;
Deff(t) = 8 * pi^2 * (dd * exp(-t/ω^2) / ω^4 + rm * F(t) / log(τ + (1 + t/Λ^2)^2));
# ---------------------------------------------------------------------------- #

# ------------------------- Struct & Operation rules ------------------------- #
struct Parameters
    x::Vector{Float64}
    w::Vector
    points::Vector{ComplexF64}
end
function +(a::Parameters, b::Parameters)::Parameters
    Parameters([a.x ; b.x], [a.w ; b.w], [a.points ; b.points])
end
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
#                                 Curve set up                                 #
# ---------------------------------------------------------------------------- #
global curvePdown
let (x_local,w_local) = gausslegendremesh(ir, uv + bigM, parabolastep, 2) 
    w_out = Vector{ComplexF64}(undef, length(w_local))
    for i in eachindex(x_local)
        w_out[i] = w_local[i] * ( 1 - im * sqrt(bigM/x_local[i]) )
    end
    global curvePdown = Parameters(x_local, w_out, [- bigM + t - 2 * im * sqrt(bigM * t) for t in x_local])
end
global curveVup
let (x_local,w_local) = gausslegendremesh(0.0, 1.0, verticlestep, 1) 
    w_out = Vector{ComplexF64}(undef, length(w_local))
    for i in eachindex(x_local)
        w_out[i] = w_local[i] * (2 * im  * sqrt((bigM + uv) * bigM))
    end
    global curveVup = Parameters(x_local, w_out, [uv + 2 * im * sqrt(bigM * (bigM + uv)) * t  for t in x_local])
end
curveP2 = curvePdown + curveVup
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #




# ------------------------------ Initialization ------------------------------ #
intQ, WintQ = gausslegendremesh(ir, uv, intQstep, 2) 
intZ ,WintZ = gausschebyshev(intZstep,2)
maxofintQ = maximum(intQ)
maxofcurveP = maximum(curveP2.x) - bigM
kernelAlist_points = Array{ComplexF64}(undef, cauchystep, intQstep, intZstep)
kernelAlist_renorm = Array{ComplexF64}(undef, intQstep, intZstep)
Gluonlist = Array{Float64}(undef, intQstep)
# kernelBlist = Array{ComplexF64}(undef, cauchystep, intQstep, intZstep)
AB_ini = [fill(complex(iniAnum...), cauchystep) ;  fill(complex(iniBnum...), cauchystep) ]
J_ini = Array{ComplexF64}(I, 2*cauchystep, 2*cauchystep)
# ---------------------------------------------------------------------------- #

# ------------------------------ Precomputation ------------------------------ #
for i in eachindex(curveP2.points)
    for j in eachindex(intQ)
        for k in eachindex(intZ)
            local k2
            p2 = curveP2.points[i]
            q2 = intQ[j]
            z = intZ[k]
            kernelAlist_points[i, j, k] = 1 + 2 * z^2 - 3 * sqrt(q2)/sqrt(p2) * z
            if i == 1
                p2 = renorm^2
                kernelAlist_renorm[j, k] = 1 + 2 * z^2 - 3 * sqrt(q2)/sqrt(p2) * z
            end
        end
    end
end

for i in eachindex(intQ)
    Gluonlist[i] = Deff(intQ[i]) 
end
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
#                                 Calc Function                                #
# ---------------------------------------------------------------------------- #
function Fab(AandB; 
    curve = curveP2, cauchystep = cauchystep,  
    Nq = intQstep, q2 = intQ, wq = WintQ, 
    Nz = intZstep, z = intZ, wz = WintZ, 
    Gluon = Gluonlist, currentQuarkMass = currentQuarkMass, 
    kernelA = kernelAlist_points, kernelA_re = kernelAlist_renorm,
    z2 = z2, z4 = z4, renorm = renorm, maxCauchyQ2 = (maxofintQ)
    )
    # * 每一次计算都要重新计算A、B在不同k²的值
    # * 得到A、B的list就可以计算Σab, 进而计算Fab
    # complexF = zeros(length(x))
    A_points = AandB[1:cauchystep]
    B_points = AandB[(1 + cauchystep):(2*cauchystep)]
    complexF = zeros(ComplexF64, 2 * cauchystep)
    Σa_renorm = complex(0.)
    Σb_renorm = complex(0.)
    a = 0
    b = 0
    # ? Calculate Σa_renorm and Σb_renorm
    for i in 1:Nq
        for j in 1:Nz
            # TODO: 这里计算重整化点的值, 将p²改为renorm²
            p2 = renorm^2
            pdotq = sqrt(p2) * sqrt(q2[i]) * z[j]
            k2 = p2 + q2[i] - 2 * pdotq
            # A = cauchyTheoremConj_old(A_points, curve.points, curve.w, k2)
            # B = cauchyTheoremConj_old(B_points, curve.points, curve.w, k2)
            A = cauchyTheoremConj_old(A_points, curve.points, curve.w, k2, maxCauchyQ2)
            B = cauchyTheoremConj_old(B_points, curve.points, curve.w, k2, maxCauchyQ2)
            # A = cauchyTheoremConj(A_points, k2)
            # B = cauchyTheoremConj(B_points, k2)
            # A = cauchyTheoremConj(A_points, k2, z = curve.points, w = curve.w, maxofintQ = maxofintQ)
            # B = cauchyTheoremConj(B_points, k2, z = curve.points, w = curve.w, maxofintQ = maxofintQ)
            samepart = q2[i] * wq[i] * wz[j] * Gluon[i] / (k2 * A^2 + B^2) 
            Σa_renorm += samepart * A * kernelA_re[i,j]
            Σb_renorm += samepart * B * 3
        end
    end
    commonpart = (2*pi)^(-3)
    Σa_renorm *= commonpart
    Σb_renorm *= commonpart
    for pstep = 1:cauchystep
        Σa = complex(0.)
        Σb = complex(0.)
        for i in 1:Nq
            for j in 1:Nz
                p2 = curve.points[pstep]
                pdotq = sqrt(p2) * sqrt(q2[i]) * z[j]
                k2 = p2 + q2[i] - 2 * pdotq
                # A = cauchyTheoremConj_old(A_points, curve.points, curve.w, k2)
                # B = cauchyTheoremConj_old(B_points, curve.points, curve.w, k2)
                A = cauchyTheoremConj_old(A_points, curve.points, curve.w, k2, maxCauchyQ2)
                B = cauchyTheoremConj_old(B_points, curve.points, curve.w, k2, maxCauchyQ2)
                # A = cauchyTheoremConj(A_points, k2)
                # B = cauchyTheoremConj(B_points, k2)
                # A = cauchyTheoremConj(A_points, k2, z = curve.points, w = curve.w, maxofintQ = maxofintQ)
                # B = cauchyTheoremConj(B_points, k2, z = curve.points, w = curve.w, maxofintQ = maxofintQ)
                samepart = q2[i] * wq[i] * wz[j] * Gluon[i] / (k2 * A^2 + B^2) 
                Σa += samepart * A * kernelA[pstep,i,j]
                Σb += samepart * B * 3
            end
        end
        Σa *= commonpart
        Σb *= commonpart
        # if pstep == 1
        #     println(Σa)
        #     println(Σb)
        # end   
        complexF[pstep] = AandB[pstep]- z2 - 4/3 * z2^2 * Σa # 这里记得修改成curve的形式
        complexF[pstep + cauchystep] = AandB[pstep + cauchystep] - currentQuarkMass - 4/3 * z2^2 * (Σb - Σb_renorm)
    end
    return complexF, Σa_renorm, Σb_renorm
end
# ---------------------------------------------------------------------------- #


# ------------------------------- Math Function ------------------------------ #
function cauchyTheoremConj_old(F::Array{Complex{Float64}}, z::Array{Complex{Float64}}, w::Array{Complex{Float64}}, x₀::Number)
    result1 = complex(0.)
    result2 = complex(0.)
    x₀ = Complex{Float64}(x₀)
    for i in eachindex(F)
        result1 += F[i] * w[i] / (z[i]-x₀)
        result1 -= conj(F[i]) * conj(w[i]) / (conj(z[i])-x₀)
        result2 += w[i] / (z[i]-x₀)
        result2 -= conj(w[i]) / (conj(z[i])-x₀)
    end
    # return result / (2*pi*im)
    return result1/result2
end

# ! 派发一个带判断的柯西积分
function cauchyTheoremConj_old(F::Array{Complex{Float64}}, z::Array{Complex{Float64}}, w::Array{Complex{Float64}}, x₀::Number, maxofintQ::Float64)
    result1 = 0. + im * 0. 
    result2 = 0. + im * 0.
    if real(x₀) > maxofintQ
        x₀ = maxofintQ + imag(x₀)
    end
    for i in eachindex(F)
        result1 += F[i] * w[i] / (z[i]-x₀)
        result1 -= conj(F[i]) * conj(w[i]) / (conj(z[i])-x₀)
        result2 += w[i] / (z[i]-x₀)
        result2 -= conj(w[i]) / (conj(z[i])-x₀)
    end
    # return result / (2*pi*im)
    return result1/result2
end

# ! 这个柯西积分优雅一点, 但不知道为啥变慢了
# function cauchyTheoremConj(F, x₀; z = curveP2.points, w  = curveP2.w, maxq2 = maxofintQ)
#     result1 = 0. + im * 0. 
#     result2 = 0. + im * 0.
#     if real(x₀) > maxq2
#         x₀ = maxq2 + imag(x₀)
#     end
#     # if real(x₀) > uv
#     #     x₀ = uv + imag(x₀)
#     # end
#     # z = curve.points
#     # w = curve.w
#     for i in eachindex(F)
#         result1 += F[i] * w[i] / (z[i]-x₀)
#         result1 -= conj(F[i]) * conj(w[i]) / (conj(z[i])-x₀)
#         result2 += w[i] / (z[i]-x₀)
#         result2 -= conj(w[i]) / (conj(z[i])-x₀)
#     end
#     # return result / (2*pi*im)
#     return result1/result2
# end

function Broydenstep(F,F_ini,J_ini,x_ini)
    # * Bad Broyden’s Method
    # 这里F生成的是Vector列向量
    # 输入的x也应该是列向量
    x_new = x_ini - J_ini * F_ini
    delta_x = x_new - x_ini
    F_new, Σa_renorm, Σb_renorm = F(x_new)
    delta_f = F_new - F_ini
    J_new = J_ini + (delta_x - J_ini*delta_f) * (delta_f') / dot(delta_f,delta_f)
    return  F_new, J_new, x_new, Σa_renorm, Σb_renorm
end
# ---------------------------------------------------------------------------- #


# ------------------------------- Old/Test loop ------------------------------ #
# * 这里是旧的循环, 已经遗弃
# for i = 1:max_iterations
#     global z2, z4, Fab_ini, J_ini, AB_ini, k19sumA, k19sumB
#     Fab_ini, J_ini, AB_ini, k19sumA, k19sumB = Broydenstep(Fab,Fab_ini,J_ini,AB_ini)
#     z2 = (-1 + sqrt(1 + 4*k19sumA)) / (2*k19sumA)
#     # z4 = (currentQuarkMass - z2^2 * k19sumB) / currentQuarkMass
#     println("i = $i ---->  ", maximum(real.(Fab_ini).^2),"   vs   ", maximum(abs.(Fab_ini)))
#     if maximum(real.(Fab_ini).^2)<0.000001
#         break
#     end
# end
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
#                                     Main                                     #
# ---------------------------------------------------------------------------- #
Fab_ini = Fab(AB_ini)[1]
while true
    global iteration
    global z2, z4, Fab_ini, J_ini, AB_ini, k19sumA, k19sumB
    iteration += 1
    # ! 这样计算delta_ab_max需要提前算
    delta_ab_max = round(maximum(abs.(J_ini * Fab_ini)), sigdigits=4)
    Fab_ini, J_ini, AB_ini, k19sumA, k19sumB = Broydenstep(Fab,Fab_ini,J_ini,AB_ini)
    z2_new = (-1 + sqrt(1 + 4*k19sumA)) / (2*k19sumA)
    delta_z2_max = round(abs(z2_new - z2), sigdigits=4)
    Fab_max = maximum(abs.(Fab_ini))
    z2 = z2_new
    # 检查收敛条件
    if Fab_max < tolerance_fab
        println("delta_ab_max = $delta_ab_max, delta_z2_max = $delta_z2_max")
        print("\033[31m Fab-Interrupt: \033[0m")
        print("The maximum value of Fab is ")
        @printf("%.6e", Fab_max)
        println(", stop iterating.")
        break
    elseif delta_ab_max < tolerance_ab && delta_z2_max < tolerance_z2
        println("delta_ab_max = $delta_ab_max, delta_z2_max = $delta_z2_max")
        print("The maximum value of Fab is ")
        @printf("%.6e.\n", Fab_max)
        print("\033[31m Delta-Interrupt: \033[0m")
        println("The solution to the equation remains essentially unchanged, stop iterating.")
        break
    elseif iteration >= max_iterations
        println("delta_ab_max = $delta_ab_max, delta_z2_max = $delta_z2_max")
        print("The maximum value of Fab is ")
        @printf("%.6e.\n", Fab_max)
        print("\033[31m End-Interrupt: \033[0m")
        println("Reached the $iteration, loop ends.")
        break
    else
        print("Error")
        @printf("%4i",iteration)
        print("  -----> ")
        @printf("%.6e\n", Fab_max)
    end
end
# ---------------------------------------------------------------------------- #
#                                                                              #
# ---------------------------------------------------------------------------- #



end # time Calculate







print_seconds(elapsed_time)

end

