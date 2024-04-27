module readdse
# 使用为Main.readdse.realA(10)
export realA, realB, uv
export ExtraA, ExtraB



using Interpolations
using MathLink
# using CSV
# using DataFrames

# # 读取CSV文件到DataFrame
# filepath = "/home/kjy/dse2/build/Output/DS_RL_0_cmplx.dat"
# df = CSV.read(filepath, DataFrame)
# 导入DelimitedFiles模块
import Base: +
using DelimitedFiles
using Gaussquad
using FastGaussQuadrature
using Dierckx
using Plots
using CSV
using DataFrames
using Interpolations
# 定义文件路径
# print(@__DIR__)
file_path = joinpath(@__DIR__,"./moduleData/DS_RL_0_cmplx.dat") ## complex

uv = 10000.
ir = 0.0001
x0_complex = 0.2
bigM = x0_complex ^ 2 
parabolastep = 100
verticlestep = 30



df = DataFrame(CSV.File(file_path; delim='\t', header=true, types=repeat([Float64], 8)))
# 假设第一部分的数据结束于第72行
df_part1 = df[1:parabolastep, :]
df_part2 = df[(parabolastep+2):end, :]

# 本打算做线性外推一下, 但用不到了
# etp_rAA1 = linear_interpolation(df_part1[!,1].^2, df_part1[!,2],extrapolation_bc=Line())
# etp_iAA1 = linear_interpolation(df_part1[!,1].^2, df_part1[!,3],extrapolation_bc=Line())
# etp_rBB1 = linear_interpolation(df_part1[!,1].^2, df_part1[!,4],extrapolation_bc=Line())
# etp_iBB1 = linear_interpolation(df_part1[!,1].^2, df_part1[!,5],extrapolation_bc=Line())

# etp_rAA2 = linear_interpolation(df_part2[!,1], df_part2[!,2],extrapolation_bc=Line())
# etp_iAA2 = linear_interpolation(df_part2[!,1], df_part2[!,3],extrapolation_bc=Line())
# etp_rBB2 = linear_interpolation(df_part2[!,1], df_part2[!,4],extrapolation_bc=Line())
# etp_iBB2 = linear_interpolation(df_part2[!,1], df_part2[!,5],extrapolation_bc=Line())


struct Parameters
    x::Vector{Float64}
    w::Vector
    points::Vector{ComplexF64}
end

function +(a::Parameters, b::Parameters)::Parameters
    Parameters([a.x ; b.x], [a.w ; b.w], [a.points ; b.points])
end

# # 弃用的完全柯西积分
# function cauchyTheorem(F, z, w, x₀)
#     result1 = 0. + im * 0. 
#     result2 = 0. + im * 0. 
#     for i in eachindex(F)
#         result1 += F[i] * w[i] / (z[i]-x₀)
#         result2 += w[i] / (z[i]-x₀)
#     end
#     # return result / (2*pi*im)
#     return result1/result2
# end


function cauchyTheoremConj(F, z, w, x₀)
    result1 = 0. + im * 0. 
    result2 = 0. + im * 0. 
    for i in eachindex(F)
        result1 += F[i] * w[i] / (z[i]-x₀)
        result1 -= conj(F[i]) * conj(w[i]) / (conj(z[i])-x₀)
        result2 += w[i] / (z[i]-x₀)
        result2 -= conj(w[i]) / (conj(z[i])-x₀)
    end
    # return result / (2*pi*im)
    return result1/result2
end

# ################################### Here is curve#################################

global curvePdown
let (x_local,w_local) = gausslegendremesh(ir, uv, parabolastep, 2) 
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
        w_out[i] = w_local[i] * (2 * im  * sqrt(uv * bigM))
    end
    global curveVup = Parameters(x_local, w_out, [- bigM + uv + 2 * im * sqrt(bigM * uv) * t  for t in x_local])
end

# global curvePdown
# let (x_local,w_local) = gausslegendremesh(ir, uv, parabolastep, 2) 
#     w_out = Vector{ComplexF64}(undef, length(w_local))
#     for i in eachindex(x_local)
#         w_out[i] = w_local[i] * ( 1 - im * sqrt(bigM / x_local[i])/2 )
#     end
#     global curvePdown = Parameters(x_local, w_out, [- bigM + t - im * sqrt(bigM * t) for t in x_local])
# end

# global curveVup
# let (x_local,w_local) = gausslegendremesh(0.0, 1.0, verticlestep, 1) 
#     w_out = Vector{ComplexF64}(undef, length(w_local))
#     for i in eachindex(x_local)
#         w_out[i] = w_local[i] * (im  * sqrt(uv * bigM))
#     end
#     global curveVup = Parameters(x_local, w_out, [- bigM + uv + im * sqrt(bigM * uv) * t  for t in x_local])
# end



curveP2 = curvePdown + curveVup

##############################
# # 这一部分测试了旧的柯西积分
# global curvePup
# let (x_local,w_local) = gausslegendremesh(ir, uv, parabolastep, 2) 
#     w_out = Vector{ComplexF64}(undef, length(w_local))
#     global curvePup = Parameters(x_local, -conj(curvePdown.w), conj(curvePdown.points))
# end

# global curveVdown
# let (x_local,w_local) = gausslegendremesh(0.0, 1.0, verticlestep, 1) 
#     w_out = Vector{ComplexF64}(undef, length(w_local))
#     global curveVdown = Parameters(x_local, -conj(curveVup.w), conj(curveVup.points))
# end

# cauchystep = 2 * verticlestep + 2 * parabolastep
# curveP2 = curvePdown + curvePup + curveVup + curveVdown
# curveConj = curvePdown + curveVup
# functest(x) =  x-3* x^2
# F = [functest(t) for t in curveP2.points];
# Fconj = [functest(t) for t in curveConj.points];
# atwhere = 9000+0.01*im
# println(cauchyTheoremConj(Fconj, curveConj.points,curveConj.w, atwhere))
# println(cauchyTheorem(F, curveP2.points,curveP2.w, atwhere))
# println(functest(atwhere))
###############################

#= 这一部分是测试那个m——hat
Nf = 4;
rm = 12/(33 - 2*Nf);
Λ = 0.234;

pp = df.pp[end-10:end]
Mass1 = df.Mass1[end-10:end]
mhat = zero(11)
mhat = (1/2 * log.(pp .^ 2 / Λ^2)) .^ rm .* Mass1

AA=Spline1D(df.pp,df.real_A1)
BB=Spline1D(df.pp,df.real_B1)
=# 

# complexAA1 = [etp_rAA1(x) + im * etp_iAA1(x) for x in curvePdown.x]
# complexBB1 = [etp_rBB1(x) + im * etp_iBB1(x) for x in curvePdown.x]
# complexAA2 = [etp_rAA2(x) + im * etp_iAA2(x) for x in curveVup.x]
# complexBB2 = [etp_rBB2(x) + im * etp_iBB2(x) for x in curveVup.x]

# complexAA = vcat(complexAA1,complexAA2)
# complexBB = vcat(complexBB1,complexBB2)

complexAA1 = df_part1[!,2] .+ im * df_part1[!,3]
complexBB1 = df_part1[!,4] .+ im * df_part1[!,5]
complexAA2 = df_part2[!,2] .+ im * df_part2[!,3] 
complexBB2 = df_part2[!,4] .+ im * df_part2[!,5]

complexAA = vcat(complexAA1,complexAA2)
complexBB = vcat(complexBB1,complexBB2)

# 测试柯西积分结果
# println(cauchyTheoremConj(complexAA, curveP2.points,curveP2.w, ir))
# println(cauchyTheoremConj(complexAA, curveP2.points,curveP2.w, 19^2))
# println(cauchyTheoremConj(complexBB, curveP2.points,curveP2.w, ir))
# println(cauchyTheoremConj(complexBB, curveP2.points,curveP2.w, 19^2))
# println(cauchyTheoremConj(complexAA, curveP2.points,curveP2.w, -0.3^2))

function createCauchyTheoremConjClosure(F, z, w)
    # 返回一个新的函数，该函数接受x₀作为输入
    return function(x₀)
        result1 = 0.0 + im * 0.0
        result2 = 0.0 + im * 0.0
        for i in eachindex(F)
            result1 += F[i] * w[i] / (z[i] - x₀)
            result1 -= conj(F[i]) * conj(w[i]) / (conj(z[i]) - x₀)
            result2 += w[i] / (z[i] - x₀)
            result2 -= conj(w[i]) / (conj(z[i]) - x₀)
        end
        # 返回计算结果
        return result1 / result2
    end
end

realA = createCauchyTheoremConjClosure(complexAA, curveP2.points,curveP2.w)
realB = createCauchyTheoremConjClosure(complexBB, curveP2.points,curveP2.w)

# 测试图像
# scatter(df_part1[!,1].^2,df_part1[!,2],xaxis=:log,yaxis=:log)
# scatter!(curvePdown.x,etp_rAA1.(curvePdown.x))


# 这里证明了辛现银取点的方式
# df_part1[!,1].^2 .- curvePdown.x
# df_part2[!,1] .- curveVup.x


function ExtraA(k,A)
    # IntA=Splines1D(k,A)
    lengthofk=length(k)
    a=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[1]][[2]]`; Adata=[k A][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    b=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[2]][[2]]`; Adata=[k A][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    return a, b
end

function ExtraB(k,B)
    # IntA=Splines1D(k,A)
    lengthofk=length(k)
    a=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[1]][[2]]`; Adata=[k B][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    b=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[2]][[2]]`; Adata=[k B][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    return a, b
end



end # module