module testpkg
using Pkg
Pkg.add("Pkg")
begin
    println("测试安装的包")
    # println("注意，必须在电脑中先安装Mathmatica or wolfram-kernel")
    Pkg.add("CSV")
    Pkg.add("DataFrames")
    Pkg.add("Dierckx")
    Pkg.add("FastGaussQuadrature")
    Pkg.add("JLD2")
    Pkg.add("MathLink")
    Pkg.add("Plots")
    Pkg.add("Interpolations")
    Pkg.add("ProgressMeter")
    Pkg.add(url = "https://github.com/kangjiayin/Gaussquad.jl")
    Pkg.add("Arpack")
end


end # module