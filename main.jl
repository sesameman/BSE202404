workdir = @__DIR__
# include(joinpath(@__DIR__,"module/testpkg.jl")) # 检测包完整性
# include(joinpath(@__DIR__,"src/equations/quark_gap.jl"))
# include(joinpath(@__DIR__,"module/inter.jl"))
include(joinpath(workdir,"module/readdse.jl"))
include(joinpath(workdir,"src/meson.jl"))
