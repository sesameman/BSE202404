mkpath("data/pseudo_BSE/meson-$now_time/")
print("创建文件--",joinpath(workdir,"data/pseudo_BSE/meson-$now_time/"),"\n")

# cp(joinpath(workdir,"config.toml"),joinpath(workdir,"data/pseudo_BSE/meson-//log-$now_time.toml"))
open(joinpath(workdir,"data/pseudo_BSE/meson-$now_time/log-$now_time.toml"), "w") do file
    TOML.print(file, dataset)
end

CSV.write(joinpath(workdir,"data/pseudo_BSE/meson-$now_time/plist.CSV"),DataFrame(p2=plist))

@showprogress for indexforp2=1:Pstep
    timetest1=time()
    global P2
    P2 = -plist[indexforp2]^2
    local psave
    psave = plist[indexforp2]
    # 分配点与权重
    include(joinpath(workdir,"src/mesonfile/setupkernel.jl"))
    # 计算kernel
    include(joinpath(workdir,"src/mesonfile/pskernel.jl"))
    # 求解函数
    include(joinpath(workdir,"src/mesonfile/solve_kernel.jl"))
    # 保存文件
    savefilepath = joinpath(workdir,"data/pseudo_BSE/meson-$now_time/")
    jldsave(joinpath(savefilepath,"p=$psave.jld2");P2, eigvals, eigvecs)
    open(joinpath(savefilepath,"output.txt"), "a") do file
        println(file,"now Pmass = ",psave)
        println(file,"eigs = ",eigvals)
    end
    # print("mesonBSA--$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
end