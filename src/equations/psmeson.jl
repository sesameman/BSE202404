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
    # 计算kernel, 求解bse
    include(joinpath(workdir,"src/mesonfile/pskernel.jl"))
    # 保存文件
    savefilepath = joinpath(workdir,"data/pseudo_BSE/meson-$now_time/")
    jldsave(joinpath(savefilepath,"p=$psave.jld2");P2, eigvals, eigvecs)
    open(joinpath(savefilepath,"output.txt"), "a") do file
        println(file, "now Pmass = ",sqrt(-P2))
        println(file, "eigs = ",eigvals)
    end
    # print("mesonBSA--$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
    if indexforp2 == 1
        global eig_results
        eig_results = EigenResult[]
    end
    
    push!(eig_results, EigenResult(P2, eigvals[1], eigvecs[:,1]))
    
    
    if Pstep != 1
        if indexforp2 == Pstep
            global normalizenum
            normalizenum = - real(eig_results[1].eigenvalues - eig_results[end].eigenvalues)/(eig_results[1].p - eig_results[end].p)
            for renum in eachindex(eig_results)
                eig_results[renum].eigenfunctions = eig_results[renum].eigenfunctions ./ normalizenum
            end
            savefilepath = joinpath(workdir,"data/pseudo_BSE/meson-$now_time/")
            @save joinpath(savefilepath,"normalizeed.jld2") eig_results
        end
    end

end