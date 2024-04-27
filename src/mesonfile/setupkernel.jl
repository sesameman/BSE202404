# 初始化变量
if dataset["mesonBSE"]["mesonmode"] == 1
    kernel=Array{ComplexF64}(undef,kstep , zstep, kstep, zstep, 4, 4)
elseif dataset["mesonBSE"]["mesonmode"] == 2
    kernel=Array{ComplexF64}(undef,kstep , zstep, kstep, zstep, 4, 4)
elseif dataset["mesonBSE"]["mesonmode"] == 3
    kernel=Array{ComplexF64}(undef,kstep , zstep, kstep, zstep, 8, 8)
elseif dataset["mesonBSE"]["mesonmode"] == 4
    kernel=Array{ComplexF64}(undef,kstep , zstep, kstep, zstep, 8, 8)
end
# A1=Array{Float64}(undef,kstep,zstep);
# B1=Array{Float64}(undef,kstep,zstep);
# A2=Array{Float64}(undef,kstep,zstep);
# B2=Array{Float64}(undef,kstep,zstep);
# D_precompute=Array{Function}(undef, kstep, zstep, kstep, zstep);
# for j=1:kstep
# for m=1:zstep
#     A1[j,m]=AA((P2/4+meshk[j]+sqrt(0*im + P2*meshk[j])*meshz[m]))
#     B1[j,m]=BB((P2/4+meshk[j]+sqrt(0*im + P2*meshk[j])*meshz[m]))
#     A2[j,m]=AA((P2/4+meshk[j]-sqrt(0*im + P2*meshk[j])*meshz[m]))
#     B2[j,m]=BB((P2/4+meshk[j]-sqrt(0*im + P2*meshk[j])*meshz[m]))
# end
# end

# for j=1:kstep
# for m=1:zstep
# for i=1:kstep
# for s=1:zstep
# D_precompute[i,s,j,m]=y->D(meshk[i]+meshk[j]-2*sqrt(meshk[j]*meshk[i])*(meshz[s]*meshz[m] + y*sqrt(1 - meshz[s]^2)*sqrt(1 - meshz[m]^2)))
# end
# end
# end
# end

# print("完成kernel准备,用时",round((time()-timetest)*100)/100,"s \n")
