#=============================================================#
#=============================================================#
# 计算大矩阵的维度
big_dim = size(kernel, 1) * dim
# 初始化大矩阵
kernelsolve = Array{ComplexF64,2}(undef, big_dim, big_dim)
# 填充大矩阵
for i in 1:size(kernel, 1)
    for j in 1:size(kernel, 2)
        # 计算当前小矩阵在大矩阵中的起始位置
        row_start = (i-1)*dim + 1
        col_start = (j-1)*dim + 1
        # 将小矩阵赋值到大矩阵的相应位置
        kernelsolve[row_start:row_start+dim-1, col_start:col_start+dim-1] = kernel[i, j]
    end
end



# kernelsolve = big_matrix
eigvals, eigvecs = eigs(kernelsolve, nev=2, which=:LM) 
println("now Pmass = ",sqrt(-P2))
println("eigs = ",eigvals)

# kernelsolve=I-kernelsolve
# solution=kernelsolve\right


# # 切比雪夫
# F1=zeros(kstep)
# F2=zeros(kstep)
# F3=zeros(kstep)
# F4=zeros(kstep)

# for u=1:4*dim
#     if u<=dim
#         u1=u
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F1[nk]+=solution[u]*weightz[nz]*2/pi
#     elseif u<=2*dim
#         u1=u-dim
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F2[nk]+=solution[u]*weightz[nz]*2/pi
#     elseif u<=3*dim
#         u1=u-2*dim
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F3[nk]+=solution[u]*weightz[nz]*2/pi
#     elseif u<=4*dim
#         u1=u-3*dim
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F4[nk]+=solution[u]*weightz[nz]*2/pi
#     end
# end

# 切比雪夫
# F1=zeros(kstep,zstep)
# F2=zeros(kstep,zstep)
# F3=zeros(kstep,zstep)
# F4=zeros(kstep,zstep)
# F5=zeros(kstep,zstep)
# F6=zeros(kstep,zstep)
# F7=zeros(kstep,zstep)
# F8=zeros(kstep,zstep)

# if (dataset["mesonBSE"]["mesonmode"] == 1) | (dataset["mesonBSE"]["mesonmode"] == 2)
#     for u=1:4*dim
#         if u<=dim
#             u1=u
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F1[nk,nz]=solution[u]
#         elseif u<=2*dim
#             u1=u-dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F2[nk,nz]=solution[u]
#         elseif u<=3*dim
#             u1=u-2*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F3[nk,nz]=solution[u]
#         elseif u<=4*dim
#             u1=u-3*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F4[nk,nz]=solution[u]
#         end
#     end
# elseif (dataset["mesonBSE"]["mesonmode"] == 3)|(dataset["mesonBSE"]["mesonmode"] == 4)
#     for u=1:8*dim
#         if u<=dim
#             u1=u
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F1[nk,nz]=solution[u]
#         elseif u<=2*dim
#             u1=u-dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F2[nk,nz]=solution[u]
#         elseif u<=3*dim
#             u1=u-2*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F3[nk,nz]=solution[u]
#         elseif u<=4*dim
#             u1=u-3*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F4[nk,nz]=solution[u]
#         elseif u<=5*dim
#             u1=u-4*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F5[nk,nz]=solution[u]
#         elseif u<=6*dim
#             u1=u-5*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F6[nk,nz]=solution[u]
#         elseif u<=7*dim
#             u1=u-6*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F7[nk,nz]=solution[u]
#         elseif u<=8*dim
#             u1=u-7*dim
#             nk=((u1-1)÷zstep+1)
#             nz=((u1-1)%zstep+1)
#             F8[nk,nz]=solution[u]
#         end
#     end
# end

# print("完成kernel求解，用时",round((time()-timetest)*100)/100,"s \n")
# print("预计误差", round(abs(B1[1,1]-F1[1]*0.003)/B1[1,1]*100000)/1000,"%\n")