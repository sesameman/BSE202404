timetest=time()
# 压平二维数组
function flat(x::Array)
    output=Array{ComplexF64}(undef,dim)
    for i=1:kstep
        for j=1:zstep
            output[j+(i-1)*zstep]=x[i,j]
        end
    end
    return output
end

function getkernel(n1,n2)
    kernelsetup = Array{ComplexF64}(undef,dim,dim)
    for j = 1:kstep
        for m = 1:zstep
            kernelsetup[:,m+(j-1)*zstep]=flat(kernel[:,:,j,m,n1,n2])
        end
    end
    return kernelsetup
end

if dataset["mesonBSE"]["mesonmode"] == 1
kernelsolve=[getkernel(1,1) getkernel(1,2) getkernel(1,3) getkernel(1,4);
             getkernel(2,1) getkernel(2,2) getkernel(2,3) getkernel(2,4);
             getkernel(3,1) getkernel(3,2) getkernel(3,3) getkernel(3,4);
             getkernel(4,1) getkernel(4,2) getkernel(4,3) getkernel(4,4)]
elseif dataset["mesonBSE"]["mesonmode"] == 2
    kernelsolve=[getkernel(1,1) getkernel(1,2) getkernel(1,3) getkernel(1,4);
                 getkernel(2,1) getkernel(2,2) getkernel(2,3) getkernel(2,4);
                 getkernel(3,1) getkernel(3,2) getkernel(3,3) getkernel(3,4);
                 getkernel(4,1) getkernel(4,2) getkernel(4,3) getkernel(4,4)]
elseif (dataset["mesonBSE"]["mesonmode"] == 3)|(dataset["mesonBSE"]["mesonmode"] == 4)
    kernelsolve=[getkernel(1,1) getkernel(1,2) getkernel(1,3) getkernel(1,4) getkernel(1,5) getkernel(1,6) getkernel(1,7) getkernel(1,8) ;
                 getkernel(2,1) getkernel(2,2) getkernel(2,3) getkernel(2,4) getkernel(2,5) getkernel(2,6) getkernel(2,7) getkernel(2,8) ;
                 getkernel(3,1) getkernel(3,2) getkernel(3,3) getkernel(3,4) getkernel(3,5) getkernel(3,6) getkernel(3,7) getkernel(3,8) ;
                 getkernel(4,1) getkernel(4,2) getkernel(4,3) getkernel(4,4) getkernel(4,5) getkernel(4,6) getkernel(4,7) getkernel(4,8) ;
                 getkernel(5,1) getkernel(5,2) getkernel(5,3) getkernel(5,4) getkernel(5,5) getkernel(5,6) getkernel(5,7) getkernel(5,8) ;
                 getkernel(6,1) getkernel(6,2) getkernel(6,3) getkernel(6,4) getkernel(6,5) getkernel(6,6) getkernel(6,7) getkernel(6,8) ;
                 getkernel(7,1) getkernel(7,2) getkernel(7,3) getkernel(7,4) getkernel(7,5) getkernel(7,6) getkernel(7,7) getkernel(7,8) ;
                 getkernel(8,1) getkernel(8,2) getkernel(8,3) getkernel(8,4) getkernel(8,5) getkernel(8,6) getkernel(8,7) getkernel(8,8) ]
end

println("now Pmass = ",sqrt(-P2))
print("eigs = ")
println(eigs(kernelsolve, nev=2, which=:LM)[1])

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