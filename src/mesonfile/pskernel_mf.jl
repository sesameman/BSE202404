# timetest = time()
# 计算
# i & s 代表外动量 k2
# j & m 代表内动量 q2

# 想到一个好方法,  循环做这个outIndex, 用反函数求出j, s,能扁平化计算


#==========================================setup_kernel==========================================#
if dataset["calcType"]["mesonmode"] == 1
    kernel = Array{Array{ComplexF64,2},2}(undef, 4, 4)
    for i in 1:4
        for j in 1:4
            kernel[i, j] = fill(Complex(0), dim, dim)
        end
    end
elseif dataset["calcType"]["mesonmode"] == 2
    kernel = Array{Array{ComplexF64,2},2}(undef, 4, 4)
    for i in 1:4
        for j in 1:4
            kernel[i, j] = fill(Complex(0), dim, dim)
        end
    end
elseif dataset["calcType"]["mesonmode"] == 3
    kernel = Array{Array{ComplexF64,2},2}(undef, 8, 8)
    for i in 1:8
        for j in 1:8
            kernel[i, j] = fill(Complex(0), dim, dim)
        end
    end
elseif dataset["calcType"]["mesonmode"] == 4
    kernel = Array{Array{ComplexF64,2},2}(undef, 8, 8)
    for i in 1:8
        for j in 1:8
            kernel[i, j] = fill(Complex(0), dim, dim)
        end
    end
end
#================================================================================================#

#==========================================calc_kernel==========================================#
Threads.@threads for i in eachindex(meshk)
    for s1 in eachindex(meshz1)
        for s2 in eachindex(meshz2)
            for j in eachindex(meshk)
                for m1 in eachindex(meshz1)
                    for m2 = eachindex(meshz2)
                        k2 = meshk[i]
                        q2 = meshk[j]
                        zk1 = meshz1[s1]
                        zk2 = meshz2[s2]
                        zq1 = meshz1[m1]
                        zq2 = meshz2[m2]
                        vec_p3 = safe_sqrt(externalP2)
                        vec_p4 = im * safe_sqrt(externalP2 - P2)
                        vec_k2 = sqrt(1 - zk2^2) * sqrt(1 - zk1^2) * sqrt(k2)
                        vec_k3 = zk2 * sqrt(1 - zk1^2) * sqrt(k2)
                        vec_k4 = zk1 * sqrt(k2)
                        vec_q3 = zq2 * sqrt(1 - zq1^2) * sqrt(q2)
                        vec_q4 = zq1 * sqrt(q2)
                        pk = vec_p3*vec_k3 + vec_p4*vec_k4
                        pq = vec_p3*vec_q3 + vec_p4*vec_q4
                        qplus2 = P2/4 + q2 + pq
                        qminus2 = P2/4 + q2 - pq
                        pqplus = pq + P2/2 
                        pqminus = pq - P2/2
                        qqplus = q2 + pq/2
                        qqminus = q2 - pq/2
                        weightzq1 = weightz1[m1]
                        weightzq2 = weightz2[m2]
                        weightq2 = weightk[j]
                        Aqplus = AA(qplus2)
                        Aqminus = AA(qminus2)
                        Bqplus = BB(qplus2)
                        Bqminus = BB(qminus2)
                        outIndex = computeIndex(i, s1, s2)
                        interIndex = computeIndex(j, m1, m2)
                        allweight = - z2^2 * weightzq1 * weightzq2 * weightq2 * q2 / (32 * pi^4) / (qplus2*Aqplus^2 + Bqplus^2) / (qminus2*Aqminus^2 + Bqminus^2)
                        # allweight = -weightzq*weightq2*q2/(16*pi^3)/((P2/4+q2+safe_sqrt(P2*q2)*zq)*Aqplus^2+Bqplus^2)/((P2/4+q2-safe_sqrt(P2*q2)*zq)*Aqminus^2+Bqminus^2)*z2^2
                        theta, weighttheta = gausslegendremesh(0.0, 2.0 *pi, zintstep, 1)
                        # theta, weighttheta = gausschebyshev(zintstep, 1)
                        for numtheta in eachindex(theta)
                            vec_q2 = sqrt(1 - zq2^2) * sqrt(1 - zq1^2) * sqrt(q2) * cos(theta[numtheta]) 
                            # vec_q2 = sqrt(1 - zq2^2) * sqrt(1 - zq1^2) * theta[numtheta] * sqrt(q2)
                            kq = vec_q2 * vec_k2 + vec_q3 * vec_k3 + vec_q4 * vec_k4
                            kqplus = kq + pk/2
                            kqminus = kq - pk/2
                            kernel[1, 1][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( -4*(Bqminus*Bqplus + Aqminus*Aqplus*(-0.25*P2 + q2)) )
                            
                            kernel[1, 2][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( 4*Aqminus*Bqplus*pqminus - 4*Aqplus*Bqminus*pqplus )
                            
                            kernel[1, 3][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( 4*pq*(Aqminus*Bqplus*qqminus - Aqplus*Bqminus*qqplus) )
                            
                            kernel[1, 4][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( 4*Aqminus*Aqplus*(pqplus*qqminus - pqminus*qqplus) )
                            
                            kernel[2, 1][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*(k2^2*(-(Aqminus*Bqplus*pqminus) + Aqplus*Bqminus*pqplus) + pk*((Aqminus*Bqplus*kqminus - Aqplus*Bqminus*kqplus)*q2 + kq*(2*Aqminus*Bqplus*(-2*kqminus + qqminus) + 2*Aqplus*Bqminus*(2*kqplus - qqplus))) + k2*((Aqminus*Bqplus*kqminus - Aqplus*Bqminus*kqplus)*pk + 2*Aqminus*Bqplus*kqminus*pq - 2*Aqplus*Bqminus*kqplus*pq - 2*Aqplus*Bqminus*kq*pqplus + Aqminus*Bqplus*pqminus*(2*kq - q2) + Aqplus*Bqminus*pqplus*q2 - 2*Aqminus*Bqplus*pq*qqminus + 2*Aqplus*Bqminus*pq*qqplus)))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[2, 2][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*(k2^2*(-2*Aqminus*Aqplus*pqminus*pqplus + P2*(-(Bqminus*Bqplus) + Aqminus*Aqplus*(-0.25*P2 + q2))) + k2*(2*Bqminus*Bqplus*kq*P2 - 2*Bqminus*Bqplus*pq^2 + 2*Aqminus*Aqplus*kqplus*pq*pqminus + 2*Aqminus*Aqplus*kqminus*pq*pqplus + 4*Aqminus*Aqplus*kq*pqminus*pqplus - Bqminus*Bqplus*P2*q2 - 2*Aqminus*Aqplus*pqminus*pqplus*q2 + Aqminus*Aqplus*(-0.25*P2 + q2)*(2*pq^2 + P2*(-2*kq + q2)) + pk^2*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + pk*(Aqminus*Aqplus*(kqplus*pqminus + kqminus*pqplus) + 2*pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2))) - 2*Aqminus*Aqplus*pq*pqplus*qqminus - 2*Aqminus*Aqplus*pq*pqminus*qqplus) + pk*(Aqminus*Aqplus*(kqplus*pqminus + kqminus*pqplus)*q2 - pk*(4*kq - q2)*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + 2*kq*(pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*(-2*kqminus + qqminus) + pqminus*(-2*kqplus + qqplus))))))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[2, 3][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*pq*(-(k2^2*(pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*qqminus + pqminus*qqplus))) + k2*(-3*Bqminus*Bqplus*pq*q2 + 3*Aqminus*Aqplus*pq*q2*(-0.25*P2 + q2) + 2*Aqminus*Aqplus*kqplus*pq*qqminus - Aqminus*Aqplus*pqplus*q2*qqminus + 2*Aqminus*Aqplus*kqminus*pq*qqplus - Aqminus*Aqplus*pqminus*q2*qqplus - 4*Aqminus*Aqplus*pq*qqminus*qqplus + pk*(kq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(kqplus*qqminus + kqminus*qqplus)) + 2*kq*(2*pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*qqminus + pqminus*qqplus))) + pk*(kq^2*(-4*Bqminus*Bqplus + 4*Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*q2*(kqplus*qqminus + kqminus*qqplus) + kq*(3*q2*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) - 4*Aqminus*Aqplus*(kqplus*qqminus + (kqminus - qqminus)*qqplus)))))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[2, 4][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*(k2^2*(-(pq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus)) + P2*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus)) + pk*(-4*kq^2*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) - pk*q2*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus) + kq*(3*Aqminus*Bqplus*pqminus*q2 + 3*Aqplus*Bqminus*pqplus*q2 + 2*(2*pk - pq)*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus))) + k2*(-3*Aqminus*Bqplus*pq*pqminus*q2 - 3*Aqplus*Bqminus*pq*pqplus*q2 + 2*Aqminus*Bqplus*pq^2*qqminus + Aqminus*Bqplus*P2*q2*qqminus + Aqplus*Bqminus*(2*pq^2 + P2*q2)*qqplus - pk^2*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus) + kq*(4*pq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) - 2*P2*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus)) + pk*(kq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) - 2*pq*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus)))))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[3, 1][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*(2*pk^2*(Aqminus*Bqplus*(kqminus - qqminus) + Aqplus*Bqminus*(-kqplus + qqplus)) + pk*(-2*Aqminus*Bqplus*kqminus*pq + 2*Aqplus*Bqminus*kqplus*pq + Aqminus*Bqplus*pqminus*(k2 - 2*kq + q2) + 2*Aqminus*Bqplus*pq*qqminus + Aqplus*Bqminus*(-(pqplus*(k2 - 2*kq + q2)) - 2*pq*qqplus)) + P2*((-(Aqminus*Bqplus*kqminus) + Aqplus*Bqminus*kqplus)*q2 + k2*(Aqminus*Bqplus*(-3*kqminus + 2*qqminus) + Aqplus*Bqminus*(3*kqplus - 2*qqplus)) + 2*kq*(Aqminus*Bqplus*(2*kqminus - qqminus) + Aqplus*Bqminus*(-2*kqplus + qqplus)))))/(3. *pk*(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[3, 2][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (8*pk^3*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + 8*pk^2*(pq*(-2*Bqminus*Bqplus + 2*Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*(kqminus - qqminus) + pqminus*(kqplus - qqplus))) + 8*pk*(Bqminus*Bqplus*pq^2 - Aqminus*Aqplus*kqplus*pq*pqminus - Aqminus*Aqplus*kqminus*pq*pqplus + Aqminus*Aqplus*pqminus*pqplus*q2 - Aqminus*Aqplus*pq^2*(-0.25*P2 + q2) + kq*(-2*Aqminus*Aqplus*pqminus*pqplus + P2*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2))) + k2*(Aqminus*Aqplus*pqminus*pqplus + P2*(-(Bqminus*Bqplus) + Aqminus*Aqplus*(-0.25*P2 + q2))) + Aqminus*Aqplus*pq*pqplus*qqminus + Aqminus*Aqplus*pq*pqminus*qqplus) + 4*P2*(-(Aqminus*Aqplus*(kqplus*pqminus + kqminus*pqplus)*q2) - 2*kq*(pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*(-2*kqminus + qqminus) + pqminus*(-2*kqplus + qqplus))) + k2*(2*pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*(-3*kqminus + 2*qqminus) + pqminus*(-3*kqplus + 2*qqplus)))))/(3. *pk*(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[3, 3][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (-4*pq*(-2*pk^2*((kq - q2)*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(kqplus*qqminus + (kqminus - 2*qqminus)*qqplus)) + pk*(-3*Bqminus*Bqplus*pq*q2 + 3*Aqminus*Aqplus*pq*q2*(-0.25*P2 + q2) + 2*Aqminus*Aqplus*kqplus*pq*qqminus - Aqminus*Aqplus*pqplus*q2*qqminus + 2*Aqminus*Aqplus*kqminus*pq*qqplus - Aqminus*Aqplus*pqminus*q2*qqplus - 4*Aqminus*Aqplus*pq*qqminus*qqplus - k2*(pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*qqminus + pqminus*qqplus)) + 2*kq*(2*pq*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(pqplus*qqminus + pqminus*qqplus))) + P2*(kq^2*(-4*Bqminus*Bqplus + 4*Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*q2*(kqplus*qqminus + kqminus*qqplus) + k2*((3*kq - 2*q2)*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) + Aqminus*Aqplus*(3*kqplus*qqminus + (3*kqminus - 4*qqminus)*qqplus)) + kq*(3*q2*(Bqminus*Bqplus - Aqminus*Aqplus*(-0.25*P2 + q2)) - 4*Aqminus*Aqplus*(kqplus*qqminus + (kqminus - qqminus)*qqplus)))))/(3. *pk*(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[3, 4][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*(-2*pk^3*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus) + 2*pk^2*(kq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) - Aqminus*Bqplus*pqminus*q2 + 2*Aqminus*Bqplus*pq*qqminus + Aqplus*Bqminus*(-(pqplus*q2) + 2*pq*qqplus)) + pk*(-2*kq*(2*pq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) + P2*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus)) + k2*(pq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) + 2*P2*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus)) + pq*(3*Aqminus*Bqplus*pqminus*q2 + 3*Aqplus*Bqminus*pqplus*q2 - 2*pq*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus))) + P2*(-(k2*(3*kq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) - 2*Aqminus*Bqplus*pqminus*q2 + 2*Aqminus*Bqplus*pq*qqminus + 2*Aqplus*Bqminus*(-(pqplus*q2) + pq*qqplus))) + kq*(4*kq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus) - 3*Aqminus*Bqplus*pqminus*q2 + 2*Aqminus*Bqplus*pq*qqminus + Aqplus*Bqminus*(-3*pqplus*q2 + 2*pq*qqplus)))))/(3. *pk*(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) )
                            
                            kernel[4, 1][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*Aqminus*Aqplus*(pqplus*(kqminus*q2 - 2*kq*qqminus) - kqplus*(pqminus*q2 + 2*(pk - pq)*qqminus) + k2*(-(pqplus*(kqminus - 2*qqminus)) + pqminus*(kqplus - 2*qqplus)) + 2*(kqminus*(pk - pq) + kq*pqminus)*qqplus))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) ) + Parameters_Xi * D_infrared(k2 + q2 - 2kq) * ( (32*Aqminus*Aqplus*(-(kqplus*pqminus) + kqminus*pqplus))/(9*k2*P2 - 9*pk^2) )
                            
                            kernel[4, 2][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*(2*Aqminus*Bqplus*kqminus*pq^2 + 2*Aqplus*Bqminus*kqplus*pq^2 - 2*Aqminus*Bqplus*kq*pq*pqminus - 2*Aqplus*Bqminus*kq*pq*pqplus - Aqminus*Bqplus*kqminus*P2*q2 - Aqplus*Bqminus*kqplus*P2*q2 + 2*Aqminus*Bqplus*kq*P2*qqminus + 2*Aqplus*Bqminus*kq*P2*qqplus + 2*pk^2*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus) + k2*(Aqminus*Bqplus*kqminus*P2 + Aqplus*Bqminus*kqplus*P2 - Aqminus*Bqplus*(pk - 2*pq)*pqminus - 2*Aqminus*Bqplus*P2*qqminus - Aqplus*Bqminus*((pk - 2*pq)*pqplus + 2*P2*qqplus)) + pk*(-2*Aqminus*Bqplus*kqminus*pq - 2*Aqplus*Bqminus*kqplus*pq + Aqminus*Bqplus*pqminus*q2 - 2*Aqminus*Bqplus*pq*qqminus + Aqplus*Bqminus*(pqplus*q2 - 2*pq*qqplus))))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) ) + Parameters_Xi * D_infrared(k2 + q2 - 2kq) * ( (-32*(Aqminus*Bqplus*(kqminus*P2 - pk*pqminus) + Aqplus*Bqminus*(kqplus*P2 - pk*pqplus)))/(9*k2*P2 - 9*pk^2) )
                            
                            kernel[4, 3][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*pq*((Aqminus*Bqplus*kqminus*pq + Aqplus*Bqminus*kqplus*pq - kq*(Aqminus*Bqplus*pqminus + Aqplus*Bqminus*pqplus))*q2 + 2*pk*(-(Aqminus*Bqplus*kqminus*q2) - Aqplus*Bqminus*kqplus*q2 + kq*(Aqminus*Bqplus*qqminus + Aqplus*Bqminus*qqplus)) + k2*(Aqminus*Bqplus*kqminus*pq + Aqplus*Bqminus*kqplus*pq - Aqminus*Bqplus*pqminus*(kq - 2*q2) - 2*Aqminus*Bqplus*pq*qqminus - Aqplus*Bqminus*(pqplus*(kq - 2*q2) + 2*pq*qqplus))))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) ) + Parameters_Xi * D_infrared(k2 + q2 - 2kq) * ( (32*pq*(Aqminus*Bqplus*(-(kqminus*pq) + kq*pqminus) + Aqplus*Bqminus*(-(kqplus*pq) + kq*pqplus)))/(9*k2*P2 - 9*pk^2) )
                            
                            kernel[4, 4][outIndex, interIndex] += weighttheta[numtheta] * D(k2 + q2 - 2 * kq) * ( (4*(Bqminus*Bqplus*kq*P2*q2 + Aqminus*Aqplus*kqplus*pq*pqminus*q2 + Aqminus*Aqplus*kqminus*pq*pqplus*q2 - 2*Aqminus*Aqplus*kq*pqminus*pqplus*q2 + Aqminus*Aqplus*kq*P2*q2*(-0.25*P2 + q2) - 2*Aqminus*Aqplus*kqplus*pq^2*qqminus + 2*Aqminus*Aqplus*kq*pq*pqplus*qqminus + Aqminus*Aqplus*kqplus*P2*q2*qqminus - 2*Aqminus*Aqplus*kqminus*pq^2*qqplus + 2*Aqminus*Aqplus*kq*pq*pqminus*qqplus + Aqminus*Aqplus*kqminus*P2*q2*qqplus - 4*Aqminus*Aqplus*kq*P2*qqminus*qqplus + 2*pk^2*(q2*(Bqminus*Bqplus + Aqminus*Aqplus*(-0.25*P2 + q2)) - 2*Aqminus*Aqplus*qqminus*qqplus) + k2*(2*Bqminus*Bqplus*pq^2 + Aqminus*Aqplus*kqplus*pq*pqminus + Aqminus*Aqplus*kqminus*pq*pqplus - 2*Bqminus*Bqplus*P2*q2 + 4*Aqminus*Aqplus*pqminus*pqplus*q2 + 2*Aqminus*Aqplus*(-0.25*P2 + q2)*(pq^2 - P2*q2) + kq*(-2*Aqminus*Aqplus*pqminus*pqplus + P2*(Bqminus*Bqplus + Aqminus*Aqplus*(-0.25*P2 + q2))) - Aqminus*Aqplus*kqplus*P2*qqminus - 4*Aqminus*Aqplus*pq*pqplus*qqminus - Aqminus*Aqplus*kqminus*P2*qqplus - 4*Aqminus*Aqplus*pq*pqminus*qqplus + 4*Aqminus*Aqplus*P2*qqminus*qqplus + pk*(-(pq*(Bqminus*Bqplus + Aqminus*Aqplus*(-0.25*P2 + q2))) + Aqminus*Aqplus*(pqplus*qqminus + pqminus*qqplus))) + pk*(-(Aqminus*Aqplus*q2*(pqplus*(2*kqminus + qqminus) + pqminus*(2*kqplus + qqplus))) + 2*kq*(-(pq*(Bqminus*Bqplus + Aqminus*Aqplus*(-0.25*P2 + q2))) + Aqminus*Aqplus*(pqplus*qqminus + pqminus*qqplus)) + pq*(-(q2*(Bqminus*Bqplus + Aqminus*Aqplus*(-0.25*P2 + q2))) + 2*Aqminus*Aqplus*(kqplus*qqminus + (kqminus + 2*qqminus)*qqplus)))))/(3. *(-(k2*P2) + pk^2)*(k2 - 2*kq + q2)) ) + Parameters_Xi * D_infrared(k2 + q2 - 2kq) * ( (8*(Bqminus*Bqplus*(-4*kq*P2 + 4*pk*pq) + Aqminus*Aqplus*(kq*(P2^2 + 8*pqminus*pqplus - 4*P2*q2) + P2*(-(pk*pq) + 4*kqplus*qqminus + 4*kqminus*qqplus) - 4*(kqplus*pq*pqminus + kqminus*pq*pqplus - pk*pq*q2 + pk*pqplus*qqminus + pk*pqminus*qqplus))))/(9*k2*P2 - 9*pk^2) )
                        end # end for Integral kernel
        
                        for i in 1:size(kernel, 1)
                            for j in 1:size(kernel, 2)
                                kernel[i, j][outIndex, interIndex] *= allweight
                            end
                        end #end for add all weight
        
                    end
                end
            end
        end
    end
end
#===============================================================================================#

#==========================================solve_kernel==========================================#
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
#=================================================================================================#


println("now Pmass = ",sqrt(-P2))
println("eigs = ",eigvals)