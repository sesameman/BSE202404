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

