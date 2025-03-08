using Printf
using DelimitedFiles, Statistics
# using FLoops
using LinearAlgebra

function binaryRead(fname)
    data = Array{Float64}(undef, filesize(fname)÷8)
    read!(fname, data)
    return data
end

function main(args)

    garr = vcat(collect(0.0:0.01:0.1), collect(0.15:0.05:0.5), collect(0.6:0.1:1)) 
    gname = map(x -> @sprintf("%0.3f", x), garr)
    
    l = parse(Int, args[1])
    N = 2^l

    nconv = parse(Int, args[2])
    s = parse(Int, args[3])
    e = parse(Float64, args[4])
    e = @sprintf("%0.1f", e)

    z_s = Array{Float64}(undef, s*nconv)
    z_g = Matrix{Float64}(undef, length(garr), length(z_s))
    z_pgi = Vector{Float64}(undef, nconv)

    for (j, g) in enumerate(gname)

        for m in 1:s
            
            eVec = binaryRead("./dataBinQREM_L=$(l)/eigvec_L=$(l)_g=$(g)_itr_$(m)_e=$e.dat")
            
           for n in 1:nconv
                psi = eVec[(n-1)*N+1:n*N]
                z_pgi[n] = -sum([abs(psi[k])^2 * log(abs(psi[k])^2) for k in 1:N])
            end 

            z_s[(m-1)*nconv+1:nconv*m] = z_pgi

        end

	print("completed g = $g z = $(mean(z_s))\n")
        flush(stdout)

        z_g[j,:] = z_s
        
    end

    writedlm("state_entropy,L=$l,nev=$nconv,e=$e.txt", z_g)

end

main(ARGS)
