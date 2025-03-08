using Printf
using DelimitedFiles, Statistics
using FLoops
using LinearAlgebra

BLAS.set_num_threads(1)

include("spee.jl")
include("spacing_ratio.jl")
include("papb.jl")
include("ipr.jl")

function binaryRead(fname)
    data = Array{Float64}(undef, filesize(fname)÷8)
    read!(fname, data)
    return data
end


function main(args)
    
    k = parse(Int, args[1])
    t = parse(Float64, args[2])
    w1 = parse(Float64, args[3])

    L=12; nconv=20; s=1000
    # L=18; nconv=50; s=100

    # barr = sort(vcat(collect(0.5:0.5:10), collect(11:1:25), collect(30:10:100), collect(120:20:200)))
    barr = unique(sort(vcat(collect(0.1:0.3:2), collect(2:0.3:5), collect(0.5:0.5:10), collect(11:1:25), collect(30:10:100), collect(120:20:200))))
    # barr = sort(vcat(collect(0.1:0.5:10), collect(10:1:19), collect(20:20:100), collect(150:50:500)));
    bname = map(x -> @sprintf("%0.1f", x), barr)
    t = @sprintf("%0.1f", t)
    w1 = @sprintf("%0.1f", w1)
    e = "0.0"

    N = L^3
    Na = div(N, 2)
    sites = collect(1:Na)
    # sites = get_sites_subsystem(L)
    
    vn_s = Array{Float64}(undef, s*nconv)
    vn_h = Matrix{Float64}(undef, length(barr), length(vn_s))

    pab_s = Array{Float64}(undef, s*nconv)
    pab_h = Matrix{Float64}(undef, length(barr), length(vn_s))

    spr_h = Matrix{Float64}(undef, length(barr), s)
    lsp_h = Matrix{Float64}(undef, length(barr), s*(nconv-1))
    
    spr_s = zeros(s)
    lsp_s = zeros(s*(nconv-1))

    ipr_s = Array{Float64}(undef, s*nconv)
    ipr_g = Matrix{Float64}(undef, length(barr), length(ipr_s))
    

    for (j, b) in enumerate(bname)

        @floop for m in 1:s

            eVec = binaryRead("./dataBinAE_PBC_L=$L/eigvec_L=$(L)_w=$(b)_w1=$(w1)_t=$(t)_k=$(k)_itr_$(m)_e=$e.dat")
            eVal = binaryRead("./dataBinAE_PBC_L=$(L)/eigval_L=$(L)_w=$(b)_w1=$(w1)_t=$(t)_k=$(k)_itr_$(m)_e=$e.dat")
            eVal2 = sort(eVal[1:nconv])

            vn_s[(m-1)*nconv+1:nconv*m] = get_sp_entropy(eVec, sites, nconv, N)
            pab_s[(m-1)*nconv+1:nconv*m] = get_papb(eVec, sites, nconv, N)
            lsp_s[(nconv-1)*(m-1)+1:(nconv-1)*m] = diff(eVal2)
	        spr_s[m] = getSpacingRatio(eVal2)
            ipr_s[(m-1)*nconv+1:nconv*m] = get_IPR(eVec, nconv, N)
          
        end

        vn_h[j,:] = vn_s
        pab_h[j,:] = pab_s
        spr_h[j, :] = spr_s
        lsp_h[j,:] = lsp_s
        ipr_g[j,:] = ipr_s

	    print("completed w = $b\n")
        flush(stdout)
                
    end

    writedlm("vn_ae_pbc,L=$L,nev=$nconv,k=$k,t=$t,w1=$w1,e=$e.txt", [vn_h barr])
    writedlm("lsp_ae_pbc,nev=$nconv,L=$L,k=$k,w1=$w1,t=$t,e=$e.txt", [lsp_h barr])
    writedlm("spr_ae_pbc,nev=$nconv,L=$L,k=$k,w1=$w1,t=$t,e=$e.txt", [spr_h barr])
    writedlm("ipr_ae_pbc,L=$L,nev=$nconv,k=$k,w1=$w1,t=$t,e=$e.txt", [ipr_g barr])
    writedlm("pa_pb_ae_pbc,L=$L,nev=$nconv,k=$k,t=$t,w1=$w1,e=$e.txt", [pab_h barr])

end

main(ARGS)
