using LinearAlgebra

function get_sp_entropy(U, sites, nconv, N)
    vn_phi = Vector{Float64}(undef, nconv)
    
    for n in 1:nconv

        psi = U[(n-1)*N+1:n*N]
        p_a = sum((psi[sites]).^2)
        p_b = 1-p_a
        
        if isapprox(p_a, 1)
            vn_phi[n] = -p_a*log2(p_a)
        else
            vn_phi[n] = -p_a*log2(p_a) - p_b*log2(p_b)
        end
	
    end

    vn_phi
end

# function get_sites_subsystem(L)
# 
#     N = L^3  # Total number of sites
#     d = L÷4
#     d_2 = L÷2
#     sites = []
# 
#     function index(x, y, z)
#         return (x - 1) + (y - 1) * L + (z - 1) * L * L + 1
#     end
# 
#     for x in d+1:d+d_2
#         for y in d+1:d+d_2
#             for z in d+1:d+d_2
#                 i = index(x, y, z)
#                 push!(sites, i)
#             end
#         end
#     end
# 
#     return sort(sites)
# end
