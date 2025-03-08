
function get_papb(U, sites, nconv, N)
    vn_phi = Vector{Float64}(undef, nconv)
    
    for n in 1:nconv

        psi = U[(n-1)*N+1:n*N]
	    p_a = sum((psi[sites]).^2)

        vn_phi[n] = p_a*(1-p_a)
		
    end

    vn_phi
end

