
function get_IPR(U, nev, N)
    ipr_pgi = Vector{Float64}(undef, nev)

    for n in 1:nev
        ipr_pgi[n] = sum((U[(n-1)*N+1:n*N]).^4)
    end

    ipr_pgi
end