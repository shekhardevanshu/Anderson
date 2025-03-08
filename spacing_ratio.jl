using Statistics, LinearAlgebra

function getSpacingRatio(eVal) #Phys. Rev. B 75, 155111

    nconv = length(eVal)
    spr_s = []
    
    for n in 2:nconv-1
        dn = eVal[n+1] - eVal[n]
        dn_1 = eVal[n] - eVal[n-1]
        
        if dn != 0 && dn_1 != 0
            push!(spr_s, ( min(dn, dn_1) / max(dn, dn_1) ))
        end

    end

    return mean(spr_s)

end
