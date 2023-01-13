


module ComplexMultiplicationFields

using Oscar
using Combinatorics
export is_cmfield, is_cmtype, all_cmtypes

# Write your package code here.
# Bogdan Adrian Dina ...
# todo: 1. write documentation about the algorithms
# 2. write examples for the methods
# 3. write also the algorithms numerically over CC.



function is_cmfield(K::NumField)
	degK = degree(K)
	
	if is_totally_real(K)
		return false
	elseif ( degK % 2 ) != 0
		return false
	end
	
	subfs = subfields(K)
	realsubfs = [ tup for tup in subfs if ( is_totally_real(tup[1]) == true ) &&  degree(tup[1]) == ( degK//2 ) ]
	
	if ( length(realsubfs) > 1 ) || ( length(realsubfs) == 0 )
		return false
	end
	
	return true
end

function totally_real_subfield(K::NumField)
	
	if is_cmfield(K) == false
		return []
	end
	
	degK = degree(K) 
	subfs = subfields(K)
	realsubfs = [ tup for tup in subfs if ( is_totally_real(tup[1]) == true ) &&  degree(tup[1]) == ( degK//2 ) ]
	
	tup = []
	if length(realsubfs) == 1
		tup = realsubfs[1]
	end
		
	return tup
end

function is_cmtype(phi)

	L = codomain( phi[1] )
	ccL = complex_conjugation(L)
	lphi = length(phi)
    for i in 1:lphi
        for j in i+1:lphi
            if ( phi[i]*ccL == phi[j] )
                return false
			end
		end
	end
	return true
end

function all_cmtypes(K::NumField)
	
	if is_cmfield(K) == false
	 	return false, []
	end	
	L = splitting_field(K.pol)
	_, phi = is_subfield(K, L)
	AutL = automorphism_list(L)
	morphisms_K_to_L = [ phi*tau for tau in AutL ]
	
	degK = degree(K) 
	tups = collect( combinations( morphisms_K_to_L , Integer(degK//2) ) )
	all_cmtypes = [ phi for phi in tups if is_cmtype(phi) == true ]

	return all_cmtypes	
end

function all_cmtypes_numerically(K::NumField)
# TODO:	
end




end
