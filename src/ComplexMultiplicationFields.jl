


module ComplexMultiplicationFields

using Oscar
using Combinatorics
export is_cmfield, is_cmtype, all_cmtypes

# Write your package code here.
# Bogdan Adrian Dina ...
# todo: 1. write documentation about the algorithms
# 2. write examples for the methods
# 3. write also the algorithms numerically over CC.


"""

"""
function compute_involution(K::NumField)

	return true
end


"""
	is_cmfield(K::NumField) --> ...
	
	Documentation: 
	Tests if K is a CM field. If so, also returns the totally real subfield K0 as well as the corresponding involution.
	
	
"""
function is_cmfield(K::NumField)
	
	if is_totally_real(K)
		return false, id_hom(K)
	elseif ( degree(K) % 2 ) != 0
		return false, id_hom(K)
	end
	
	subfs = subfields(K)
	realsubfs = [ tup for tup in subfs if ( is_totally_real(tup[1]) == true ) &&  degree(tup[1]) == ( degree(K)//2 ) ]
	
	if ( length(realsubfs) > 1 ) || ( length(realsubfs) == 0 )
		return false, id_hom(K)
	end
	
	#TODO: find element corresponding to complex conjugation on K and return it
	# Write for that a function. 
	#auts = automorphism_list(K)
	ccK = compute_involution(K)
	
	# computes the totally real sub
	K0_tup = totally_real_subfield(K)
	
	return true, ccK, K0_tup
end

"""
	totally_real_subfield(K::NumField) --> Tuple{AnticNumberField, NfToNfMor}
	
	Documentation: 
	Given a cm field K, returns its totally real subfield K0 and an embedding iota: K0 --> K.
	If K is not a cm field, returns `false` and the identity map on K.
"""
function totally_real_subfield(K::NumField)
	
	if is_cmfield(K) == false
		return false, id_hom(K)
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

"""
	is_cmtype(...) --> Bool
	Documentation:
	...
"""
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


"""
	all_cmtypes(K::NumField, primitive = true) --> 
	Documentation:
	...
"""
function all_cmtypes(K::NumField, primitive = true)
	
	@assert is_cmfield(K) "error: K is not a cm field ..."
	
	L = splitting_field(K.pol)
	_, phi = is_subfield(K, L)
	AutL = automorphism_list(L)
	morphisms_K_to_L = [ phi*tau for tau in AutL ]
	
	degK = degree(K) 
	tups = collect( combinations( morphisms_K_to_L , Integer(degK//2) ) )
	all_cmtypes = [ phi for phi in tups if is_cmtype(phi) == true ]

	return all_cmtypes	
end



end
