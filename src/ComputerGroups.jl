module ComputerGroups

using Combinatorics
using Primes

include("structures.jl")
export Group, FactorGroup, element_order, element_power, create_from_generators

include("property_functions.jl")
export find_subgroups, find_inverse, is_normal, find_normal_subgroups, is_abelian, find_equivalences, find_factor_group, find_center, find_commutator, find_derived_series, is_solvable, find_inn, is_isomorphic, conjugacy_classes

# --- Example: Dihedral group of order 10 --- #
# d10_elems = reshape([(i,j) for i in 0:4, j in 0:1], 10)
# function d10_op(a, b)
#     if a[2] == 1 
#         ((a[1]+(5-b[1])) % 5, (a[2]+b[2]) % 2)
#     else
#         ((a[1]+b[1]) % 5, (a[2]+b[2]) % 2)
#     end
# end

# --- Example: S3 --- #
# s3_elems = [Dict(1 => i[1], 2 => i[2], 3 => i[3]) for i in permutations((1,2,3))]
# function s3_op(a, b)
#     result = Dict()
#     for i in 1:3 
#         result[i] = a[b[i]]
#     end
#     return result
# end

end # module
