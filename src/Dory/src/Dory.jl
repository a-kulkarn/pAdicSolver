module Dory


## NOTE: It is possible to resolve the impending namespace conflict for the names of linear algebra functions
# via conditional importing (see https://github.com/mikeInnes/Requires.jl). For now, we put the burden of use on the user if they want to use both packages at the same time.

## In general, the issue is a Julia developer desicion. See https://discourse.julialang.org/t/function-name-conflict-adl-function-merging/10335/30

# Functions with the same name as LinearAlgebra functions.
#
#import LinearAlgebra.eigvecs
#import LinearAlgebra.svd
#
#

export broadcast, iterate, collect, matrix, rectangular_solve, my_nullspace, eigen, eigvecs, eigspaces, charpoly, MyEigen

export /, valuation, abs, modp, test_rings, rand, rand_padic_int, random_test_matrix, padic_qr, inverse_iteration, iseigenvector, singular_values

using Hecke, Distributed
import Hecke: Generic.Mat, Generic.MatElem, nmod_mat, charpoly


## Export the namespace of Hecke for use after Dory

exclude = [ :AbstractAlgebra, :CoerceMap, :CoerceMap, :Hecke, :Nemo, :RealField, :ResidueRingPolyMap,
            :_hnf_modular, :call, :can_solve, :cols, :den, :factor, :factors, :inverse, :isid, :num,
            :parseint, :qq, :random_SMatSLP, :rows, :strongequal, :upper_triangular, :wedderburn_decomposition,
            :window, :xgcd, :zz ]



for i in names(Hecke)
    i in exclude && continue
    eval(Meta.parse("import Hecke." * string(i)))
    eval(Expr(:export, i))
end


## End Export ##

include("dory_matrix.jl")
include("padic_util.jl")

function __init__()

    if myid() == 1
        println("")
        print("You've found \n")
        printstyled("DORY", color = :red)
        println()
        print("Version")
        printstyled(" $VERSION_NUMBER ", color = :green)
        print("... \n ... which comes with absolutely no warranty whatsoever")
        println()
        println("(c) 2019 by Avinash Kulkarni")
        println()
    else
        println("DORY $VERSION_NUMBER ...")
    end

end

#####
# Version Number
####
global const VERSION_NUMBER = "0.0"




end
