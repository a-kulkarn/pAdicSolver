push!(LOAD_PATH, pwd()*"/src/")
using Revise
using HeckeExt

## Some test data.
Qp,F7 = test_rings()
A = random_test_matrix(Qp)

