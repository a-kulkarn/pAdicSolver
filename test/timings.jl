###############################################################################################
#
#    Timings
#
###############################################################################################
#
# In a timing run, we keep track of the average time to solve a system, as well as the
# average number of solutions of a polynomial system.
#
# (At least, it seems like a good idea to do this.)

using pAdicSolver, Random

# Enforce that tests are identical.
Random.seed!(123)


##############################
# Affine plane curve timings
##############################

Qp = PadicField(31,10)
R, (x1,x2) = PolynomialRing(Qp, 2)
num_samp = 10

failed_test_count = 0
for d in [2,4,6,10,15]
    tot_time = 0
    tot_sols = 0
    
    for i=1:num_samp
        local P = random_square_system(R, d)

        t0 = time()
        local sol = solve_affine_system(P)
        t1 = time()
        tot_time += t1-t0
        tot_sols += size(sol, 1)
    end

    avg_time = tot_time/num_samp
    avg_sols = tot_sols/num_samp

    @info "Affine: ambient dimension = 2, polynomial degrees = $d" avg_time avg_sols
end


##############################
# Projective quadrics timings.
##############################

for n in [2,4,5]
    R, vars = PolynomialRing(Qp, n+1)
    tot_time  = 0
    tot_sols = 0
    
    for i=1:num_samp
        local P = random_square_system(R, 2, homogeneous=true)

        t0 = time()
        local sol = solve_projective_system(P)
        t1 = time()
        tot_time += t1-t0
        tot_sols += size(sol, 1)
    end

    avg_time = tot_time/num_samp
    avg_sols = tot_sols/num_samp
    @info "Projective, ambient dimension = $n, polynomial degrees = 2" avg_time avg_sols
end
