#### main.sage
load('NTRU_lattices.sage')
load('GA_auxiliary_functions.sage')

B = ntru_lattice(n_local = 10, df_local = 5, dg_local = 4, dphi_local = 12, q_local = 256)
n = B.nrows()
#B = B.LLL()
B_inv = B**-1


def GA_SVP_NTRU(B):
    print("starting new try")
    n = B.nrows()
    B_inv = B**-1
    popsize = 2*n
    pop = matrix(ZZ, popsize, n, 0)
    for i in range(popsize):
        pop[i] = vector(random.choice((-1,0,1), size=n, replace=True))
    cost_values = vector(QQ, (0 for _ in range(popsize)))
    for i in range(popsize):
        cost_values[i] = cost(pop[i], B_inv)
    # order pop by increasing cost_values and sort cost_value
    pop = matrix([val for (_, val) in sorted(zip(cost_values, pop), key=lambda x: x[0])])
    cost_values = vector(sorted(cost_values))
    if cost_values[0] == 0:
        return pop[0]
    inv_cost = vector( (1/x for x in cost_values) )
    
    iterations_without_improvement = 0
    smallest_cost = cost_values[0]
    
    random_slice = 0 #popsize // 2
    
    # initialize matrix that will contain each newly generated population
    new_pop = matrix(ZZ, popsize, n, 0)
    while iterations_without_improvement < 2*n:
        
        inv_cost = vector( (1/x for x in cost_values) )
        for q in range(0 ,popsize - 1 - random_slice):
            # select pair proportional to fitness
            #i, j = random.choice(range(popsize), size=2, p=inv_cost/sum(inv_cost), replace=False)
            # cross, mutation and localsearch
            #new_pop[q] = localsearch(mutation(crossover(pop[i], pop[j])), B_inv)
            # select pair proportional to fitness
            #i, j = random.choice(range(popsize), size=2, p=inv_cost/sum(inv_cost), replace=False)
            j = random.choice(range(popsize), size =1, p=inv_cost/sum(inv_cost), replace=False)
            new_pop[q] = localsearch(mutation(copy(pop[j])), B_inv)
            #new_pop[q] = localsearch(mutation(crossover(pop[i], pop[j])), B_inv)
        
        for q in range(popsize - 1 - random_slice, popsize - 1):
            new_pop[q] = localsearch(vector(random.choice((-1,0,1), size=n, replace=True)), B_inv)
        
        new_pop[popsize-1] = pop[0]
        pop = copy(new_pop)
        # order pop by increasing cost_values and sort cost_values
        for i in range(popsize):
            cost_values[i] = cost(pop[i], B_inv)
        pop = matrix([val for (_, val) in sorted(zip(cost_values, pop), key=lambda x: x[0])])
        cost_values = vector(sorted(cost_values))
        if cost_values[0] == 0:
            return pop[0]
        
        new_smallest_cost = cost_values[0]
        if new_smallest_cost == smallest_cost:
            iterations_without_improvement += 1
            print('_', iterations_without_improvement, end='')
        else:
            print('\n', pop[0], cost_values[0], 'advancing took ', 1+iterations_without_improvement, ' iterations.')
            iterations_without_improvement = 0
        
        smallest_cost = new_smallest_cost
        
    return vector(ZZ, (0 for _ in range(popsize)))
    

### testing examples    
import time
start = time.time()
print('n = ', n)
v = GA_SVP_NTRU(B)
while min(v) == max(v):
    v = GA_SVP_NTRU(B)
    
print('\nSVP vector =', v)
print('SVP coordinates =', v*B_inv)
end = time.time(); print('seconds passed: ', round(end - start, 2))
