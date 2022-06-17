# GGA

def fitness_vector(v):
    """
    Return the fitness of the vector ``v``
    
    INPUT:
        - ``v`` -- a binary vector
    """
    f = 0
    
    for i in v:
        if i != 0:
            f += 1
    
    return f

def fitness_matrix(M):
    """
    Return the minimum fitness of the rows of ``M``
    
    INPUT:
        - ``M`` -- a binary matrix
    """
    min_f = Infinity
    
    for row in M:
        f = fitness_vector(row)
        
        if f < min_f:
            min_f = f
    
    return min_f

def permutation_matrix(x):
    """
    Return the matrix of permutations associated with ``x``
    
    INPUT:
        - ``x`` -- a permutation
    """
    M = matrix(GF(2), len(x), len(x))
    
    for j in range(len(x)):
        i = x[j] - 1
        M[i,j] = 1
    
    return M

def fitness_permutation(M, x):
    """
    Return the the minimum fitness of the rows of the 
    reduced row echelon form of the matrix ``M`` 
    permuted by ``x``
    
    INPUT:
        - ``M`` -- a binary matrix
        - ``x`` -- a permutation applied to ``M``
    """
    P = permutation_matrix(x)
    
    return fitness_matrix((M * P).rref())

def random_sol(n):
    """
    Return a random vector of length ``n``
    
    INPUT:
        - ``n`` -- length
    """
    sol = []
    
    while len(sol) != n:
        r = randint(1, n)
        
        if r not in sol:
            sol = sol + [r]
            
    return sol

def crossover(p1, p2):
    """
    Return the descendant of the cross of ``p1`` and ``p2``
    
    INPUT:
        - ``p1`` -- a vector of length ``n``
        - ``p2`` -- a vector of length ``n``
    """
    return Permutation(p1) * Permutation(p2), Permutation(p2) * Permutation(p1)

def mutation(p):
    """
    Return a mutated vector
    
    INPUT:
        - ``p`` -- a vector
    """
    n = len(p)
    k = randint(1, n - 1)
    k1 = randint(0, k - 1)
    k2 = randint(k, n - 1)
    
    p[k1], p[k2] = p[k2], p[k1]

    return p

def GGA(G, N, pc, max_reinit, min_fitness):
    """
    Return a permutation of the ``G`` whose reduced row echelon form
    has a row with weight ``min_fitness``
    
    INPUT:
        - ``G`` -- a matrix made up of the public key of a McEliece 
          cryptosystem and a row with the encrypted message
        - ``N`` -- population size
        - ``pc`` -- crossover probability
        - ``max_reinit`` -- maximum number of evaluations of the 
          solutions that do not produce improvement in the fitness 
          of the best solution found.
        - ``min_fitness`` -- minimun fitness
    """
    t = 0
    max_t = 500000
    Pt = []
    reinit = 0
    n = G.ncols()
    
    # Initialize the Population P(t)
    for i in range(0, N):
        sol = random_sol(n)
        Pt = Pt + [sol]
    
    Pt = Matrix(Pt)
    
    # Evaluate
    fitness_Pt = []
    for Pti in Pt:
        fitness_Pt = fitness_Pt + [fitness_permutation(G, Pti)]
        
    while(t < max_t and min(fitness_Pt) > min_fitness):
        #print(min(fitness_Pt))
        Pt_sig = []
        fitness_Pt_sig = []
        
        # Binary tournament selection
        parents = []
        for i in range(0, N):            
            ind1 = Pt[randint(0, N-1)]
            ind2 = Pt[randint(0, N-1)]
            
            if (fitness_permutation(G, ind1) >= fitness_permutation(G, ind2)):
                parents = parents + [ind2]
            else:
                parents = parents + [ind1]
        
        parents = Matrix(parents)
        
        # Generation
        n_crossovers = pc * N/2 # number of crossovers

        for i in range(0, floor(N/2)):            
            if n_crossovers > 0:
                c1, c2 = crossover(parents[2*i], parents[2*i + 1])
                n_crossovers -= 1
            else:
                c1 = mutation(copy(parents[2*i]))
                c2 = mutation(copy(parents[2*i + 1]))
            
            # Update
            Pt_sig = Pt_sig + [c1]
            Pt_sig = Pt_sig + [c2]
        
        if N % 2 != 0:
            Pt_sig = Pt_sig + [parents[N-1]]
        
        Pt_sig = Matrix(Pt_sig)
        
        # Evaluate
        for Pt_sigi in Pt_sig:
            fitness_Pt_sig = fitness_Pt_sig + [fitness_permutation(G, Pt_sigi)]
        
        # No improvement
        if min(fitness_Pt_sig) > min(fitness_Pt):
            worst = fitness_Pt_sig.index(max(fitness_Pt_sig))
            best = fitness_Pt.index(min(fitness_Pt))
            Pt_sig[worst] = Pt[best]
            fitness_Pt_sig[worst] = fitness_Pt[best]
            reinit += 1
        else:
            reinit = 0
        
        # Restart
        if reinit >= max_reinit:
            reinit = 0
            Pt_sig = []
            for i in range(0, N-1):
                Pt_sig = Pt_sig + [random_sol(n)]
            
            best = fitness_Pt.index(min(fitness_Pt))
            Pt_sig = Pt_sig + [Pt[best]]
            Pt_sig = Matrix(Pt_sig)
            
            # Evaluate
            fitness_Pt_sig = []
            for Pt_sigi in Pt_sig:
                fitness_Pt_sig = fitness_Pt_sig + [fitness_permutation(G, Pt_sigi)]
        
        fitness_Pt = fitness_Pt_sig
        Pt = Pt_sig
        t += 1
    
    best = fitness_Pt.index(min(fitness_Pt))
    
    return Pt[best], fitness_permutation(G, Pt[best])


# CHC

def distance(x, y):
    """
    Return de Hamming distance between ``x`` and ``y``
    
    INPUT:
        - ``x`` -- a vector of length ``n``
        - ``y`` -- a vector of length ``n``
    """
    assert(len(x) == len(y))
    
    n = len(x)
    count = 0
    for i in range(0, n):
        if x[i] != y[i]:
            count += 1
            
    return count

def actualizar(P, tau):
    """
    Update the crossover decrement
    
    INPUT:
        - ``P`` -- population
        - ``tau`` -- crossover rate
    """
    d = 0
    max_dist = 0
    count = 0
    
    for i in range(0, P.nrows()):
        for j in range(i, P.nrows()):
            dist = distance(P[i], P[j])
            
            if dist > max_dist:
                max_dist = dist
            
            d += dist
            count += 1
    
    d = d * 1.0/count
    dec = tau * max_dist
    
    return d, dec

def CHC(G, N, tau, min_fitness):
    """
    Return a permutation of the ``G`` whose reduced row echelon form
    has a row with weight ``min_fitness``
    
    INPUT:
        - ``G`` -- a matrix made up of the public key of a McEliece 
          cryptosystem and a row with the encrypted message
        - ``N`` -- population size
        - ``tau`` -- crossover rating
        - ``min_fitness`` -- minimun fitness
    """
    t = 0
    max_t = 500000
    Pt = []
    reinit = 0
    n = G.ncols()
    
    # Initialize the Population P(t)
    for i in range(0, N):
        sol = random_sol(n)    
        Pt = Pt + [sol]
    
    Pt = Matrix(Pt)
    
    # Evaluate
    fitness_Pt = []   
    for Pti in Pt:
        fitness_Pt = fitness_Pt + [fitness_permutation(G, Pti)]
    
    # Distance
    d, dec = actualizar(Pt, tau)
        
    while(t < max_t and min(fitness_Pt) > min_fitness):
        #print(min(fitness_Pt))
        Ct = []
        Pt_sig = []
        fitness_Ct = []
        fitness_Pt_sig = []
        improvement = False
        
        # Parent selection
        parents = []
        for i in range(0, N):            
            p = randint(0, N-1)
            parents = parents + [Pt[p]]
        
        parents = Matrix(parents)
        
        # Generation
        for i in range(0, N/2):            
            if distance(parents[2*i], parents[2*i + 1]) < d:
                c1, c2 = crossover(parents[2*i], parents[2*i + 1])
                Ct = Ct + [c1]
                Ct = Ct + [c2]
        
        Ct = Matrix(Ct)
        
        # Evaluate
        for Cti in Ct:
            fitness_Ct = fitness_Ct + [fitness_permutation(G, Cti)]
        
        if Ct.nrows() == 0:
            Pt_sig = Pt
            fitness_Pt_sig = fitness_Pt
        else:
            fitness_Pt_aux = fitness_Pt
            fitness_Ct_aux = fitness_Ct
            
            for i in range(0, N):
                best_Pt = min(fitness_Pt_aux)
                best_Ct = min(fitness_Ct_aux)

                if best_Pt < best_Ct:
                    index = fitness_Pt_aux.index(best_Pt)
                    Pt_sig = Pt_sig + [Pt[index]]
                    fitness_Pt_sig = fitness_Pt_sig + [best_Pt]
                    fitness_Pt_aux[index] = Infinity
                else:
                    index = fitness_Ct_aux.index(best_Ct)
                    Pt_sig = Pt_sig + [Ct[index]]
                    fitness_Pt_sig = fitness_Pt_sig + [best_Ct]
                    fitness_Ct_aux[index] = Infinity
                    improvement = True
        
            Pt_sig = Matrix(Pt_sig)
        
        # No improvement
        if not improvement:
            d = d - dec
            
            if d <= 0:
                Pt_sig = []
                for i in range(0, N-1):
                    Pt_sig = Pt_sig + [random_sol(n)]

                best = fitness_Pt.index(min(fitness_Pt))
                Pt_sig = Pt_sig + [Pt[best]]
                Pt_sig = Matrix(Pt_sig)
                
                # Distance
                d, dec = actualizar(Pt, tau)
        
        fitness_Pt = fitness_Pt_sig
        Pt = Pt_sig
        t += 1
    
    best = fitness_Pt.index(min(fitness_Pt))
    
    return Pt[best], fitness_permutation(G, Pt[best])