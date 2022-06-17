def get_weight(c):
    """
    Return the weight of ``c``
    
    INPUT:
        - ``c`` -- a vector
    """
    w = 0
    
    for elem in c:
        if elem != 0:
            w += 1
    
    return w

class McEliece:
    r"""
    Implementation of McEliece cryptosystem.
    
    INPUT:
    - ``n`` -- length
    - ``q`` -- Dimension of the finite field `GF(2^q)`
    - ``g``-- a monic polynomial with coefficients in
    a finite field `\GF{2^q}`
    """
    def __init__(self, n, q, g):
        """
        Initialize.
        """
        # Goppa code
        F = GF(2)
        L = GF(2^q)
        defining_set = get_defining_set(n, g, L)
        C = Goppa(defining_set, g, F)
        G = C.get_generator_matrix()
        k = G.nrows()
        n = G.ncols()
    
        # random binary nonsingular matrix
        S = matrix(F, k, [choice(F.list()) for i in range(k^2)])
        
        while rank(S) < k:
            i = randint(0, k-1)
            j = randint(0, k-1)
            S[i,j] = choice(F.list())
        
        # random permutation matrix
        columns = list(range(n))
        P = matrix(F, n)
        
        for i in range(n):
            l = randint(0, len(columns)-1)
            j = columns[l]
            P[i,j] = 1
            columns.remove(j)
        
        self._goppa = C
        self._G = G
        self._t = floor(g.degree()/2)
        self._k = k
        self._n = n
        self._S = S
        self._P = P
        self._public_key = S * self._G * P
        
    def __repr__(self):
        """
        Representation of a McEliece cryptosystem.
        """
        return "McEliece cryptosystem over {}".format(self._goppa)
        
    def get_S(self):
        """
        Return the random binary nonsingular matrix that is part of the private key
        """
        return self._S
    
    def get_G(self):
        """
        Return the generating matrix associated with the Goppa code that is part of the private key
        """
        return self._G
    
    def get_P(self):
        """
        Return the random permutation matrix that is part of the private key
        """
        return self._P
    
    def get_public_key(self):
        """
        Return the public key
        """
        return self._public_key

    def encrypt(self, m):
        """
        Return a chipertext
        
        INPUT:
        - ``m`` -- a plaintext to encrypt
        """
        # random errors
        e = vector([F(0) for i in range(self._n)])
        
        while get_weight(e) != self._t:
            i = randint(0, len(e) - 1)
            e[i] = choice(F.list())
        
        # chipertext
        c = m * self._public_key + e
        
        return c
    
    def decrypt(self, c):
        """
        Return the plain text associated with the cryptogram ``c``
        
        INPUT:
        - ``c`` -- a chiphertext
        """
        word = c * self._P^(-1)
        D = GoppaDecoder(self._goppa)
        word = D.decode_to_message(word)
        message = word * (self._S)^(-1)
        
        return message