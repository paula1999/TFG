from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
from sage.coding.linear_code import AbstractLinearCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder

def get_defining_set(n, pol, field):
    """
    Return a tuple of ``n`` distinct elements of ``field``
    that are not roots of ``pol`` 
    
    INPUT:
        - ``n`` -- length
        - ``pol`` -- a monic polynomial with coefficients in `field`
        - ``field`` -- finite field
    """
    defining_set = []
    
    while len(defining_set) != n:
        aux = field.random_element()
        if pol(aux) != 0 and aux not in defining_set:
            defining_set = defining_set + [aux]
    
    return defining_set

def random_word(n, field):
    """
    Return a random vector of length ``n`` belonging to the 
    finite field ``field``
    
    INPUT:
        - ``n`` -- length
        - ``field``-- finite field
    """
    v = []

    while len(v) != n:
        v = v + [choice(field.list())]

    word = vector(field, v)
    
    return word

def random_error(n, num_errors, field):
    """
    Return a random vector of length ``n`` with maximum weight 
    ``num_errors`` belonging to the finite field ``field``
    
    INPUT:
        - ``n`` -- length
        - ``num_errors`` -- maximum weight
        - ``field`` -- finite field
    """
    v = [0] * n
    i = 0
    
    while i < num_errors:
        m = randint(0, len(v)-1)
        v[m] = choice(field.list())
        i += 1

    e = vector(field, v)
    
    return e

class Goppa(AbstractLinearCode):
    r"""
    Implementation of Goppa codes.
    
    INPUT:
    - ``field`` -- finite field on which `self` is defined.
    - ``generating_pol`` -- a monic polynomial with coefficients in
    a finite field `\GF{p^m}` extending from `field`
    - ``defining_set`` -- tuple of `n` distinct elements of `\GF{p^m}`
    that are not roots of `generating_pol`
    """
    def __init__(self, defining_set, generating_pol, field):
        """
        Initialize.
        """
        
        if not generating_pol.is_monic():
            raise ValueError("ERROR. Generating polynomial isn't monic")
        
        for gamma in defining_set:
            if generating_pol(gamma) == 0:
                raise ValueError("ERROR. Defining elements are roots of generating polynomial")
        
        self._field = field
        self._field_L = generating_pol.base_ring()

        if not self._field.is_field() or not self._field.is_finite():
            raise ValueError("ERROR. Generating polynomial isn't definied over a finite field")
        
        self._length = len(defining_set)
        self._generating_pol = generating_pol
        self._defining_set = defining_set
        
        super(Goppa, self).__init__(self._field, self._length, "GoppaEncoder", "GoppaDecoder")
        
        if self.get_dimension() == 0:
            raise ValueError("ERROR. Code dimension is null")
    
    def _repr_(self):
        """
        Representation of a Goppa code.
        """
        return "[{}, {}] Goppa code".format(self._length, self.get_dimension())
    
    def get_generating_pol(self):
        """
        Return the generating polynomial of ``self``.
        """ 
        return self._generating_pol
    
    def get_defining_set(self):
        """
        Return the defining set of ``self``.
        """ 
        return self._defining_set
    
    def get_parity_pol(self):
        """
        Return the parity polynomial of ``self``.
        """
        parity_pol = list()
        
        for elem in self._defining_set:
            parity_pol.append((self._generating_pol.parent().gen() - elem).inverse_mod(self._generating_pol))
        
        return parity_pol
    
    def get_parity_check_matrix(self):
        """
        Return a parity check matrix of ``self``.
        """
        V, from_V, to_V = self._field_L.vector_space(self._field, map = True)

        parity = self.get_parity_pol()

        vector_L = []
        for i in range(len(parity)):
            vector_L = vector_L + parity[i].list()

        vector_F = []
        for i in range(len(vector_L)):
            vector_F = vector_F + to_V(vector_L[i]).list()

        matriz_F = matrix(self._length, vector_F)
        matriz_F = matriz_F.T
        
        return matriz_F
    
    def get_generator_matrix(self):
        """
        Return a generador matrix of ``self``.
        """
        H = self.get_parity_check_matrix()
        G = transpose(H).left_kernel().basis_matrix()
        
        return G
    
    def get_dimension(self):
        """
        Return the dimension of the code.
        """
        
        return rank(self.get_generator_matrix())
    
class GoppaEncoder(Encoder):
    r"""
    Encoder for Goppa codes
    
    INPUT:
    - ``code`` -- code associated with the encoder
    """
    def __init__(self, code):
        """
        Initialize.
        """
        super(GoppaEncoder, self).__init__(code)
        self._generator_matrix = self.code().get_generator_matrix()
        
    def _repr_(self):
        """
        Representation of a encoder for a Goppa code
        """
        return "Encoder for {}".format(self.code())
    
    def get_generator_matrix(self):
        """
        Return a generador matrix of the code
        """
        return self._generator_matrix
    
    def encode (self, m):
        """
        Return a codeword
        
        INPUT:
        - ``m`` -- a vector to encode
        """
        return m * self.get_generator_matrix()
    


class GoppaDecoder(Decoder):
    r"""
    Decoder for Goppa codes
    
    INPUT:
    - ``code`` -- code associated with the decoder
    """
    
    def __init__(self, code):
        """
        Initialize
        """
        super(GoppaDecoder, self).__init__(code, code.ambient_space(), "GoppaDecoder")
                
        self._generating_pol = self.code().get_generating_pol()
        self._defining_set = self.code().get_defining_set()
    
    def _repr_(self):
        """
        Representation of a decoder for a Goppa code
        """
        return "Decoder for {}".format(self.code())
    
    def get_syndrome(self, c):
        """
        Return the syndrome polynomial
        
        INPUT:
        - ``c`` -- a element of the input space of ``self``
        """   
        h = self.code().get_parity_pol()
        
        syndrome = 0
        
        for i in range(len(h)):
            syndrome = syndrome + c[i]*h[i]
            
        return syndrome
    
    def get_generating_pol(self):
        """
        Return the generating polynomial
        """
        return self._generating_pol
    
    def decode_to_code(self, word):
        r"""
        Corrects the errors in ``word`` and returns a codeword
        
        INPUT:
        - ``word`` -- a codeword of ``self``
        """
        i = 0

        # Step 1
        S = self.get_syndrome(word)
        
        if S == 0:
            return word

        # Step 2
        r_prev = self.get_generating_pol()
        t = floor(self.get_generating_pol().degree()/2)
        r_i = S
        U_prev = 0
        U_i = 1
        
        # Steps 3 and 4
        while r_i.degree() >= t or r_prev.degree() < t:
            (q, r) = r_prev.quo_rem(r_i)
            aux_r_i = r_i
            aux_U_i = U_i
            r_i = r
            U_i = q * U_i + U_prev
            r_prev = aux_r_i
            U_prev = aux_U_i  
            i += 1
        
        # Step 5
        # make sigma monic     
        sigma = U_i / U_i.coefficients()[-1]
        
        eta = (-1)^i * r_i / U_i.coefficients()[-1]
        
        # roots of sigma are the locations of the errors
        roots_loc = []
        for root in sigma.roots():
            if g(root[0] != 0):
                roots_loc = roots_loc + [self._defining_set.index(root[0])]
        
        error = [0] * len(self._defining_set)
        x = self.get_generating_pol().parent().gen()
        sigma_diff = sigma.diff(x)
        
        for i in range(0, len(roots_loc)):
            error[roots_loc[i]] = eta.subs(x = sigma.roots()[i][0]) / (sigma_diff.subs(x = sigma.roots()[i][0]))
        
        x = word - vector(self.code()._field, error)
        
        return x

    def decode_to_message(self, word):
        r"""
        Decode ``word`` to the message space.
        
        INPUT:
        - ``word`` -- a codeword of ``self``
        """
        word = self.decode_to_code(word)
        message = self.code().get_generator_matrix().solve_left(word)
        return message

Goppa._registered_encoders["GoppaEncoder"] = GoppaEncoder
Goppa._registered_decoders["GoppaDecoder"] = GoppaDecoder