from playFHE.math.poly import Polynomial
from playFHE.ckks.ciphertext import Ciphertext
from playFHE import util
from math import e, pi
import numpy as np


class CKKS:
    def __init__(self, N, P, q0, delta, L, sec_dist):
        self.m = 2*N
        self.xi = e ** (2 * pi * 1j / self.m)  # Used for encoding
        self.sigma_R_basis = self.create_sigma_R_basis() # used for encoding
        self.sigma_R_basis_T = util.transpose(self.sigma_R_basis) # precompute transposition

        self.P = P #Used to scale evaluation key to avoid large error terms

        self.q0 = q0 # First mod
        self.delta = delta # Scale mod
        self.L = L # Number of available levels


        self.sec_dist = sec_dist
        """Distribution used for secret key term s (all error terms e, e0, e1, v, u automatically use discrete gaussian)
        Source: https://eprint.iacr.org/2024/463.pdf """
    
    def qL(self):
        """Return the full modulus"""
        return self.q0 * self.delta**self.L


    def vandermonde(self):
        """Computes the Vandermonde matrix from an m-th root of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        n = self.m // 2
        matrix = []
        for i in range(n):
            root = self.xi ** (2 * i + 1)
            row = []
            for j in range(n):
                row.append(root ** j)
            matrix.append(row)
        return matrix

    ########

    def create_sigma_R_basis(self):
        """Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        return util.transpose(self.vandermonde())

    def compute_basis_coordinates(self, z):
        """Computes the coordinates of a vector with respect to the orthogonal lattice basis
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        output = [(util.vdot(z, b) / util.vdot(b, b)).real for b in self.sigma_R_basis]
        return output
        
    def sigma_R_discretization(self, z):
        """Projects a vector on the lattice using coordinate wise random rounding.
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        coordinates = self.compute_basis_coordinates(z)
        
        rounded_coordinates = util.coordinate_wise_random_rounding(coordinates)
        y = util.matmul(self.sigma_R_basis_T, rounded_coordinates)
        y = [yy[0] for yy in y] # "Flattens" the vector from the from [[1], [2], [3]] to [1, 2, 3]
        return y
    
    ##########
    
    def sigma_inverse(self, vec):
        """Encodes the vector b in a polynomial using an M-th root of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        van = self.vandermonde()
        #coeffs = util.gaussian_elimination(van, vec)
        coeffs = np.linalg.solve(van, vec)
        p = Polynomial(coeffs)
        return p
    
    def sigma(self, p):
        """Decodes a polynomial by applying it to the M-th roots of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        outputs = []
        n = self.m // 2
        for i in range(n):
            root = self.xi ** (2 * i + 1)
            output = p.solve(root)
            outputs.append(output)
        return outputs
    
    ###########
    
    def pi(self, z):
        """Projects a vector of H into C^{N/2}
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        n = self.m // 4
        return z[:n]
    
    def pi_inverse(self, z):
        """Expands a vector of C^{N/2} by expanding it with its complex conjugate
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        z_conjugate = z[::-1]
        z_conjugate = [x.conjugate() for x in z_conjugate]
        return z + z_conjugate
    
    
    ############
    
    
    def encode(self, vec):
        """Encodes a vector by expanding it first to H, scale it, project it on the lattice of sigma(R), and performs sigma inverse.
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""

        vec = [complex(vv, 0) for vv in vec] #Convert list of real numbers to complex numbers
        pi_z = self.pi_inverse(vec)
        scaled_pi_z = [self.delta * z for z in pi_z] #self.delta * pi_z
        rounded_scaled_pi_z = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scaled_pi_z)

        coef = [round(x.real) for x in p.coeffs]
        p = Polynomial(coef, self.qL())
        return p
        
    
    def decode(self, p):
        """Decodes a polynomial by removing the scale, evaluating on the roots, and project it on C^(N/2)
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        rescaled_p = Polynomial([c / self.delta for c in p.coeffs])
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        #Extract real parts of complex vector
        pi_z = [c.real for c in pi_z]
        return pi_z
    
    
    ######################################################


    def keygen(self):
        """Function to generate public key, private key
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        match self.sec_dist:
            case 1:
                s_c = util.sample_uniform_coeffs(self.m // 2, self.qL())
            case 2:
                s_c = util.sample_gaussian_coeffs(self.m // 2)
            case 3:
                s_c = util.sample_uniform_ternary_coeffs(self.m // 2)
            case _:
                s_c = util.sample_uniform_ternary_coeffs(self.m // 2)
        s = Polynomial(s_c, self.qL())
        a = Polynomial(util.sample_uniform_coeffs(self.m // 2, self.qL()), self.qL())
        e = Polynomial(util.sample_gaussian_coeffs(self.m // 2), self.qL())
        b = a*s
        b = b * -1
        b = b + e
        
        sk = (1, s)
        pk = (b, a)
        
        return pk, sk
    

    def evkeygen(self, sk):
        """ Function to generate evaluation key for Ciphertext-Ciphertext multiplication
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        ap = Polynomial(util.sample_uniform_coeffs(self.m // 2, self.P * self.qL()), self.P * self.qL())
        ep = Polynomial(util.sample_gaussian_coeffs(self.m // 2), self.P * self.qL())

        bp = ap * sk[1]
        bp = bp * -1
        bp = bp + ep
        p_term = (sk[1] * sk[1] * self.P)
        bp = bp + p_term
        evk = (bp, ap)
        return evk
    
    def encrypt(self, m, pk):
        """Encrypts a previously encoded plaintext into a ciphertext using the public key
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        #v = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        v = Polynomial([1]*(self.m//2), self.qL()) # This polynomial with 1s as coefficients is a placeholder until I understand the purpose of the term "v"
        e0 = Polynomial(util.sample_gaussian_coeffs(self.m // 2), self.qL())
        e1 = Polynomial(util.sample_gaussian_coeffs(self.m // 2), self.qL())
        
        #c = ((v, v) * self.pk) + (m + e0, e1)
        c0 = (v * pk[0]) + m + e0
        c1 = (v * pk[1]) + e1
        
        c = Ciphertext(c0, c1, self.P, self.q0, self.delta, self.L)
        return c
    
    def decrypt(self, c, sk):
        """Decrypts a ciphertext using the secret key and returns a plaintext polynomial
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        return c.b + (c.a * sk[1])


