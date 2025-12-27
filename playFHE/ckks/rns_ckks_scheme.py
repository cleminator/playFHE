from __future__ import annotations

from playFHE.math.poly import Polynomial
from playFHE.math.rns_poly import RNSPolynomial
from playFHE.ckks.rns_ciphertext import RNSCiphertext
from playFHE.ckks.ckks_scheme import CKKS
from playFHE import util
from playFHE.util import Number
from math import e, pi, prod
import numpy as np

class RNSKey:
    def __init__(self, b: RNSPolynomial | int, a: RNSPolynomial):
        self.a = a
        self.b = b

    def __getitem__(self, index):
        if index == 0:
            return self.b
        elif index == 1:
            return self.a
        else:
            raise IndexError("Index out of range")

    def __len__(self):
        return 2


class RNSPublicKey(RNSKey):
    def __init__(self, b: RNSPolynomial, a: RNSPolynomial):
        super().__init__(b, a)

class RNSPrivateKey(RNSKey):
    def __init__(self, s: RNSPolynomial):
        super().__init__(1, s)

class RNSMultKey(RNSKey):
    def __init__(self, b: RNSPolynomial, a: RNSPolynomial):
        super().__init__(b, a)


class RNSCKKS(CKKS):

    B: list[int]
    C: list[int]
    m: int
    xi: complex
    sigma_R_basis: list[list[Number]]
    sigma_R_basis_T: list[list[Number]]
    k: int
    p: int
    B: list[int]
    L: int
    C: list[int]
    q: int
    q0: int
    roots: dict[int, int]
    sec_dist: int


    def __init__(self, N, p, k, q0, q, L, sec_dist):#C, q, sec_dist):
        self.B = []
        self.C = []
        self.m = 2*N
        self.xi = e ** (2 * pi * 1j / self.m)  # Used for encoding
        #print("Create sigma_R_basis")
        self.sigma_R_basis = self.create_sigma_R_basis()  # used for encoding
        self.sigma_R_basis_T = util.transpose(self.sigma_R_basis)  # precompute transposition

        self.k = k
        self.p = 2**p
        #print("Create B")
        self.B = self.generate_basis(p, k)

        self.L = L
        #print("Create C")
        self.C = self.generate_basis(q0, 1)
        self.C += self.generate_basis(q, L) #self.generate_base_C(q0, q, L)
        self.q = 2**q # All moduli q_i should be as close to this q as possible, to reduce the approximation error during rescaling
        self.q0 = 2**q0
        # q will be used for the scaling during encoding and decoding

        self.roots = {} #Will contain a map of modulus->root for each q and p
        #print("Create roots")
        self.generate_roots()
        self.sec_dist = sec_dist

    def qL(self) -> int:
        """Return the full modulus"""
        return prod(self.C)

    def generate_roots(self):
        for q in self.C:
            self.roots[q] = util.find_2nth_root_of_unity(self.m // 2, q)
        for p in self.B:
            self.roots[p] = util.find_2nth_root_of_unity(self.m // 2, p)

    def generate_basis(self, bitsize: int, length: int) -> list[int]:
        b = []
        next_prime = 2**bitsize
        for i in range(length):
            while True:
                next_prime = util.find_next_prime(next_prime)
                if util.mod(next_prime, self.m) == 1 and next_prime not in self.B and next_prime not in self.C:
                    b.append(next_prime)
                    next_prime += 1
                    break
                else:
                    next_prime += 1
        return b


    def sigma_inverse(self, vec: list[Number]) -> list[Number]:
        """Encodes the vector b in a polynomial using an M-th root of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        van = self.vandermonde()
        coeffs = util.gaussian_elimination(van, vec)
        #p = Polynomial(coeffs)
        #p = RNSPolynomial(self.B, self.C, [x.real for x in coeffs])
        return [c.real for c in coeffs]


    def encode(self, vec: list[float]) -> RNSPolynomial:
        """Encodes a vector by expanding it first to H, scale it, project it on the lattice of sigma(R), and performs sigma inverse.
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        vec = [complex(vv, 0) for vv in vec]  # Convert list of real numbers to complex numbers
        pi_z = self.pi_inverse(vec)
        scaled_pi_z = [self.q * z for z in pi_z]  # self.delta * pi_z
        rounded_scaled_pi_z = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scaled_pi_z)

        coef = [round(x) for x in p]
        #p = Polynomial(coef, self.qL())

        p = RNSPolynomial(self.B, self.C, self.roots, False, coef)

        return p

    def decode(self, p: RNSPolynomial) -> list[float]:
        """Decodes a polynomial by removing the scale, evaluating on the roots, and project it on C^(N/2)
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        #rescaled_p = Polynomial([c / self.delta for c in p.coeffs])
        pcopy = RNSPolynomial(p.B[:], p.C[:], p.roots, p.ntt_domain)
        pcopy.set_limbs(p.limbs[:])
        if pcopy.ntt_domain:
            pcopy.convert_NTT_to_RNS()
        rescaled_p = Polynomial([c / self.q for c in pcopy.get_coeffs()])
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        # Extract real parts of complex vector
        pi_z = [c.real for c in pi_z]
        return pi_z

    ###########################

    def keygen(self) -> tuple[RNSPublicKey, RNSPrivateKey]:
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
        #s = RNSPolynomial(self.B, self.C, self.roots, coeffs=s_c)s
        #a = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_uniform_coeffs(self.m // 2, self.qL()))
        #print("Sampling")
        s = RNSPolynomial(self.B, self.C, self.roots, coeffs=[1]*(self.m // 2))
        a = RNSPolynomial(self.B, self.C, self.roots, coeffs=[1]*(self.m // 2))
        e = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_gaussian_coeffs(self.m // 2))
        #print("Mult poly")
        b = a * s
        #print("mult const")
        b = b * -1
        b = b + e

        #sk = (1, s)
        #pk = (b, a)

        sk = RNSPrivateKey(s)
        pk = RNSPublicKey(b, a)

        return pk, sk

    def encrypt(self, m: RNSPolynomial, pk: RNSPublicKey) -> RNSCiphertext:
        """Encrypts a previously encoded plaintext into a ciphertext using the public key
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        # v = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        v = RNSPolynomial(self.B, self.C, self.roots, coeffs=[1] * (self.m // 2))  # This polynomial with 1s as coefficients is a placeholder until I understand the purpose of the term "v"
        e0 = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_gaussian_coeffs(self.m // 2))
        e1 = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_gaussian_coeffs(self.m // 2))

        # c = ((v, v) * self.pk) + (m + e0, e1)
        c0 = (v * pk[0]) + m + e0
        c1 = (v * pk[1]) + e1

        c = RNSCiphertext(c0, c1, self.B, self.C, self.p, self.q0, self.q, self.roots) #c0, c1, self.P, self.q0, self.delta, self.L)
        return c

    def decrypt(self, c: RNSCiphertext, sk: RNSPrivateKey) -> RNSPolynomial:
        """Decrypts a ciphertext using the secret key and returns a plaintext polynomial
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        return c.b + (c.a * sk[1])
