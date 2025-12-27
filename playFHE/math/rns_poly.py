from __future__ import annotations

from copy import copy

from playFHE.math import ntt
from playFHE import util
from playFHE.util import Number
from playFHE.math.poly import Polynomial

from typing import Union, TYPE_CHECKING

if TYPE_CHECKING:
    from playFHE.ckks.rns_ciphertext import RNSCiphertext


class RNSLimb:

    coeffs: list[int]
    q: int
    root: int

    def __init__(self, coeffs: list[int], q: int, root: int):
        self.coeffs = coeffs
        self.q = q
        self.root = root
    def __copy__(self):
        return RNSLimb(self.coeffs[:], self.q, self.root)
    def __getitem__(self, item: int):
        return self.coeffs[item]
    def __setitem__(self, key: int, value: int):
        self.coeffs[key] = value
    def __iter__(self):
        for c in self.coeffs:
            yield c
    def __len__(self):
        return len(self.coeffs)
    def __str__(self) -> str:
        strng = "".join(str(self.coeffs))
        strng += " (Modulus: " + str(self.q) + ", Root: " + str(self.root) + ")"
        return strng



class RNSPolynomial(Polynomial):
    # B: RNS base of P [p0, p1, ..., pk-1] (len(B) = k)
    # C: RNS base of Q [q0, q1, ..., qL] (len(C) = L+1)
    # limbs: list of RNS limbs, each limb i is a list of coeffs mod q_i

    B: list[int]
    C: list[int]
    roots: dict[int, int]
    limbs: list[RNSLimb]
    ntt_domain: bool

    def __init__(self, B: list[int], C: list[int], roots: dict[int, int], is_ntt: bool = False, coeffs: list[int] | None = None):
        self.B = B[:]
        self.C = C[:]
        self.roots = roots
        self.limbs = []

        self.n = 0
        self.ntt_domain = is_ntt
        if coeffs is not None:
            self.create_limbs(coeffs[:])


    def __str__(self) -> str:
        strng = "".join(str(self.get_coeffs()))
        strng += " (Modulus = " + str(self.get_modulus()) + ")"
        if self.ntt_domain:
            strng += " (NTT)"
        return strng


    def get_coeffs(self) -> list[int] | None:

        if len(self.C) == 1:
            return self.limbs[0].coeffs
        if len(self.C) == 0:
            raise Exception("No limbs exist")

        # x = l1*m2*n2 + l2*m1*n1

        coeffs = [0]*len(self.limbs[0])
        l1 = self.limbs[0]
        q1 = self.C[0]

        for i in range(1, len(self.C)):
            _, m1, m2 = util.extGCD(q1, self.C[i])
            for j in range(len(l1)):
                coeffs[j] = util.mod(l1[j] * m2 * self.C[i] + self.limbs[i][j] * m1 * q1, self.get_modulus())
            l1 = coeffs
            q1 = q1 * self.C[i]
        return coeffs

    def get_modulus(self) -> int:
        q = 1
        for l in self.limbs:#range(len(self.C)):
            q *= l.q#self.C[i]
        return q

    ###########################################################################

    def create_limbs(self, coeffs: list[int]):
        # Generate list of limbs from base (Q) and coefficients
        self.n = len(coeffs[:])
        for i in range(len(self.C)):
            #print("Generating limb", i)
            l = RNSLimb([c % self.C[i] for c in coeffs[:]], self.C[i], self.roots[self.C[i]])
            self.limbs.append(l) #[c % self.C[i] for c in coeffs[:]])

        if not self.ntt_domain:
            self.convert_RNS_to_NTT()
            self.ntt_domain = True

    def set_limbs(self, limbs: list[RNSLimb]) -> RNSPolynomial:
        self.n = len(limbs[0])
        self.limbs = [copy(l) for l in limbs[:]] #limbs[:]]
        if not self.ntt_domain:
            self.convert_RNS_to_NTT()
        return self

    ###########################################################################

    ## These functions will transform each of the limbs to and from NTT domain

    def convert_RNS_to_NTT(self):
        for l in self.limbs:
            #print("Converting to NTT limb")
            l.coeffs = ntt.fast_ntt(l.coeffs, l.q, l.root)
        self.ntt_domain = True

    def convert_NTT_to_RNS(self):
        for l in self.limbs:
            l.coeffs = ntt.fast_intt(l.coeffs, l.q, l.root)
        self.ntt_domain = False

    ###########################################################################

    def conv(self, base1, base2):
        # Convert limbs from basis1 (e.g. B) to basis2 (e.g. C)
        pass

    def mod_up(self):
        pass

    def mod_down(self):
        pass

    def solve(self, x: Number) -> Number:
        """Simple function to solve the polynomial for x (used in decoding procedure); no optimizations performed"""
        result = 0
        for i, coeff in enumerate(self.get_coeffs()):
            result += coeff * (x ** i)
        return result

    ######################################

    def rescale(self):
        for l in self.limbs:
            limb = []
            for i in range(self.n):
                lp = l[i] - self.limbs[-1][i]
                lp = lp * util.findMultInv(self.limbs[-1].q, l.q)
                lp = util.mod(lp, l.q)
                limb.append(lp)
            l.coeffs = limb

        self.limbs.pop()
        self.C.pop()


    def mod_reduction(self):
        self.limbs.pop()
        self.C.pop()

    ######################################

    def __add__(self, other: Union[RNSPolynomial, "RNSCiphertext"]) -> Union[RNSPolynomial, "RNSCiphertext"]:
        from playFHE.ckks.rns_ciphertext import RNSCiphertext
        if isinstance(other, RNSCiphertext):
            return other + self

        if self.n != other.n:
            raise Exception("Polynomials need to have equal ring dimension")


        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, lb, q in zip(self.limbs, other.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] + lb[i], q)
            limbs.append(RNSLimb(lc, la.q, la.root))#lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def __sub__(self, other: RNSPolynomial) -> RNSPolynomial:
        if self.n != other.n:
            raise Exception("Polynomials need to have equal ring dimension")

        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, lb, q in zip(self.limbs, other.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] - lb[i], q)
            limbs.append(RNSLimb(lc, la.q, la.root))#lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def scalar_mult(self, scalar: int | float) -> RNSPolynomial:
        limbs = []
        for l, q in zip(self.limbs, self.C):
            limbs.append(RNSLimb([util.mod(scalar * c, q) for c in l], l.q, l.root))
        return RNSPolynomial(self.B, self.C, self.roots, self.ntt_domain).set_limbs(limbs)

    def negacyclic_convolution(self, poly1: RNSPolynomial, poly2: RNSPolynomial) -> RNSPolynomial:
        limbs = []
        for la, lb in zip(poly1.limbs, poly2.limbs):
            n = len(la)
            m = len(lb)

            if n != m:
                raise ValueError("Negacyclic convolution requires polynomials of the same length.")

            result = [0] * n

            for i in range(self.n):
                result[i] = la[i] * lb[i]

            limbs.append(RNSLimb([util.mod(r, la.q) for r in result], la.q, la.root))
        return RNSPolynomial(poly1.B, poly1.C, poly1.roots, poly1.ntt_domain).set_limbs(limbs)

    def __mul__(self, other: Union[RNSPolynomial, int, float, "RNSCiphertext"]) -> Union[RNSPolynomial, "RNSCiphertext"]:
        """Multiplication of two polynomials modulo (X^N +1); not optimized"""
        from playFHE.ckks.rns_ciphertext import RNSCiphertext
        if isinstance(other, RNSPolynomial):
            return self.negacyclic_convolution(self, other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.scalar_mult(other)
        elif isinstance(other, RNSCiphertext):
            return other * self
        else:
            return NotImplemented
