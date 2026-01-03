from __future__ import annotations

import copy

from playFHE.math import ntt
from playFHE import util
from playFHE.util import Number
from playFHE.math.poly import Polynomial

from math import prod

from typing import Union, TYPE_CHECKING

if TYPE_CHECKING:
    from playFHE.ckks.rns_ciphertext import RNSCiphertext


class RNSLimb:

    coeffs: list[int]
    q: int
    root: int

    def mod(self):
        self.coeffs = [util.mod(c, self.q) for c in self.coeffs]

    def __init__(self, coeffs: list[int], q: int, root: int):
        self.coeffs = coeffs
        self.q = q
        self.root = root
        self.mod()
    def __copy__(self):
        return RNSLimb(self.coeffs[:], self.q, self.root)
    def __getitem__(self, item: int):
        return self.coeffs[item]
    def __setitem__(self, key: int, value: int):
        self.coeffs[key] = value
        self.mod()
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
    # D: RNS base of P*Q -> B+C (len(D) = k + L+1)
    # limbs: list of RNS limbs, each limb i is a list of coeffs mod q_i

    B: list[int]
    C: list[int]
    roots: dict[int, int]
    limbs: list[RNSLimb]
    ntt_domain: bool

    def __init__(self, B: list[int], C: list[int], roots: dict[int, int], is_ntt: bool = False, coeffs: list[int] | None = None):
        self.B = B[:]
        self.C = C[:]
        self.D = self.B + self.C
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
        self.limbs = [copy.copy(l) for l in limbs[:]] #limbs[:]]
        if not self.ntt_domain:
            self.convert_RNS_to_NTT()
        return self

    ###########################################################################

    ## These functions will transform each of the limbs to and from NTT domain

    def convert_RNS_to_NTT(self):
        if not self.ntt_domain:
            for l in self.limbs:
                #print("Converting to NTT limb")
                l.coeffs = ntt.fast_ntt(l.coeffs, l.q, l.root)
            self.ntt_domain = True

    def convert_NTT_to_RNS(self):
        if self.n:
            for l in self.limbs:
                l.coeffs = ntt.fast_intt(l.coeffs, l.q, l.root)
            self.ntt_domain = False

    ###########################################################################

    def conv(self, limbs: list[RNSLimb], base1: list[int], base2: list[int]) -> list[RNSLimb]:
        # Convert limbs from basis1 (e.g. B) to basis2 (e.g. C)
        res: list[RNSLimb] = []

        q = base1[:] # C
        p = base2[:] # B
        l = len(q)
        k = len(p)

        Q = prod(q)

        qhat = []
        qhat_inv = []

        for j in range(l):
            qhat.append(Q // q[j])  # Product of all q except qj
            qhat_inv.append(util.findMultInv(qhat[j], q[j]))

        for i in range(k):
            limb_coeffs = []

            for c in range(self.n):
                a = 0
                for j in range(l):
                    a_term = util.mod(limbs[j][c] * qhat_inv[j], q[j])
                    a_term = a_term * qhat[j]
                    a += a_term
                limb_coeffs.append(util.mod(a, p[i]))

            res.append(RNSLimb(limb_coeffs, p[i], self.roots[p[i]]))

        return res

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

    def mod(self):
        for l in self.limbs:
            l.mod()

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
        self.mod()


    def mod_reduction(self):
        self.limbs.pop()
        self.C.pop()
        self.mod()

    ######################################

    def add_poly(self, other: RNSPolynomial) -> RNSPolynomial:
        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, lb, q in zip(self.limbs, other.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] + lb[i], q)
            limbs.append(RNSLimb(lc, la.q, la.root))  # lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def add_scalar(self, scalar: int | float) -> RNSPolynomial:
        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, q in zip(self.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] + scalar, q)
            limbs.append(RNSLimb(lc, la.q, la.root))  # lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def __add__(self, other: RNSPolynomial | int | float) -> RNSPolynomial:

        if self.n != other.n:
            raise Exception("Polynomials need to have equal ring dimension")

        if isinstance(other, RNSPolynomial):
            return self.add_poly(other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.add_scalar(other)
        else:
            return NotImplemented

    ###

    def sub_poly(self, other: RNSPolynomial) -> RNSPolynomial:
        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, lb, q in zip(self.limbs, other.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] - lb[i], q)
            limbs.append(RNSLimb(lc, la.q, la.root))  # lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def sub_scalar(self, scalar: int | float) -> RNSPolynomial:
        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, q in zip(self.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] - scalar, q)
            limbs.append(RNSLimb(lc, la.q, la.root))  # lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def __sub__(self, other: RNSPolynomial | int | float) -> RNSPolynomial:

        if self.n != other.n:
            raise Exception("Polynomials need to have equal ring dimension")

        if isinstance(other, RNSPolynomial):
            return self.sub_poly(other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.sub_scalar(other)
        else:
            return NotImplemented

    ###

    def mult_scalar(self, scalar: int | float) -> RNSPolynomial:
        limbs = []
        for l, q in zip(self.limbs, self.C):
            limbs.append(RNSLimb([util.mod(scalar * c, q) for c in l], l.q, l.root))
        return RNSPolynomial(self.B, self.C, self.roots, self.ntt_domain).set_limbs(limbs)

    def mult_poly(self, other: RNSPolynomial) -> RNSPolynomial:
        p1 = copy.deepcopy(self)
        p2 = copy.deepcopy(other)
        if not p1.ntt_domain:
            p1.convert_RNS_to_NTT()
        if not p2.ntt_domain:
            p2.convert_RNS_to_NTT()

        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, lb, q in zip(self.limbs, other.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] * lb[i], q)
            limbs.append(RNSLimb(lc, la.q, la.root))  # lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)


    def __mul__(self, other: Union[RNSPolynomial, int, float, "RNSCiphertext"]) -> Union[RNSPolynomial, "RNSCiphertext"]:
        """Multiplication of two polynomials modulo (X^N +1); not optimized"""

        if isinstance(other, RNSPolynomial):
            return self.mult_poly(other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.mult_scalar(other)
        else:
            return NotImplemented
