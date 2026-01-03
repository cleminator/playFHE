from __future__ import annotations

from playFHE.math.poly import Polynomial
import copy

from playFHE.ckks.ciphertext import Ciphertext
from playFHE.math.rns_poly import RNSPolynomial


class RNSCiphertext:

    b: RNSPolynomial
    a: RNSPolynomial
    B: list[int]
    C: list[int]
    p: int
    q0: int
    q: int
    k: int
    l: int
    roots: dict[int, int]

    def __init__(self, b: RNSPolynomial, a: RNSPolynomial, B: list[int], C: list[int], p: int, q0: int, q: int, roots: dict[int, int]): #P, q0, delta, L):
        self.b = b
        self.a = a

        self.B = B[:]
        self.C = C[:]
        self.p = p
        self.q0 = q0
        self.q = q
        self.k = len(B)
        self.l = len(C) - 1

        self.roots = roots


    def __str__(self) -> str:
        return "Ciphertext (q0: " + str(self.q0) + ", q: " + str(self.q) + ", l: " + str(self.l) + ")"

    ############################################################

    def rescale(self):
        """Performs the rescaling for each of the involved polynomials"""
        self.b.rescale()
        self.a.rescale()
        self.l -= 1

    def mod_reduction(self):
        """Performs modular reduction for each involved polynomial"""
        self.b.mod_reduction()
        self.a.mod_reduction()
        self.l -= 1

    ############################################################

    def add_plain(self, other: RNSPolynomial) -> RNSCiphertext:
        """Addition of ciphertext + plaintext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 4.1)"""
        # self: ciphertext; other: constant
        res = copy.deepcopy(self)

        res.b = res.b + other
        return res

    def add_ciph(self, other: RNSCiphertext) -> RNSCiphertext:
        """Addition of ciphertext + ciphertext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""

        res = copy.deepcopy(self)

        res.b = res.b + other.b
        res.a = res.a + other.a
        return res

    def __add__(self, other: RNSCiphertext | RNSPolynomial) -> RNSCiphertext:
        """Overloaded operator which chooses the correct addition"""
        if isinstance(other, RNSCiphertext):
            return self.add_ciph(other)
        elif isinstance(other, RNSPolynomial):
            return self.add_plain(other)
        else:
            return NotImplemented

    ################

    def sub_plain(self, other: RNSPolynomial) -> RNSCiphertext:
        """Subtraction of ciphertext - plaintext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 4.1)"""
        # self: ciphertext; other: constant
        res = copy.deepcopy(self)

        res.b = res.b - other
        return res

    def sub_ciph(self, other: RNSCiphertext) -> RNSCiphertext:
        """Subtraction of ciphertext - ciphertext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""

        res = copy.deepcopy(self)

        res.b = res.b - other.b
        res.a = res.a - other.a
        return res

    def __sub__(self, other: RNSCiphertext | RNSPolynomial) -> RNSCiphertext:
        """Overloaded operator which chooses the correct subtraction"""
        if isinstance(other, RNSCiphertext):
            return self.sub_ciph(other)
        elif isinstance(other, RNSPolynomial):
            return self.sub_plain(other)
        else:
            return NotImplemented

    ################