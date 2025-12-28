from __future__ import annotations

from playFHE.math.poly import Polynomial
import copy

from typing import Union, TYPE_CHECKING

if TYPE_CHECKING:
    from playFHE.ckks.ckks_scheme import MultKey

class Ciphertext:

    b: Polynomial
    a: Polynomial
    P: int
    q0: int
    delta: int
    l: int

    def __init__(self, b: Polynomial, a: Polynomial, P: int, q0: int, delta: int, L: int):
        self.b = b
        self.a = a
        
        self.P = P
        
        self.q0 = q0
        self.delta = delta
        self.l = L

    def __str__(self) -> str:
        return "Ciphertext (q0: " + str(self.q0) + ", delta: " + str(self.delta) + ", l: " + str(self.l) + ")"
            
    ############################################################
    
    def rescale(self):
        """Performs the rescaling for each of the involved polynomials"""
        self.b.rescale(self.delta)
        self.a.rescale(self.delta)
        self.l -= 1

    def mod_reduction(self):
        """Performs modular reduction for each involved polynomial"""
        self.b.mod_reduction(self.delta)
        self.a.mod_reduction(self.delta)
        self.l -= 1
    
    ############################################################
    
    def add_constant(self, other: Polynomial | int | float) -> Ciphertext:
        """Addition of ciphertext + plaintext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 4.1)"""
        # self: ciphertext; other: constant
        b = self.b + other
        a = self.a
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
        
    def add_ciph(self, other: Ciphertext) -> Ciphertext:
        """Addition of ciphertext + ciphertext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        print("Adding b")
        b = self.b + other.b
        #print(self.b.coeffs)
        #print(other.b.coeffs)
        #print(b.coeffs)

        print("Adding a")
        a = self.a + self.a
        print("done")
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
    
    def __add__(self, other: Ciphertext | Polynomial | int | float) -> Ciphertext:
        """Overloaded operator which chooses the correct addition"""
        if isinstance(other, Ciphertext):
            return self.add_ciph(other)
        elif isinstance(other, Polynomial) or isinstance(other, int) or isinstance(other, float):
            return self.add_constant(other)
        else:
            return NotImplemented

    ################
    
    def sub_constant(self, other: Polynomial | int | float) -> Ciphertext:
        """Subtraction of ciphertext - plaintext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 4.1)"""
        b = self.b - other
        a = self.a
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
        
    def sub_ciph(self, other: Ciphertext):
        """Subtraction of ciphertext - ciphertext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        b = self.b - other.b
        a = self.a - other.a
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
    
    def __sub__(self, other: Polynomial | Ciphertext | int | float) -> Ciphertext:
        """Overloaded operator which chooses the correct subtraction"""
        if isinstance(other, Ciphertext):
            return self.sub_ciph(other)
        elif isinstance(other, Polynomial) or isinstance(other, int) or isinstance(other, float):
            return self.sub_constant(other)
        else:
            return NotImplemented

    ################

    def mult_constant(self, other: Polynomial) -> Ciphertext:
        """Multiplication of ciphertext * plaintext
        Source: https://eprint.iacr.org/2016/421.pdf (Section 4.1)"""
        b = self.b * other
        a = self.a * other
        cmult = Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
        cmult.rescale()
        return cmult
    
    def mult_ciph(self, other: tuple[Ciphertext, "MultKey"]) -> Ciphertext:
        """Multiplication of ciphertext * ciphertext with an evaluation key; Syntax is ciph3 = ciph1 * (ciph2, evk)
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""

        c1 = copy.deepcopy(self)
        c2 = copy.deepcopy(other[0])
        evk = copy.deepcopy(other[1])
        
        # Check if levels are different between both ciphertexts
        while True:
            if c1.l > c2.l:
                c1.mod_reduction()
            elif c1.l < other[0].l:
                c2.mod_reduction()
            else:
                break

        d0 = c1.b * c2.b
        d1 = c1.a * c2.b + c2.a * c1.b
        d2 = c1.a * c2.a

        # cmult = (d0, d1) + round((d2 * evk) / P) mod ql

        b_with_evk = evk[0] * d2
        a_with_evk = evk[1] * d2
        b_with_evk.rescale(self.P) #.coeffs = [round(c / self.P) for c in b_with_evk.coeffs]
        a_with_evk.rescale(self.P) #.coeffs = [round(c / self.P) for c in a_with_evk.coeffs]

        cmult0 = d0 + b_with_evk
        cmult1 = d1 + a_with_evk
        
        cmult = Ciphertext(cmult0, cmult1, c1.P, c1.q0, c1.delta, c1.l)
        cmult.rescale()
        return cmult
    
    def __mul__(self, other: Union[tuple[Ciphertext, "MultKey"], Polynomial]) -> Ciphertext:
        """Overloaded operator which chooses the correct multiplication"""
        # First option: Ciphertext-Ciphertext Multiplication
        #     Other needs to be a tuple of (Ciphertext, evk)
        #     Example: cmult = c1 * (c2, evk)
        # Second option: Ciphertext-Constant Mult
        #     Other needs to be an int or float
        if isinstance(other, tuple):
            if len(other) == 2 and isinstance(other[0], Ciphertext):
                return self.mult_ciph(other)
            else:
                return NotImplemented
        elif isinstance(other, Polynomial):
            return self.mult_constant(other)
        else:
            return NotImplemented

