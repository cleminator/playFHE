from __future__ import annotations

from playFHE import util
from playFHE.util import Number

from typing import Union, TYPE_CHECKING
import copy

if TYPE_CHECKING:
    from playFHE.ckks.ciphertext import Ciphertext

class Polynomial:
    """
    This class represents a polynomial of the integer polynomial ring Zq[X]/(X^N + 1)
    The coefficients are integers modulo q, and are represented as a list with the highest order coefficient at index 0
    """
    coeffs: list[int | float]
    q: int | None
        
    def __init__(self, coeffs: list[int | float], q: None | int = None):
        self.q = q #Optional modulus (encoded plaintexts do not have moduli, but ciphertexts do
        self.n = len(coeffs) #Number of coefficients

        self.coeffs = coeffs
        self.mod()
    
    def __str__(self) -> str:
        """Return a string representation of the polynomial."""
        strng = "".join(str(self.coeffs))
        if self.q:
            strng += " (q = " + str(self.q) + ")"
        return strng
    
    
    ##################################################
    
    
    def mod(self):
        if self.q is not None:
            self.coeffs = [util.mod(c, self.q) for c in self.coeffs]
    
    def rescale(self, ql: int):
        """Operation to get rid of the extra scaling factor after multiplying two encoded polynomials
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3"""
        self.coeffs = [(c // ql) for c in self.coeffs]
        self.q //= ql
        self.mod()
    
    def mod_reduction(self, ql: int):
        """Operation to scale down the modulus without scaling down the coefficients; used to even moduli of two ciphertexts on different levels before multiplication
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3 "Homomorphic Operations of Ciphertexts at different levels"
        """
        self.q //= ql
        self.mod()
    

    ##################################################

    def solve(self, x: Number) -> Number:
        """Simple function to solve the polynomial for x (used in decoding procedure); no optimizations performed"""
        result = 0
        for i, coeff in enumerate(self.coeffs):
            result += coeff * (x ** i)
        return result

    ##################################################

    def add_poly(self, other: Polynomial) -> Polynomial:
        #print(self.coeffs)
        #print(other.coeffs)
        #print(self.q)
        #print(other.q)
        max_len = max(self.n, other.n)
        result = [0] * max_len

        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 + c2
        res = Polynomial(result, self.q)
        res.mod()
        return res

    def add_scalar(self, scalar: int | float) -> Polynomial:
        res = copy.deepcopy(self)
        res.coeffs = [c + scalar for c in res.coeffs]
        res.mod()
        return res

    def __add__(self, other: Polynomial | int | float) -> Polynomial:
        """Adding coefficients of two polynomials"""
        if isinstance(other, Polynomial):
            return self.add_poly(other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.add_scalar(other)
        else:
            return NotImplemented

    def sub_poly(self, other: Polynomial) -> Polynomial:
        max_len = max(self.n, other.n)
        result = [0] * max_len

        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 - c2
        res = Polynomial(result, self.q)
        res.mod()
        return res

    def sub_scalar(self, scalar: int | float) -> Polynomial:
        res = copy.deepcopy(self)
        res.coeffs = [c - scalar for c in res.coeffs]
        res.mod()
        return res

    def __sub__(self, other: Polynomial | int | float) -> Polynomial:
        """Subtracting coefficients of two polynomials"""
        if isinstance(other, Polynomial):
            return self.sub_poly(other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.sub_scalar(other)
        else:
            return NotImplemented
    
    def mult_scalar(self, scalar: float | int) -> Polynomial:
        """Multiplication of each coefficient with a constant"""
        res = copy.deepcopy(self)
        res.coeffs = [scalar * c for c in res.coeffs]
        res.mod()
        return res

    def mult_poly(self, poly: Polynomial) -> Polynomial:
        return self.negacyclic_convolution(self, poly)
    
    def __mul__(self, other: Polynomial | int | float) -> Polynomial | int | float:
        """Multiplication of two polynomials modulo (X^N +1); not optimized"""
        if isinstance(other, Polynomial):
            return self.mult_poly(other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.mult_scalar(other)
        else:
            return NotImplemented


    def negacyclic_convolution(self, poly1: Polynomial, poly2: Polynomial) -> Polynomial:
        """
        Perform negacyclic convolution of two polynomials; Not yet FTT-based
        """
        n = len(poly1.coeffs)
        m = len(poly2.coeffs)

        if n != m:
            raise ValueError("Negacyclic convolution requires polynomials of the same length.")

        if poly1.q != poly2.q:
            raise ValueError("Polynomial multiplication requires the same modulus q for both polynomials")

        # Initialize the result coefficients
        result = [0] * n

        # Perform direct computation of negacyclic convolution
        for i in range(n):
            for j in range(n):
                index = (i + j) % n
                value = poly1.coeffs[i] * poly2.coeffs[j]

                if i + j >= n:
                    result[index] -= value  # Negacyclic wraparound
                else:
                    result[index] += value

                if poly1.q is not None or poly2.q is not None:
                    result[index] = util.mod(result[index], poly1.q)

        res = Polynomial(result, poly1.q)
        res.mod()

        return res
    
    

