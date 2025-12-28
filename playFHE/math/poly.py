from __future__ import annotations

from playFHE import util
from playFHE.util import Number

from typing import Union, TYPE_CHECKING

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
        if q is not None:
            self.coeffs = [self.mod(c) for c in coeffs]
        else:
            self.coeffs = coeffs

    
    def __str__(self) -> str:
        """Return a string representation of the polynomial."""
        strng = "".join(str(self.coeffs))
        if self.q:
            strng += " (q = " + str(self.q) + ")"
        return strng
    
    
    ##################################################
    
    
    def mod(self, val: int) -> int:
        return util.mod(val, self.q)
    
    def rescale(self, ql: int):
        """Operation to get rid of the extra scaling factor after multiplying two encoded polynomials
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3"""
        self.coeffs = [(c // ql) for c in self.coeffs]
        self.q //= ql
    
    def mod_reduction(self, ql :int):
        """Operation to scale down the modulus without scaling down the coefficients; used to even moduli of two ciphertexts on different levels before multiplication
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3 "Homomorphic Operations of Ciphertexts at different levels"
        """
        self.q //= ql
    
    

    
    ##################################################

    def solve(self, x: Number) -> Number:
        """Simple function to solve the polynomial for x (used in decoding procedure); no optimizations performed"""
        result = 0
        for i, coeff in enumerate(self.coeffs):
            result += coeff * (x ** i)
        return result

    def add(self, other: Polynomial) -> Polynomial:
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
            if self.q:
                result[i] = self.mod(result[i])  # % self.q
        return Polynomial(result, self.q)



    def __add__(self, other: Union[Polynomial, "Ciphertext"]) -> Union[Polynomial, "Ciphertext"]:
        """Adding coefficients of two polynomials"""
        from playFHE.ckks.ciphertext import Ciphertext
        if isinstance(other, Ciphertext):
            return other + self
        
        return self.add(other)

    def sub(self, other: Polynomial) -> Polynomial:
        max_len = max(self.n, other.n)
        result = [0] * max_len

        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 - c2
            if self.q:
                result[i] = self.mod(result[i])  # % self.q
        return Polynomial(result, self.q)

    def __sub__(self, other: Polynomial) -> Polynomial:
        """Subtracting coefficients of two polynomials"""

        return self.sub(other)
    
    def scalar_mult(self, scalar: float | int) -> Polynomial:
        """Multiplication of each coefficient with a constant"""
        cfs = [scalar * c for c in self.coeffs]
        return Polynomial(cfs, self.q)
    
    def __mul__(self, other: Union[Polynomial, int, float, "Ciphertext"]) -> Union[Polynomial, "Ciphertext"]:
        """Multiplication of two polynomials modulo (X^N +1); not optimized"""
        from playFHE.ckks.ciphertext import Ciphertext
        if isinstance(other, Polynomial):
            return self.negacyclic_convolution(self, other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.scalar_mult(other)
        elif isinstance(other, Ciphertext):
            return other * self
        else:
            return NotImplemented
            #raise Exception("Only Poly-Poly or Poly-int multiplication is possible")


    def negacyclic_convolution(self, poly1: Polynomial, poly2: Polynomial) -> Polynomial:
        """
        Perform negacyclic convolution of two polynomials; Not yet FTT-based
        """
        n = len(poly1.coeffs)
        m = len(poly2.coeffs)

        if n != m:
            raise ValueError("Negacyclic convolution requires polynomials of the same length.")

        # Initialize the result coefficients
        result = [0] * n

        # Perform direct computation of negacyclic convolution
        for i in range(n):
            for j in range(n):
                index = (i + j) % n
                value = util.mod(poly1.coeffs[i] * poly2.coeffs[j], self.q)
                if i + j >= n:
                    result[index] -= value  # Negacyclic wraparound
                else:
                    result[index] += value

        return Polynomial(result, self.q)
    
    

