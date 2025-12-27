from playFHE.math.poly import Polynomial
import copy

from playFHE.ckks.ciphertext import Ciphertext


class RNSCiphertext(Ciphertext):
    def __init__(self, b, a, B, C, p, q0, q, roots): #P, q0, delta, L):
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


    def __str__(self):
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
