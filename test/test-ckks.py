from playFHE.ckks import ckks_scheme
from timeit import default_timer as timer

print("=== CKKS TEST ===\n")

N = 2**2
delta = 2**40
q0 = 2**45
L = 4
P = 2**(50+4*30)

start = timer()
ckks = ckks_scheme.CKKS(N, P, q0, delta, L, 3)
print("Setup duration:", timer() - start)

start = timer()
pk, sk = ckks.keygen()
evk = ckks.evkeygen(sk)
print("Keygen duration:", timer() - start)

#m1 = [0.1, 0.2, 0.3, 0.4]
#m2 = [0.03, 0.04, 1, 1]
m1 = [0.2]*(N//2)
m2 = [0.1]*(N//2)

print("m1: ", m1)
print("m2: ", m2)
print("\n")

start = timer()
p1 = ckks.encode(m1)
print("Encoding duration", timer() - start)
p2 = ckks.encode(m2)
print("p1: ", p1)
print("p2: ", p2)

start = timer()
c1 = ckks.encrypt(p1, pk)
print("Encryption duration", timer() - start)
c2 = ckks.encrypt(p2, pk)

print("c1: ", c1)
print("c2: ", c2)

c_add = c1 + c2
c_add_const = c1 + p2

c_sub = c1 - c2
c_sub_const = c1 - p2
c_mult_const = c1 * p2
c_mult = c1 * (c2, evk)

print(c_mult)

print("")
print("e_add", ckks.decode(ckks.decrypt(c_add, sk)))
print("e_add_const", ckks.decode(ckks.decrypt(c_add_const, sk)))
print("e_sub", ckks.decode(ckks.decrypt(c_sub, sk)))
print("e_sub_const", ckks.decode(ckks.decrypt(c_sub_const, sk)))
print("e_mult", ckks.decode(ckks.decrypt(c_mult, sk)))
print("e_mult_const", ckks.decode(ckks.decrypt(c_mult_const, sk)))

