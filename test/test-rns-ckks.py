from playFHE.ckks import rns_ckks_scheme

print("=== RNS CKKS TEST ===\n")


N = 2**3
q = 35
q0 = 40
L = 2
k = 1
#P = 2**90
p = 38


print("Setup")
rnsckks = rns_ckks_scheme.RNSCKKS(N, p, k, q0, q, L, 3)
print("Keygen")
pk, sk = rnsckks.keygen()


m1 = [0.1]*(N//2)
m2 = [0.03]*(N//2)

print("Encoding m1")
p1 = rnsckks.encode(m1)
print("Encoding m2")
p2 = rnsckks.encode(m2)


print("\n")

print("p1", p1)
print("p2", p2)

p3 = p1 + p1
print("p3", p3)

print("m'1:", rnsckks.decode(p1))
print("m'2:", rnsckks.decode(p2))
print("\n")
print("p3=p1+p1:", rnsckks.decode(p3))

p4 = p3 * p2
p4.rescale()
print("p4=p3*p2:", rnsckks.decode(p4))

p5 = p3 * p1
p5.rescale()
print("p5=p3*p1:", rnsckks.decode(p5))

c1 = rnsckks.encrypt(p1, pk)
c2 = rnsckks.encrypt(p2, pk)
c3 = c1 + c2
c4 = c1 - c2

print("\nc1:", c1)
print("m1'", rnsckks.decode(rnsckks.decrypt(c1, sk)))
print("m2'", rnsckks.decode(rnsckks.decrypt(c2, sk)))

print("m3'", rnsckks.decode(rnsckks.decrypt(c3, sk)))
print("m4'", rnsckks.decode(rnsckks.decrypt(c4, sk)))
