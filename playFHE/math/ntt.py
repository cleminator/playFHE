from playFHE import util


# Schoolbook version implemented after https://eprint.iacr.org/2024/585.pdf

##################

def ntt_psi(coeffs, n, q, root2n):#, q, root2n):
    a_hat = []

    for j in range(0, n):
        #print("j:", j)
        aj = 0
        for i in range(0, n):
            #print("i:", i)
            aj += util.mod_exp(root2n, 2 * i * j + i, q) * coeffs[i]#root2n**(2*i*j+i) * coeffs[i]
        a_hat.append(util.mod(aj, q))
    return a_hat#Polynomial(a_hat, q)

def intt_psi(coeffs, n, q, root2n):#, q, root2n):
    #print("START")
    a = []
    inv_n = util.findMultInv(n, q)
    inv_root2n = util.findMultInv(root2n, q)

    for i in range(0, n):
        #print("i:", i)
        ai = 0
        for j in range(0, n):
            #print("j:", j)
            tmp = util.mod_exp(inv_root2n, 2 * i * j + i, q)
            ai += util.mod((tmp * coeffs[j]), q)
        ai %= q
        ai *= inv_n
        a.append(ai % q)
    return a#Polynomial(a, a_hat.q)

#############################

def bit_reverse(x, n):
    return int(''.join(reversed(bin(x)[2:].zfill(n))), 2)

def generate_psi_rev(n, q, psi):
    num_bits = n.bit_length() - 1  # Number of bits needed to represent n-1
    psi_powers = [(pow(psi, i, q)) for i in range(n)]  # Powers of ψ modulo q
    psi_rev = [0] * n  # Initialize the Ψrev table

    # Populate Ψrev with bit-reversed indices
    for i in range(n):
        reversed_index = bit_reverse(i, num_bits)
        psi_rev[reversed_index] = psi_powers[i]

    return psi_rev

###########

# Fast NTT/INTT implemented after algorithms 1 and 2 in https://eprint.iacr.org/2016/504.pdf

def fast_ntt(a, q, psi):
    a = a[:]
    n = len(a)
    t = n
    psi_rev = generate_psi_rev(n, q, psi)
    m = 1
    while m < n:
        t = t//2
        for i in range(m):
            j1 = 2 * i * t
            j2 = j1 + t - 1
            S = psi_rev[m+i]
            for j in range(j1, j2+1):
                U = a[j]
                V = a[j+t] * S
                a[j] = util.mod(U + V, q)
                a[j+t] = util.mod(U - V, q)
        m *= 2
    return a


def fast_intt(a, q, psi):
    a = a[:]
    n = len(a)
    t = 1
    psi_inv_rev = generate_psi_rev(n, q, util.findMultInv(psi, q))
    inv_n = util.findMultInv(n, q)

    m = n
    while m > 1:
        j1 = 0
        h = m//2
        for i in range(h):
            j2 = j1 + t - 1
            S = psi_inv_rev[h + i]
            for j in range(j1, j2+1):
                U = a[j]
                V = a[j+t]
                a[j] = util.mod(U + V, q)
                a[j+t] = util.mod((U - V) * S, q)
            j1 = j1 + 2*t
        t = 2*t
        m //= 2

    for j in range(n):
        a[j] = util.mod(a[j] * inv_n, q)

    return a
