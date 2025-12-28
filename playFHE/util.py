import random
from math import floor, ceil

Number = float | complex

def transpose(matrix: list[list[Number]]) -> list[list[Number]]:
    rows = len(matrix)
    cols = len(matrix[0])

    transposed = [[] for _ in range(cols)]

    for row in matrix:
        for j, val in enumerate(row):
            transposed[j].append(val)

    return transposed

def vdot(a: list[Number], b: list[Number]) -> Number:
    if len(a) != len(b):
        raise ValueError("Inputs must have the same length.")
    return sum((x.conjugate() if isinstance(x, complex) else x) * y for x, y in zip(a, b))

def matmul(a: list[list[Number]], b: list[list[Number]] | list[Number]) -> list[list[Number]]:
    if len(a[0]) != len(b):
        raise ValueError("Incompatible matrix dimensions for multiplication")
    if not isinstance(b[0], list):
        b = [[bb] for bb in b]
    result = [[0] * len(b[0]) for _ in range(len(a))]
    for i in range(len(a)):  # Iterate over rows of A
        for j in range(len(b[0])):  # Iterate over columns of B
            for k in range(len(b)):  # Iterate over rows of B (columns of A)
                result[i][j] += a[i][k] * b[k][j]

    return result


def coordinate_wise_random_rounding(coordinates: list[float]) -> list[int]:
    """Rounds coordinates randomly.
    Source: https://web.eecs.umich.edu/~cpeikert/pubs/toolkit.pdf"""

    rounded_coordinates = []

    for c in coordinates:
        r = c - floor(c)  # Get integral rest (digits after decimal point)
        rand = random.random()  # Random float between 0 and 1

        if rand <= r:
            rounded_coordinates.append(ceil(c))
        else:
            rounded_coordinates.append(floor(c))

    return rounded_coordinates

def gaussian_elimination(matrix: list[list[Number]], vector: list[Number]) -> list[complex]:
    """Simple solver for linear equation systems using gaussian elimination; not optimized"""
    n = len(matrix)

    # Ensure the vector supports complex numbers
    vector = [complex(v) for v in vector]

    # Forward elimination
    for i in range(n):
        # Pivot selection (partial pivoting)
        max_row = i + max(range(n - i), key=lambda k: abs(matrix[i + k][i]))
        if matrix[max_row][i] == 0:
            raise ValueError("Matrix is singular or nearly singular.")

        # Swap rows in the matrix and the vector
        matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
        vector[i], vector[max_row] = vector[max_row], vector[i]

        # Make the diagonal element 1 and eliminate below
        for j in range(i + 1, n):
            factor = matrix[j][i] / matrix[i][i]
            for k in range(i, n):
                matrix[j][k] -= factor * matrix[i][k]
            vector[j] -= factor * vector[i]

    # Back substitution
    solution = [0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = sum(matrix[i][j] * solution[j] for j in range(i + 1, n))
        solution[i] = (vector[i] - sum_ax) / matrix[i][i]

    return solution

def mod(val: int, mod: int) -> int:
    """Modulo operation with representation between -q/2 and q/2
    Source: https://eprint.iacr.org/2016/421.pdf Section 2.1 Basic"""
    z_mod_q = val % mod
    if z_mod_q > mod / 2:
        z_mod_q -= mod
    return z_mod_q

def findMultInv(a: int, n: int) -> int:
    """Determine the multiplicative inverse of a number a mod n"""
    # Compute Extended Euclidian algorithm
    # a*x = 1 mod n
    # 1 = a*x + n*y = gcd(a, n)
    gcd, x, _ = extGCD(a, n)
    return x % n

def extGCD(a: int, b: int) -> tuple[int, int, int]:
    """Extended euclidian algorithm, used mainly for determining mult inverse
    ax + by = gcd(a, b)
    """
    if b == 0:
        return a, 1, 0
    else:
        gcd, x1, y1 = extGCD(b, a % b)
        x = y1
        y = x1 - (a // b) * y1
        return gcd, x, y
        
def GCD(a: int, b: int) -> int:
    """Greatest common divisor"""
    gcd, _, _ = extGCD(a, b)
    return gcd


def totient(m: int) -> int:
    """Calculate Euler's Totient function φ(m) manually."""
    count = 0
    for i in range(1, m):
        if GCD(i, m) == 1:
            count += 1
    return count

def is_prime(n: int) -> bool:
    """Simple test whether a number is prime; not optimized"""
    if n <= 1:
        return False
    if n == 2:
        return True
    if n%2 == 0:
        return False
    
    for i in range(3, int(n**0.5)+1, 2):
        if n % i == 0:
            return False
    return True


def find_next_prime(n: int) -> int:
    """Determines the next prime number after the number n; not optimized"""
    p = n
    if is_prime(p):
        return p
    # Go to the next odd number after n
    if n % 2 == 0:
        p+=1
    else:
        p+=2
    
    #Now loop through every odd number to check for primality
    while True:
        if is_prime(p):
            return p
        p+=2

def mod_exp(base: int, exp: int, q: int) -> int:
    """ Simple modular exponentiation function"""
    if not q:
        raise Exception("No modulus defined!")
    result = 1
    while exp > 0:
        if exp % 2 == 1:  # If exp is odd
            result = mod(result * base, q)
        base = mod(base * base, q)
        exp //= 2
    return result

def find_primitive_nth_root_of_unity(n: int, q: int) -> int:
    while True:
        x = random.randint(1, q - 1)
        g = mod_exp(x, (q - 1) / n, q)
        if mod_exp(g, n / 2, q) != 1:
            return g

def find_2nth_root_of_unity(N: int, q: int) -> int:
    root2n =  find_primitive_nth_root_of_unity(2*N, q)
    return root2n


def sample_uniform_coeffs(n: int, q: int) -> list[int]:
    """Uniformly sample coefficients from Z_q"""
    return [random.randint(round(-(q-1)/2), round((q-1)/2)) for i in range(0, n)]

def sample_gaussian_coeffs(n: int, sigma: float = 3.2) -> list[int]:
    """Discrete gaussian distribution with variance sigma^2"""
    return [abs(round(random.gauss(0, sigma**2))) for i in range(0, n)]

def sample_uniform_ternary_coeffs(n: int) -> list[int]:
    """Uniform ternary distribution"""
    return [random.choice([-1, 0, 1]) for i in range(0, n)]

def sample_sparse_ternary_coeffs(n: int, p: float = 0.5) -> list[int]:
    coeffs = []
    for _ in range(n):
        r = random.random()
        if r < p / 2:
            coeffs.append(1)
        elif r < p:
            coeffs.append(-1)
        else:
            coeffs.append(0)
    return coeffs
