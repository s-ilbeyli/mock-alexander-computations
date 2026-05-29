from itertools import permutations
import sympy as sp


def permanent(matrix):
    """
    Compute the permanent of a square matrix.

    The permanent is defined by

        per(A) = sum_{sigma in S_n} prod_i A[i, sigma(i)].

    This is the determinant formula without the sign of the permutation.
    """
    A = sp.Matrix(matrix)
    n, m = A.shape

    if n != m:
        raise ValueError("Permanent is defined only for square matrices.")

    total = 0

    for sigma in permutations(range(n)):
        term = 1
        for i in range(n):
            term *= A[i, sigma[i]]
        total += term

    return sp.expand(sp.simplify(total))


def check_polynomial(name, matrix, expected):
    """
    Compute the permanent of a matrix and compare it with an expected polynomial.
    """
    result = permanent(matrix)
    expected = sp.expand(sp.simplify(expected))
    difference = sp.expand(sp.simplify(result - expected))

    print(f"\n{name}")
    print("-" * len(name))
    print("Permanent:")
    print(result)
    print("Expected:")
    print(expected)

    if difference == 0:
        print("Status: OK")
    else:
        print("Status: NOT OK")
        print("Difference:")
        print(difference)

    return result


def check_multiplicativity(name, connected_matrix, first_matrix, second_matrix):
    """
    Check whether

        per(connected_matrix) = per(first_matrix) * per(second_matrix).

    This is the computational form of the multiplicativity formula.
    """
    connected = permanent(connected_matrix)
    first = permanent(first_matrix)
    second = permanent(second_matrix)
    expected = sp.expand(sp.simplify(first * second))
    difference = sp.expand(sp.simplify(connected - expected))

    print(f"\n{name}")
    print("-" * len(name))
    print("Permanent of connected sum:")
    print(connected)
    print("Product of permanents:")
    print(expected)

    if difference == 0:
        print("Status: OK")
    else:
        print("Status: NOT OK")
        print("Difference:")
        print(difference)

    return connected, expected


W = sp.symbols("W")


# ============================================================
# Matrices for the examples
#
# Rows: crossings
# Columns: unstarred regions
# Entries: local weights of the corresponding quadrant
# ============================================================

M_K1 = [
    [W**-1, 1, 0],
    [-W, 1, W**-1],
    [0, 1, -W],
]

M_K2 = [
    [W**-1 + 1, -W, 0, 0],
    [1, W**-1, -W, 0],
    [1, 0, W**-1, -W],
    [1, 0, 0, W**-1],
]

M_K3 = [
    [1, W**-1],
    [1 + W**-1, -W],
]

M_K1K2 = [
    [W**-1, 1, 0, 0, 0, 0, 0],
    [-W, 1, W**-1, 0, 0, 0, 0],
    [0, 1, -W, 0, 0, 0, 0],
    [0, 0, 0, W**-1, 0, 0, 1],
    [0, 0, 0, -W, W**-1, 0, 1],
    [0, 0, 0, 0, -W, W**-1, 1],
    [0, 0, 0, 0, 0, -W, 1 + W**-1],
]

M_K3K2 = [
    [W**-1, 1, 0, 0, 0, 0],
    [-W, 1 + W**-1, 0, 0, 0, 0],
    [0, 1, W**-1, 0, 0, 1],
    [0, 1, -W, W**-1, 0, 1],
    [0, 1, 0, -W, W**-1, 1],
    [0, 1, 0, 0, -W, 1 + W**-1],
]


# ============================================================
# Individual mock Alexander polynomials
# ============================================================

P_K1 = permanent(M_K1)
P_K2 = permanent(M_K2)
P_K3 = permanent(M_K3)

print("\nIndividual mock Alexander polynomials")
print("-------------------------------------")
print("nabla(K_1*; W) =", P_K1)
print("nabla(K_2*; W) =", P_K2)
print("nabla(K_3*; W) =", P_K3)


# ============================================================
# Multiplicativity checks
# ============================================================

P_K1K2, expected_K1K2 = check_multiplicativity(
    "K_1* # K_2*",
    M_K1K2,
    M_K1,
    M_K2,
)

P_K3K2, expected_K3K2 = check_multiplicativity(
    "K_3* # K_2*",
    M_K3K2,
    M_K3,
    M_K2,
)


# ============================================================
# Optional: print LaTeX output for the paper
# ============================================================

print("\nLaTeX output")
print("------------")
print(r"\nabla(K_{1*};W) =", sp.latex(P_K1))
print(r"\nabla(K_{2*};W) =", sp.latex(P_K2))
print(r"\nabla(K_{3*};W) =", sp.latex(P_K3))
print(r"\nabla(K_{1*}\#K_{2*};W) =", sp.latex(P_K1K2))
print(r"\nabla(K_{3*}\#K_{2*};W) =", sp.latex(P_K3K2))