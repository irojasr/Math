# %%
from admcycles import *
from math import factorial

# %%
def pslambdaclass(d, g, n):
    L = lambdaclass(d, g, n)
    i = 1
    while i < d + 1:
        # Construct V: [g - i, 1, 1, ..., 1]
        V = [g - i] + [1] * i

        # Construct H:
        # [[1, 2, ..., n + i], [n + i + 1], [n + i + 2], ..., [n + i + i]]
        H = [list(range(1, n + i + 1))] + [[n + i + k] for k in range(1, i + 1)]

        # Construct E: [(1, i+1), (2, i+2), ..., (i, i+i)]
        E = [(n+k,n+i+k) for k in range(1, i + 1)]

        # Construct A: [lambdaclass(i, g - i, i), fundclass(1, 1), ..., fundclass(1, 1)]
        A = [lambdaclass(d-i, g - i, n + i)] + [fundclass(1, 1)] * i

        # Accumulate into L
        L += (1/factorial(i))*StableGraph(V, H, E).boundary_pushforward(A) 
        ## Curiously factorial can be removed and Matt's identity still holds

        i += 1

    return L


# %%
def pslambdaclass_iterative(d, g, n):
    """
    Iterative pslambdaclass that works for general n.
    Reuses lists between iterations but initializes them using n.
    """
    L = lambdaclass(d, g, n)

    # i = 1 initial setup (use general formula with n)
    i = 1
    V = [g - i] + [1] * i
    H = [list(range(1, n + i + 1))] + [[n + i + k] for k in range(1, i + 1)]
    E = [(n + k, n + i + k) for k in range(1, i + 1)]
    A = [lambdaclass(d - i, g - i, n + i)] + [fundclass(1, 1)] * i

    L += (1 / factorial(i)) * StableGraph(V, H, E).boundary_pushforward(A)

    # subsequent iterations: update incrementally but maintain correct labeling using n
    for i in range(2, d + 1):
        # Update V
        V[0] = g - i
        V.append(1)

        # Append the new label to the main block H[0]
        H[0].append(n + i)

        # Rebuild singleton leg blocks to reflect the new i:
        # H[1:] = [[n+i+1], [n+i+2], ..., [n+2*i]]
        H[1:] = [[n + i + k] for k in range(1, i + 1)]

        # Recompute E for current i (cheap)
        E = [(n + k, n + i + k) for k in range(1, i + 1)]

        # Update A (first entry changes; then append fundclasses)
        A[0] = lambdaclass(d - i, g - i, n + i)
        A = [A[0]] + [fundclass(1, 1)] * i

        L += (1 / factorial(i)) * StableGraph(V, H, E).boundary_pushforward(A)

    return L


# %%
from itertools import product

def OtherSide(d,g,n):
    """
    Build the 'other side' sum described:
      result starts as fundclass(g,n)
      for i = 1..d:
        build V,H,E
        let G = StableGraph(V,H,E)
        expand prod_{k=1..i} (psi_{n+k} - psi_{n+k+i})
        for each bit-vector b in {0,1}^i:
            coef = (-1)^{sum(b)}
            first_entry = product_{k: b_k==1} psiclass(n+k, g-i, n+i)
            tail_entries = [ psiclass(1,1,1) if b_k==0 else fundclass(1,1) for k=1..i ]
            A = [ first_entry ] + tail_entries
            result += coef * norm * G.boundary_pushforward(A)
    Returns result.
    """

    result = fundclass(g, n)

    for i in range(1, d + 1):
        # Construct V: [g - i, 1, 1, ..., 1]
        V = [g - i] + [1] * i

        # Construct H: [[1,2,..., n+i], [n+i+1], [n+i+2], ..., [n+i+i]]
        H = [list(range(1, n + i + 1))] + [[n + i + k] for k in range(1, i + 1)]

        # Construct E: [(n+1,n+i+1), (n+2,n+i+2), ..., (n+i, n+2i)]
        E = [(n + k, n + i + k) for k in range(1, i + 1)]

        # Build StableGraph once
        G = StableGraph(V, H, E)

        # Normalization factor (preserve previous behavior)
        norm = 1 / factorial(i)

        # Expand the product into 2^i monomials via bit-vectors
        terms = list(product((0, 1), repeat=i))  # each term is a tuple of 0/1 length i

        # For each bit-vector build the corresponding A and accumulate
        for idx, bits in enumerate(terms, start=1):
            # coefficient sign = (-1)^{# of ones}
            coef = (-1) ** (i+sum(bits))

            # Build first entry: product of psiclasses for every bit==1
            # If no bits==1, first_entry is the multiplicative identity '1' by default.
            first_entry = fundclass(g-i,n+i)
            any_psi_in_first = False
            for k, b in enumerate(bits, start=1):
                if b == 1:
                    any_psi_in_first = True
                    # multiply into first_entry
                    # assumes psiclass(...) returns an object supporting multiplication
                    first_entry = first_entry * psiclass(n + k, g - i, n + i)

            # If your ring does not accept integer 1 as multiplicative identity,
            # replace the '1' above with the appropriate identity object, e.g. unitclass()

            # Build tail entries: psiclass(1,1,1) when bits[k]==0, else fundclass(1,1)
            tail = []
            for k, b in enumerate(bits, start=1):
                if b == 0:
                    tail.append(psiclass(1, 1, 1))
                else:
                    tail.append(fundclass(1, 1))

            # Full A vector: first entry then the i leg-entries
            A = [first_entry] + tail

            # Accumulate into result
            result += coef * norm * G.boundary_pushforward(A)

    return result


# %%
def total_lambda(g, n):
    """Compute (1 + λ₁ + λ₂ + ... + λ_g) using pslambdaclass."""
    result = fundclass(g, n)  # this represents '1'
    for d in range(1, g + 1):
        result += pslambdaclass_iterative(d, g, n)
    return result


def total_dual_lambda(g, n):
    """Compute (1 - λ₁ + λ₂ - λ₃ + ... + (-1)^g λ_g) using pslambdaclass."""
    result = fundclass(g, n)  # start with '1'
    for d in range(1, g + 1):
        sign = (-1) ** d
        result += sign * pslambdaclass_iterative(d, g, n)
    return result

# %%
g=2
n=1
##Check Matthew relation on M_(g,n)ps
TotPS=total_lambda(g,n)*total_dual_lambda(g,n)
MattIsZero = (TotPS-OtherSide(g,g,n)).is_zero()
print(MattIsZero)

# %%



