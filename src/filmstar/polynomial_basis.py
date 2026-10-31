"""Small polynomial helpers shared by planar and spatial patch solvers."""


def monomials(degree, coordinates):
    first, second = coordinates
    return tuple(first**(n-b)*second**b for n in range(degree+1)
                 for b in range(n+1))
