# fermat_determinants.py
#
# This script computes the determinant data attached to the two normalized
# Fermat sixers.
#
# Mathematical input:
#   the normalized matrices for the two orbit representatives
#
#       O_54 and O_18
#
#   over K = Q(omega), where omega^2 + omega + 1 = 0.
#
# Computational task:
#   for each orbit, compute the determinants of the ten pairwise differences
#   among
#
#       0, I_2, M_2, M_3, M_4.
#
# Output:
#   - the determinant list for the orbit of size 54;
#   - the determinant set, up to repetition, for the orbit of size 54;
#   - the determinant list for the orbit of size 18;
#   - the determinant set, up to repetition, for the orbit of size 18.
#
# These outputs are used to compute the image of the determinant square-class
# character
#
#       delta_K : PGL_2(K) -> K^*/(K^*)^2.


from dataclasses import dataclass


@dataclass(frozen=True)
class Eisenstein:
    """
    Element a + b*omega of Z[omega], where omega^2 + omega + 1 = 0.

    We store such an element as the pair (a,b).
    """
    a: int
    b: int = 0

    def __add__(self, other):
        other = to_eisenstein(other)
        return Eisenstein(self.a + other.a, self.b + other.b)

    def __neg__(self):
        return Eisenstein(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-to_eisenstein(other))

    def __mul__(self, other):
        other = to_eisenstein(other)

        # Since omega^2 = -1 - omega, we have
        #
        #   (a+b omega)(c+d omega)
        #   = ac + (ad+bc)omega + bd omega^2
        #   = (ac-bd) + (ad+bc-bd)omega.
        #
        return Eisenstein(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b
        )

    def __repr__(self):
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            if self.b == 1:
                return "omega"
            if self.b == -1:
                return "-omega"
            return f"{self.b}*omega"
        sign = "+" if self.b > 0 else "-"
        coeff = abs(self.b)
        if coeff == 1:
            return f"{self.a}{sign}omega"
        return f"{self.a}{sign}{coeff}*omega"


def to_eisenstein(x):
    """
    Convert an integer to an Eisenstein integer.
    """
    if isinstance(x, Eisenstein):
        return x
    return Eisenstein(x, 0)


zero = Eisenstein(0)
one = Eisenstein(1)
omega = Eisenstein(0, 1)


def mat_sub(A, B):
    """
    Difference of two 2 x 2 matrices over Z[omega].
    """
    return [
        [A[i][j] - B[i][j] for j in range(2)]
        for i in range(2)
    ]


def det(A):
    """
    Determinant of a 2 x 2 matrix over Z[omega].
    """
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def determinant_data(matrices):
    """
    Given the list

        [0, I_2, M_2, M_3, M_4],

    compute the determinants det(A-B) for all ten unordered pairs.
    """
    values = []

    for i in range(len(matrices)):
        for j in range(i + 1, len(matrices)):
            values.append(det(mat_sub(matrices[i], matrices[j])))

    return values


# The zero matrix and the identity matrix.

Z = [
    [zero, zero],
    [zero, zero]
]

Id = [
    [one, zero],
    [zero, one]
]


# Normalized matrices for the Fermat orbit of size 54.

M54_2 = [
    [zero, one + omega],
    [omega, -one]
]

M54_3 = [
    [-omega, omega],
    [zero, -one - omega]
]

M54_4 = [
    [one + omega, one],
    [zero, omega]
]


# Normalized matrices for the Fermat orbit of size 18.

M18_2 = [
    [zero, -one - omega],
    [-omega, one]
]

M18_3 = [
    [zero, one],
    [-one, one]
]

M18_4 = [
    [zero, omega],
    [one + omega, one]
]


# Compute determinant lists.

data_54 = determinant_data([Z, Id, M54_2, M54_3, M54_4])
data_18 = determinant_data([Z, Id, M18_2, M18_3, M18_4])


print("orbit 54 determinant list:")
print(data_54)
print("orbit 54 determinant set:")
print(sorted(set(data_54), key=lambda x: (x.a, x.b)))

print()

print("orbit 18 determinant list:")
print(data_18)
print("orbit 18 determinant set:")
print(sorted(set(data_18), key=lambda x: (x.a, x.b)))
