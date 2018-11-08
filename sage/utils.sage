######################
#  SOME BASIC TOOLS  #
######################
def print_separator():
    print "========"

def matrix_truncate(A, d):
    # input: polynomial matrix A and truncation order d
    # action: A = A mod X^d
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            A[i,j] = A[i,j].truncate(d)

def matrix_shift(A, d):
    # input: polynomial matrix A and degree d
    # action: A = A * X^d if d>=0, A//X^-d if d<0
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            A[i,j] = A[i,j].shift(d)

load('./subroutines/approximant-basis.sage')
load('./subroutines/matrix-division-with-remainder.sage')
