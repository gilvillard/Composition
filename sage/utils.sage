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

##################
#  HERMITE FORM  #
##################

def hermite_form(mat, transformation=False):
    # Returns the column-wise, upper triangular Hermite form of 'mat'.
    # If "transformation=True", also return the corresponding unimodular
    # transformation.
    # Unfortunately we cannot directly rely on the Hermite form implementation
    # of SageMath; the latter computes the row-wise, upper triangular form.
    # Taking transposed matrices and row/column permutations does not give the
    # desired form.
    # Instead we use shifted column basis reduction with a "Hermite" shift,
    # obtaining a triangular matrix but non-normalized. This is enough to give
    # us the degrees of the diagonal entries of the Hermite form, from which
    # we can obtain the Hermite form by another shifted column basis reduction.
    # Reference: [Labahn - Neiger - Zhou, J. Complexity, 2017]

    dim = mat.nrows()

    # 1. Compute shifted reduced form with the Hermite shift
    # -->  matrix column-wise equivalent to 'mat', which is upper triangular up
    # to column permutation
    hermite_shift = [mat.degree() * dim * i for i in range(dim)]
    hnf = mat.reduced_form(shifts=hermite_shift,row_wise=False)

    # 2. Deduce the degrees of the diagonal entries of the Hermite form
    pivot_index,pivot_degree = hnf.leading_positions(return_degree=True,shifts=hermite_shift,row_wise=False)
    # permute these degrees in the right order
    # (the pivot index of the Hermite form is [0,...,dim-1])
    pivot_degree = list(zip(*sorted([(pivot_index[i],pivot_degree[i]) for i in range(dim)]))[1])

    # 3. Run a new column basis reduction using these degrees in the shift
    shift = [mat.degree()*dim - pivot_degree[i] for i in range(dim)]
    # FIXME: here we use transpose instead of row_wise=false because of a bug
    # in Sage 8.3 when asking to return the transformation
    hnf,trans = mat.transpose().reduced_form(shifts=shift,transformation=True)
    hnf = hnf.transpose()
    trans = trans.transpose()

    # 4. Multiply by the inverse of the leading matrix of the obtained form
    inv_lmat = hnf.leading_matrix(shifts=shift,row_wise=False).inverse()
    hnf = hnf * inv_lmat
    trans = trans * inv_lmat

    if transformation:
        return hnf,trans
    else:
        return hnf

load('./subroutines/approximant-basis.sage')
load('./subroutines/system-solve.sage')
load('./subroutines/matrix-division-with-remainder.sage')
