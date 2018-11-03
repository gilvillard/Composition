###############################
#  MINIMAL APPROXIMANT BASIS  #
###############################

# When order=1 (constant matrix): rely on kernel
def approximant_basis1(mat, shift):
    ##Input:
    #   * m x n matrix of field elements 'mat'
    #   * degree shift 'shift' (list of n integers)
    ##Output:
    #   * the shift-Popov approximant basis for 'mat' at order 1
    #         (working column-wise: right approximants)
    #   * the shift-column degree of this basis
    ##Reference:  [Giorgi-Jeannerod-Villard, ISSAC 2003, Algorithm M-Basis]
    n = mat.ncols()
    cdeg = copy(shift)
    # Find the permutation which stable-sorts the shift in nondecreasing order
    perm = Permutation(list(zip(*sorted( [(cdeg[i],i+1) for i in range(n)] ))[1]))
    # Permute the rows of mat accordingly
    matt = mat.with_permuted_columns(perm)
    ker = matt.right_kernel(basis="pivot").matrix().T
    pivind = []  # the kth element in this list will be the pivot index in column k of ker
    for j in range(ker.ncols()):
        i = ker.nrows()-1
        while i>=0 and ker[i,j] == 0:
            i -= 1
        pivind.append(i) # note that there must be a nonzero coefficient in ker[:,j]
    PolyRing.<x> = mat.base_ring()[]
    appbas = Matrix(PolyRing, n, n)
    for j in range(n):
        try:
            ind = pivind.index(j)
            appbas[:,j] = copy(ker[:,ind])
        except ValueError:
            appbas[j,j] = x
            cdeg[perm[j]-1] += 1

    # permute back the matrix
    appbas.permute_rows_and_columns(perm.inverse(),perm.inverse())
    return appbas,cdeg

# When order is small: apply mbasis1 iteratively
def approximant_basis_iter(pmat, order, shift):
    ##Input:
    #   * m x n univariate polynomial matrix 'mat'
    #   * approximation order 'order': a positive integer
    #   * degree shift 'shift' (list of n integers)
    ##Output:
    #   * a 'shift'-ordered weak Popov approximant basis for 'pmat'
    #   at order 'order'
    #   * the shift-column degree 'cdeg' of 'appbas'
    ##Reference:  [Giorgi-Jeannerod-Villard, ISSAC 2003, Algorithm M-Basis]
    n = pmat.ncols()
    cdeg = copy(shift)
    (pR,x) = pmat.base_ring().objgen()
    appbas = Matrix.identity(pR, n, n)
    residual = copy(pmat)
    for k in range(order):
        appbas1,cdeg = approximant_basis1(residual(0),cdeg)
        appbas = appbas * appbas1
        residual = residual * appbas1
        matrix_shift(residual, -1)
        matrix_truncate(residual, order-k-1)
    return appbas,cdeg


# When order is large: divide-and-conquer via polynomial matrix multiplication
def approximant_basis_dnc(pmat, order, shift):
    ##Input:
    #   * m x n univariate polynomial matrix 'pmat'
    #   * approximation order 'order': a positive integer
    #   * degree shift 'shift' (list of n integers)
    ##Output:
    #   * a shift-ordered weak Popov approximant basis 'appbas'
    #     for pmat at order 'order'
    #   * the shift-column degree 'cdeg' of 'appbas'
    # Reference: [Giorgi-Jeannerod-Villard, ISSAC 2003, Algorithm PM-Basis]
    if order<=32:
        return approximant_basis_iter(pmat,order,shift)
    else:
        order1 = order//2
        order2 = order - order1
        # first recursive call
        residual = copy(pmat)
        matrix_truncate(residual, order1)
        appbas1,cdeg = approximant_basis_dnc(residual, order1, shift)
        # residual = (x^{-order1} * pmat * appbas1) mod x^order2
        residual = pmat*appbas1
        matrix_shift(residual, -order1)
        matrix_truncate(residual, order2)
        # second recursive call
        appbas2,cdeg = approximant_basis_dnc(residual, order2, cdeg)
        return appbas1*appbas2, cdeg

# Popov approximant basis for uniform shift
def approximant_basis(pmat, order, shift=None):
    ##Input:
    #   * m x n univariate polynomial matrix 'pmat'
    #   * approximation order 'order': a positive integer
    #   * degree shift 'shift' (list of n integers)
    #     (shift==None is interpreted as uniform shift)
    ##Output:
    #   * the shift-Popov approximant basis 'appbas' for pmat at order 'order'
    n = pmat.ncols()
    if shift==None:
        shift=[0]*n
    # compute approximant basis for uniform shift
    appbas,cdeg = approximant_basis_dnc(pmat, order, shift)
    # deduce pivot degree of Popov form and use them as shifts
    popov_shift = [shift[i]-cdeg[i] for i in range(n)]
    # compute approximant basis for this "Popov shift"
    appbas,cdeg = approximant_basis_dnc(pmat, order, popov_shift)
    # right-multiply by inverse of leading matrix
    lmat = appbas.leading_matrix(shifts=popov_shift, row_wise=False)
    appbas = appbas * lmat.inverse()
    return appbas
