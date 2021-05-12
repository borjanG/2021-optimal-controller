import numpy as np
from scipy.linalg import eig
from scipy.optimize import minimize
from scipy.sparse import rand

def laplacian(n):
    """
    Finite difference laplacian in 1d.
    """
    A = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            if i==j: 
                A[i][j] = -2/pow(n-1,2)
            elif abs(i-j)==1: 
                A[i][j] = 1/pow(n-1,2)
    return A

def wave(n_): 
    """
    Finite difference wave operator in 1d.
    """

    # The parameter n is the size of the laplacian within the wave 
    # operator [0 Id; Delta 0].
    # The size of the wave matrix operator is n_=2*n
    
    n = n_//2
    A = np.zeros((n_, n_))
    mat = laplacian(n)
    
    for i in range(n_):
        for j in range(n_):
            if i in list(range(n)) and j in list(range(n, n_)) and abs(i-j)==n:
                A[i][j] = 1 
            elif i in list(range(n, n_)) and j in list(range(n)):
                A[i][j] = mat[i-n][j]
    return A

def p_(k, n, flag): 
    """
    The matrices p_k(A).
    """

    # We distinguish the laplacian vs wave
    if flag: 
        A = laplacian(n)
    else: 
        # Ici, n = n_wave. Donc c'est la dimension du systeme, donc
        # 2*dimension(b)
        A = wave(n)

    # We need to recover the coefficients of the characteristic
    # polynomial of A
    from sympy import Matrix
    M = Matrix(A)
    a_ = M.charpoly().all_coeffs()

    # We aim to construct p_k(A), which reads
    # p_k(A) = A^{n-k} + \sum_{j=1}^{n-k} a_j*A^{n-k-j} if k<n-1
    # and p_n(A) = Id. 
    __ = np.zeros((n, n))
    if k == n: 
        return np.eye(n)
    else: 
        for j in range(n-k): 
            __ += float(a_[j+1])*np.linalg.matrix_power(A, n-k-(j+1))
    return np.linalg.matrix_power(A, n-k) + __

def eigenvalue(b):
    """
    The objective functional to be maximized: \lambda_1(P(b)P(b)^*) 
    where the matrix P(b)P(b)^* = \sum_{k=1}^n p_k bb^* p_k^*
    with p_k = A^{n-k} + \sum_{j=1}^{n-k} a_j A^{n-k-j} whenever k<= n-1 
    and p_k = Id whenever k=n.
    
    Equivalently, we can minimize |P^-1(b)|^2_2.
    """
    
    # Flag = 0 -> wave; else heat.
    flag = 1

    n = len(b)
    if flag:
        mat = np.zeros((n, n))
        for k in range(n): 
            mat[k] = np.matmul(p_(k+1, n, flag), b)
    else:
        # In the case of the wave equation, the full control operator
        # bb takes the form bb = [0 b].T, where b is of the same dimension
        # as the laplacian within A.
        n_wave = 2*n
        mat = np.zeros((n_wave, n_wave))
        bb = np.zeros(n_wave)
        bb[n:] = b

        # We construct the matrix P(b)
        for k in range(n_wave): 
            mat[k] = np.matmul(p_(k+1, n_wave, flag), bb)

    #ev, ew = eig(np.matmul(mat.T, mat))    
    #return -sorted(ev, key=float)[0]

    return np.linalg.norm(np.linalg.inv(mat.T).T, 2)**2

def sphere(b):
    return sum([pow(b[k], 2) for k in range(len(b))])

""" def kalman(b):
    n = len(b)
    A = laplacian(n)
    M = np.zeros((n,n))
    for k in range(n):
        M[k] = np.matmul(np.linalg.matrix_power(A, k), b)
    return np.linalg.matrix_rank(M.T) """

n = 5
b0 = np.random.rand(n)

from scipy.optimize import NonlinearConstraint
nonlinear_constraint = NonlinearConstraint(sphere, 1, 1)
##nonlinear_constraint = [NonlinearConstraint(sphere, 1, 1), NonlinearConstraint(kalman, n, n)]

res = minimize(eigenvalue, b0, method='trust-constr', constraints=nonlinear_constraint)
print(res)
##print('Initial guess = ', b0)