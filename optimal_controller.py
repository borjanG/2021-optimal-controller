#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 21:28:07 2021

@author: borjangeshkovski
"""

import numpy as np
from sympy import Matrix
from scipy.linalg import eig


def laplacian(n):
    """ Finite difference Dirichlet Laplacian (second derivative) in 1d.

    Arguments: 
        n: integer designating the number of interior discretization points
    
    Returns:
        A: list comprehension designating the finite difference second derivative
    """

    A = np.zeros((n, n))
    h = 1 / (n + 1)

    for i in range(n):
        for j in range(n):
            if i == j:
                A[i][j] = -2 / pow(h, 2)
            elif abs(i - j) == 1:
                A[i][j] = 1 / pow(h, 2)
    return A


def wave_operator(n_):
    """ Finite difference wave operator in 1d.

    Arguments: 
        n_: integer designating 2*(number of interior discretization points)
    
    Returns:
        A: list comprehension designating the finite difference wave operator [0 Id; Laplacian 0]
    """

    # The size of the wave matrix operator is n_=2*n.
    n = n_ // 2
    A = np.zeros((n_, n_))
    mat = laplacian(n)

    for i in range(n_):
        for j in range(n_):
            if i in list(range(n)) and j in list(range(n, n_)) and abs(i - j) == n:
                A[i][j] = 1
            elif i in list(range(n, n_)) and j in list(range(n)):
                A[i][j] = mat[i - n][j]
    return A


def convection_operator(n):
    """ Finite difference convection-diffusion operator.

    Arguments: 
        n: integer
    
    Returns: 

    """

    A_ = laplacian(n)
    h = 1 / (n + 1)

    A = 1 / (2 * h) * (np.zeros((n, n)) - np.diag(np.ones(n - 1), -1)
                       + np.diag(np.ones(n - 1), 1)) + A_

    return A


def p_(k, n, model):
    """ The matrices p_k(A) defined as
        p_k(A) = A^{n-k} + sum_{j=1}^{n-k} a_j*A^{n-k-j}
    if k<n-1, and 
        p_n(A) = Id. 
    
    Arguments: 
        k: integer denoting the column in P(b)
        n: integer denoting the dimension of the system, and thus the rank of P(b)
        model: a string key designating the model
    
    Returns:
        A numpy array designating the matrix p_k(A)
    """

    models = {"heat": laplacian,
              "wave": wave_operator,
              "convection": convection_operator}
    A = models[model](n)

    # We need to recover the coefficients of the characteristic
    # polynomial of A
    M = Matrix(A)
    a_ = M.charpoly().all_coeffs()

    # We aim to construct p_k(A) 
    __ = np.zeros((n, n))
    if k == n:
        return np.eye(n)
    else:
        for j in range(n - k):
            __ += float(a_[j + 1]) * np.linalg.matrix_power(A, n - k - (j + 1))
    return np.linalg.matrix_power(A, n - k) + __


def eigenvalue(b):
    """ The objective functional to be maximized: 
        lambda_1(P(b)P(b)^*)
    where the matrix 
        P(b)P(b)^* = sum_{k=1}^n p_k bb^* p_k^*
    with 
        p_k = A^{n-k} + sum_{j=1}^{n-k} a_j A^{n-k-j}
    whenever k<= n-1 and 
        p_k = Id 
    whenever k=n. Equivalently, we can minimize |P^-1(b)|^2_2.

    Arguments: 
        b: numpy array designating a controller vector.

    Returns: 
        A float designating the value of the eigenvalue
    """

    model = "convection"
    inverse = False

    n = len(b)
    if model == "heat" or model == "convection":
        mat = np.zeros((n, n))
        for k in range(n):
            mat[k] = np.matmul(p_(k + 1, n, model), b)
    else:
        # In the case of the wave equation, the full control operator
        # takes the form bb = [0 b].T, where b is of the same dimension
        # as the laplacian within A.
        n_wave = 2 * n
        mat = np.zeros((n_wave, n_wave))
        bb = np.zeros(n_wave)
        bb[n:] = b

        #  We construct the matrix P(b)
        for k in range(n_wave):
            mat[k] = np.matmul(p_(k + 1, n_wave, model), bb)

    if inverse:
        return np.linalg.norm(np.linalg.inv(mat.T).T, 2) ** 2
    ev, _ = eig(np.matmul(mat, mat.T))
    return -np.real(sorted(ev, key=float)[0])


def sphere(b):
    """ The map f(b) = sum(b[k]^2)
    """

    return sum([pow(b[k], 2) for k in range(len(b))])


if __name__ == "__main__":

    from scipy.optimize import minimize
    from scipy.optimize import dual_annealing, shgo, differential_evolution, NonlinearConstraint, Bounds
    import time

    METHODS = ["genetic_alg", "dual_annealing, trust-constr", "shgo"]
    METHOD = "genetic_alg"
    N = 3

    start = time.time()

    if METHOD == "genetic_alg":
        nlc = NonlinearConstraint(sphere, 1, 1)
        bounds = Bounds([-1.] * N, [1.] * N)
        result = differential_evolution(eigenvalue,
                                        bounds,
                                        maxiter=1000,
                                        constraints=nlc,
                                        x0=np.array([0.4, -0.35, 0.9]),
                                        seed=1)
        print(result)
        print(sphere(result.x))

    elif METHOD == "shgo":
        cons = ({'type': 'eq', 'fun': sphere})
        bounds = [(-1, 1), ] * N
        res = shgo(eigenvalue, bounds, iters=3, constraints=cons)
        print(res)

    elif METHOD == "dual_annealing":
        lw = [-1] * N
        up = [1] * N
        ret = dual_annealing(lambda x: eigenvalue(x) - 50 * (sphere(x) - 1),
                             bounds=list(zip(lw, up)))
        print(ret.x)
        print(ret.fun)

    else:
        b0 = np.random.rand(N)
        from scipy.optimize import NonlinearConstraint

        nonlinear_constraint = NonlinearConstraint(sphere, 1, 1)
        res = minimize(eigenvalue,
                       b0,
                       method='trust-constr',
                       constraints=nonlinear_constraint)
        print(res)
    end = time.time()
    print(f"Elapsed time: {end - start} seconds")