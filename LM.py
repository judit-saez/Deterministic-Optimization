# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 14:05:14 2021

@author: judit
"""


#### MatrixTools ####
def transpose(matrix):
    """
    Returns a transpose of a matrix.
        :param matrix: The matrix to be transposed
 
        :return: The transpose of the given matrix
    """
    Tmatrix = [[matrix[j][i] for j in range(len(matrix))] 
               for i in range(len(matrix[0]))]
    return Tmatrix
 
def zeros_matrix(rows, cols):
    """
    Creates a matrix filled with zeros.
        :param rows: the number of rows the matrix should have
        :param cols: the number of columns the matrix should have
 
        :return: list of lists that form the matrix
    """
    M = []
    while len(M) < rows:
        M.append([])
        while len(M[-1]) < cols:
            M[-1].append(0.0)
 
    return M
    
def matrix_multiply(A, B):
    """
    Returns the product of the matrix A * B
        :param A: The first matrix - ORDER MATTERS!
        :param B: The second matrix
 
        :return: The product of the two matrices
    """
    #Ensure A & B dimensions are correct for multiplication
    rowsA = len(A)
    colsA = len(A[0])
    rowsB = len(B)
    colsB = len(B[0])
    if colsA != rowsB:
        raise ArithmeticError(
            'Number of A columns must equal number of B rows.')
 
    #Store matrix multiplication in a new matrix
    C = zeros_matrix(rowsA, colsB)
    for i in range(rowsA):
        for j in range(colsB):
            total = 0
            for ii in range(colsA):
                total += A[i][ii] * B[ii][j]
            C[i][j] = total
 
    return C
    
def multiply_matrices(list):
    """
    Find the product of a list of matrices from first to last
        :param list: The list of matrices IN ORDER
 
        :return: The product of the matrices
    """
    #Start matrix product using 1st matrix in list
    matrix_product = list[0]
 
    #Loop thru list to create product
    for matrix in list[1:]:
        matrix_product = matrix_multiply(matrix_product, matrix)
 
    return matrix_product


def as_scalar(one_element_matrix):
    """
    Converts [[num]] to scalar: num

        :param one_element_matrix: object of type [[num]]
 
        :return: num as a number
    """
    result = one_element_matrix[0][0]
    return result

def scalar_times_matrix(scalar, A):
    """
    Multiply a matrix by a scalar
        :param scalar: any number (float, int...)
        :param A: A matrix
    
        :return: The resultant matrix once A has been multiplied by the scalar
    """
    rowsA = len(A)
    colsA = len(A[0])
    C = zeros_matrix(rowsA, colsA)
    for i in range(rowsA):
        for j in range(colsA):
           C[i][j] = scalar * A[i][j] 
    return C

def matrix_addition(A, B):
    """
    Adds two matrices and returns the sum
        :param A: The first matrix
        :param B: The second matrix
 
        :return: Matrix sum
    """
    #Ensure dimensions are valid for matrix addition
    rowsA = len(A)
    colsA = len(A[0])
    rowsB = len(B)
    colsB = len(B[0])
    if rowsA != rowsB or colsA != colsB:
        raise ArithmeticError('Matrices are NOT the same size.')
 
    #Create a new matrix for the matrix sum
    C = zeros_matrix(rowsA, colsB)
 
    #Perform element by element sum
    for i in range(rowsA):
        for j in range(colsB):
            C[i][j] = A[i][j] + B[i][j]
 
    return C

def matrix_subtraction(A, B):
    """
    Subtracts matrix B from matrix A and returns difference
        :param A: The first matrix
        :param B: The second matrix
 
        :return: Matrix difference
    """
    #Ensure dimensions are valid for matrix subtraction
    rowsA = len(A)
    colsA = len(A[0])
    rowsB = len(B)
    colsB = len(B[0])
    if rowsA != rowsB or colsA != colsB:
        raise ArithmeticError('Matrices are NOT the same size.')
 
    #Create a new matrix for the matrix difference
    C = zeros_matrix(rowsA, colsB)
 
    #Perform element by element subtraction
    for i in range(rowsA):
        for j in range(colsB):
            C[i][j] = A[i][j] - B[i][j]
 
    return C



#### Functions to get inverse of a matrix ####
def get_matrix_minor(m,i,j):
    """
    :return: minor of a matrix
    
    """
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def get_matrix_determinant(m):
    """
    :return: determinant of a 2x2 matrix
    
    """
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*get_matrix_determinant
        (get_matrix_minor(m,0,c))
    return determinant

def get_matrix_inverse(m):
    """
    :return: inverse of a 2x2 matrix
    
    """
    determinant = get_matrix_determinant(m)
    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = get_matrix_minor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * get_matrix_determinant(minor))
        cofactors.append(cofactorRow)
    cofactors = transpose(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors


#### Method to diagonalize a matrix ####
def jacobi_method(A):
    """
    :return: Diagonalized matrix of a 2x2 symmetric matrix
    
    """
    #special case for 2x2 symmetric matrix:
    a12 = A[0][1]
    a11 = A[0][0]
    a22 = A[1][1]
    alpha = (a22 - a11)/(2*a12)
    if alpha < 0:
        sign_alpha = -1
    else:
        sign_alpha = 1
        
    tan_theta = -alpha+sign_alpha*(alpha**2+1)**(1/2)
    Adiag = zeros_matrix(2, 2)
    Adiag[0][0] = a11 - tan_theta*a12
    Adiag[1][1] = a22 + tan_theta*a12
    return Adiag


#### LMA Functions ####
def evalJ(x1):
    """
    Evaluates the Jacobian matrix of the Rosenbrock function at point (x1, x2)
        :param x1: The first coordinate of the point 
 
        :return: The evaluated Jacobian matrix
    """
    J = [[-1,        0],
         [-10*x1*2, 10]]
    return J

def evalF(x1,x2):
    """
    Evaluates the partial functions of the Rosenbrock's at point (x1, x2)
        :param x1: The first coordinate of the point 
        :param x2: The second coordinate of the point
        
        :return: The evaluated matrix
    """
    F = [[1-x1],
         [10*(x2-x1**2)]]
    return F

def get_max_diag_A(A):
    """
        :param A: The matrix
        
        :return: the maximum element in the diagonal
    """
    rowsA = len(A)
    colsA = len(A[0])
    diag_list = []
    for i in range(rowsA):
        for j in range(colsA):
            if i==j:
                diag_list.append(A[i][j])
    max_diag_A = max(diag_list)
    return max_diag_A

#### Auxiliar Functions ####
def get_length(vector):
    """
        :param vector: The vector
        
        :return: the length of this vector
    """
    length = ((vector[0][0])**2+(vector[1][0])**2)**(1/2)
    return length

def get_h_lm(g, A, mu, Adiag):
    """
    Computes the parameter h_lm

    """
    identity = [[1,0],[0,1]]
    matrix = matrix_addition(A, scalar_times_matrix(mu, Adiag))
    h_lm = multiply_matrices([get_matrix_inverse(matrix), scalar_times_matrix((-1), g)])
    return h_lm

def evalRosenbrock(x1,x2):
    """
    Evaluates the Rosenbrock function at point (x1, x2)
        :param x1: The first coordinate of the point 
        :param x2: The second coordinate of the point
 
        :return: The result of the Rosenbrock function
    """
    f = 100*(x2 - x1*x1)*(x2 - x1*x1) + (1 - x1)*(1 - x1)
    return f

def compute_ro(x, x_new, h_lm, mu, g):
    """
    Computes the parameter h_lm

    """
    ro_num = evalRosenbrock(x[0][0], x[1][0]) - evalRosenbrock(x_new[0][0], x_new[1][0])
    h_lm_T = transpose(h_lm)
    mu_h_lm = scalar_times_matrix(mu, h_lm)
    difference = matrix_subtraction(mu_h_lm, g)
    ro_den = 1/2*as_scalar(multiply_matrices([h_lm_T, difference]))
    ro = ro_num/ro_den
    return ro



#### Levenberg-Marquardt algorithm ####
#### Main ####
# Initialize 
x = [[-1.5], [-1]] 
print("INFORMATION: \n")
print("Initial Point:")
print("Coord. x1: ", x[0][0], ", Coord. x2: ",  x[1][0], "\n")

J  = evalJ(x[0][0])
JT = transpose(J)
F = evalF(x[0][0], x[1][0])
A = multiply_matrices([JT,J])
g = multiply_matrices([JT, F])
tau = 0.0001 #si no funciona canviar la tau

max_diag_A = get_max_diag_A(A)
mu = tau * max_diag_A
epsilon = 0.0001
nu = 2
k = 0
kmax = 100

# LMA Main Loop
while not get_length(g) <= epsilon and k < kmax:
    k = k + 1
    Adiag = jacobi_method(A)
    h_lm = get_h_lm(g, A, mu, Adiag)
    if get_length(h_lm) <= epsilon*(get_length(x)+epsilon):
        break
    else:
        x_new = matrix_addition(x, h_lm)
        ro = compute_ro(x, x_new, h_lm, mu, g)
        if ro > 0:
            x = x_new
            A = multiply_matrices([transpose(evalJ(x[0][0])), evalJ(x[0][0])])
            g = multiply_matrices([transpose(evalJ(x[0][0])), evalF(x[0][0], x[1][0])])
            mu = mu * max([1/3, 1-(2*ro-1)**3])
            nu = 2
        else:
            mu = mu * nu
            nu = 2 * nu

            
    #Uncomment this line to see x at every iteration 
    #print(x)
       
    #Uncomment this line to see the value of the Rosenbrock functions at every 
    #iteration (its value should gradually decrease)
    #print(evalRosenbrock(x[0][0], x[1][0]))
    
    
print("Final Point:")    
print("Coord. x1: ", "{0:.7f}".format(x[0][0]), ", Coord. x2: ",  
      "{0:.7f}".format(x[1][0]), "\n")


print("Number of iterations:") 
print(k)
    

















