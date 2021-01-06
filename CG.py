# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 17:04:00 2020

@author: judit
"""




#### MatrixTools #############################################################
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

  
    



#### CGD Functions ###########################################################
def evalA(x1, x2):
    """
    Evaluates the Hessian matrix of the Rosenbrock function at point (x1, x2)
        :param x1: The first coordinate of the point 
        :param x2: The second coordinate of the point
 
        :return: The evaluated Hessian matrix
    """
    A = [[1200**x1-400*x2+2, -400*x1],
         [-400*x1,            200   ]]
    return A

def evalG(x1, x2):
    """
    Evaluates the gradient vector of the Rosenbrock function at point (x1, x2)
        :param x1: The first coordinate of the point 
        :param x2: The second coordinate of the point
 
        :return: The evaluated gradient vector
    """
    grad = [[(400*x1*x1*x1 - 400*x1*x2 + 2*x1 - 2)],
                  [200*(x2-x1**2)]]
    return grad

def evalRosenbrock(x1,x2):
    """
    Evaluates the Rosenbrock function at point (x1, x2)
        :param x1: The first coordinate of the point 
        :param x2: The second coordinate of the point
 
        :return: The result of the Rosenbrock function
    """
    f = 100*(x2 - x1*x1)*(x2 - x1*x1) + (1 - x1)*(1 - x1)
    return f


#### Auxiliar Functions ######################################################
def compute_alpha(x, d):
    """
 
    :return: alpha value

    """
    alpha_num = as_scalar(multiply_matrices([transpose(evalG(x[0][0], 
                                                             x[1][0])), d]))
    alpha_den = as_scalar(multiply_matrices([transpose(d), 
                                             evalA(x[0][0], x[1][0]), d]))
    alpha = alpha_num/alpha_den*(-1)

    return alpha
    
def compute_beta(r, r_old, delta_old):
    """
 
    :return: beta value

    """
    beta_num = as_scalar(multiply_matrices([transpose(r), 
                                            matrix_subtraction(r, r_old)]))
    beta_den = as_scalar(delta_old)
    beta = beta_num/beta_den
    return beta





#### Main ####################################################################
# Initialize 
i = 0
x = [[-1.5], [-1]] 
print("INFORMATION: \n")
print("Initial Point:")
print("Coord. x1: ", x[0][0], ", Coord. x2: ",  x[1][0], "\n")


r = scalar_times_matrix((-1), evalG(x[0][0], x[1][0]))
d = r
delta_new = multiply_matrices([transpose(r), r])
delta0 = delta_new


i_max = 1000
tol = 0.0001

# CGD Main Loop
while i < i_max and abs(as_scalar(delta_new)) > abs(tol*tol):
    alpha = compute_alpha(x, d)
    x = matrix_addition(x, scalar_times_matrix(alpha, d))
    r_old = r
    r = evalG(x[0][0], x[1][0])
    delta_old = delta_new
    delta_new = multiply_matrices([transpose(r), r])
    beta = compute_beta(r, r_old, delta_old)
    if beta > 0:
        beta = beta
    elif 0 > beta:
        beta = 0
    d = matrix_addition(r, scalar_times_matrix(beta, d))
    i = i + 1 
    
    #Uncomment this line to see x at every iteration 
    #print(x)
    
    #Uncomment this line to see the gradient at every iteration
    #print(r)
    
    #Uncomment this line to see the value of the Rosenbrock functions at every 
    #iteration (its value should gradually decrease)
    #print(evalRosenbrock(x[0][0], x[1][0]))
    
    
print("Final Point:")    
print("Coord. x1: ", "{0:.7f}".format(x[0][0]), ", Coord. x2: ",  
      "{0:.7f}".format(x[1][0]), "\n")

print("Final Gradient:")
print("grad.  x1: ", "{0:.7f}".format(r[0][0]), ", grad.  x2: ",
      "{0:.7f}".format(r[1][0]), "\n")

print("Number of iterations:") 
print(i)
    

















