# Throw any square matrix of any dimensions at the program here:
matrix = [
    [12, 13, 32, 64, 11, 23, 32],
    [12, 34, 23, 43, 23, 12, 21],
    [16, 72, 79, 86, 94, 75, 42],
    [70, 68, 79, 21, 25, 27, 33],
    [37, 45, 65, 50, 11, 29, 15],
    [25, 35, 11, 10, 19, 92, 14],
    [69, 98, 41, 23, 54, 33, 12],
] # det(matrix) = -416561511505

# Use the (1, 1)*(2, 2) - (0, 1)*(1, 0) formula to calculate a 2x2 matrix
# Then multiplies that with the entry it is a minor of, which is an essential part of the algorithm used
def two_x_two_determinant(matrix):
    """ Calculates the determinant of a 2x2 matrix and multplies it with the entry it is a minor of """
    return (matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]) * matrix[-1]


# is_modified is a boolean value
# is_modified is set to False if the matrix is an usual 2D list of entries
# is_modified is set to True when the last element of the list (matrix[-1]) is the entry of the original matrix that the matrix passed into the function is a minor of
def dimensions(matrix, is_modified):
    """ Calculates the dimension of the matrix passed as an argument to the function """
    rows = len(matrix)-1 if is_modified else len(matrix)    # omits the last entry in list when is_modified is True
    columns = len(matrix[0])
    return rows, columns


# Calculate the cofactor and return a number which is the cofactor of entry (i, j) times the entry itself
# The cofactor is multiplied with the entry here because of how the algorithm of this program works
# Which is original but not necessarily effective or fast
# Also the algorithm is slightly different than the algorithm used in traditional manual matrix solving on paper
def cofactor_multiple(matrix, is_modified, i, j):
    """ Returns the cofactor of the entry, and the entry of matrix at (i, j) as last element """
    # if the matrix is a modified list, the entry at (i, j) needs the multiplied with already existing last element of the provided matrix
    factor = matrix[i][j] if not is_modified else matrix[i][j] * matrix[-1]
    # tradional indices start at (1, 1)
    i += 1
    j += 1
    return ((-1) ** (i+j)) * factor


# Calculate not only the minor, but also the entry*cofactor using the cofactor_multiple function
# This entry*cofactor is returned as the last element of a list while the elements preceding are essentially the rows of the required minor
def minor_multiple(matrix, is_modified, i, j):
    """ Returns the minor, and the cofactor_multiple of the entry at index (i, j) as last element """
    minor_matrix = []
    row_indices = range(len(matrix)-1) if is_modified else range(len(matrix))    # omitting last element if modified
    # creating a copy of the argument "matirx" and breaking the alias
    for row_index in row_indices:
        row = matrix[row_index]
        minor_row = []
        for entry in row:
            minor_row.append(entry)
        minor_matrix.append(minor_row)
    
    if is_modified:
        minor_matrix.append(matrix[-1])

    # deleting all entries from the row and column that pass through entry (i, j)
    del minor_matrix[i]    
    for l in range(len(minor_matrix)):
        if not is_modified:
            del minor_matrix[l][j]

        else:
            if l != len(minor_matrix)-1:
                del minor_matrix[l][j]

    # setting the last element of list to be the entry the returned matrix is a minor of
    if is_modified:
        minor_matrix[-1] = cofactor_multiple(matrix, is_modified, i, j)
        return minor_matrix
    
    else:
        return minor_matrix + [cofactor_multiple(matrix, is_modified, i, j), ]


# Uses the "original algorithm" here
def det(matrix):
    """ Calculates the determinant of provided matrix """
    modification = False
    minor_list = [matrix]    # list of all the modified 2x2 minor matrices of the given matrix
    # keep looping until all minors are 2x2
    while dimensions(minor_list[0], modification)[0] > 2:
        if len(minor_list) == 1:
            minor_list.pop()
            for j_index in range(len(matrix[0])):
                # find the minor of row 0 and column j of matrix[0], and put it inside minor_list
                minor_list.append(minor_multiple(matrix, False, 0, j_index))

            modification = True
        
        else:
            for minor_matrix_index in range(len(minor_list)):
                minor_matrix = minor_list[minor_matrix_index]
                broken_minor_matrices = []    # list of further broken down matrices
                for j_index in range(len(minor_matrix[0])):
                    # find the minor of row 0 and column j of matrix_matrix[0], and put it inside broken_minor_matrices
                    broken_minor_matrix = minor_multiple(minor_matrix, True, 0, j_index)
                    broken_minor_matrices.append(broken_minor_matrix)

                # replacing each bigger minor with smaller minors (where the "smaller" minors are the minors of the entries in the 1st row of the bigger matrix)
                minor_list[minor_matrix_index] = broken_minor_matrices[0]
                for broken_minor_matrix_index in range(1, len(broken_minor_matrices)):
                    minor_list.append(broken_minor_matrices[broken_minor_matrix_index])

    determinant = 0
    # find determinant of each 2x2 minor_multiple and take summation, which returns our desired determinant
    for minor_matrix in minor_list:
        determinant += two_x_two_determinant(minor_matrix)

    return determinant


# Calculate det(matrix)
determinant = det(matrix)


# Print the provided matrix and its determinant
for row_index in range(len(matrix)):
    row = matrix[row_index]
    midline = (dimensions(matrix, 0)[0]-1) // 2
    print("    |" if row_index != midline else "det |", end=" ")
    for entry in row:
        print(entry, end=" ")
        
    print("|" if row_index != midline else "| = "+str(determinant))
        
