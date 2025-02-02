from copy import deepcopy
from NumericalMethods import GaussSeidel
from Gauss_Elim import MakeDiagDom
import warnings


def is_diagonally_dominant_by_row(matrix):
    """Check if the matrix is diagonally dominant by rows."""
    num_rows = len(matrix)
    num_cols = len(matrix[0])

    for i in range(num_rows):  # Iterate over rows
        row_sum = 0
        for j in range(num_cols):  # Iterate over columns
            if i != j:  # Skip the diagonal element
                row_sum += abs(matrix[i][j])
        if abs(matrix[i][i]) < row_sum:
            return False
    return True


def is_diagonally_dominant_by_column(matrix):
    """Check if the matrix is diagonally dominant by columns."""
    num_rows = len(matrix)
    num_cols = len(matrix[0])

    if num_rows != num_cols:
        raise ValueError("Matrix is not square, cannot check diagonal dominance by columns.")

    for j in range(num_cols):  # Iterate over columns
        column_sum = 0
        for i in range(num_rows):  # Iterate over rows
            if i != j:  # Skip the diagonal element
                column_sum += abs(matrix[i][j])
        if abs(matrix[j][j]) < column_sum:
            return False
    return True


def is_diagonally_dominant(matrix):
    """Check if the matrix is diagonally dominant by both row and column."""
    # Check diagonal dominance by row
    if not is_diagonally_dominant_by_row(matrix):
        return False

    # Check diagonal dominance by column
    try:
        if not is_diagonally_dominant_by_column(matrix):
            return False
    except ValueError:
        print("Error: Matrix is not square.")
        return False

    return True


def main():
    # Suppress warnings for specific non-critical diagonal dominance warnings
    warnings.filterwarnings("ignore", message="Column .* could not be made diagonally dominant.")

    # System 1
    Aaug1 = [
        [3, 1, -1, 2],
        [1, 4, 1, 12],
        [2, 1, 2, 10]
    ]
    # Ensure diagonal dominance (but allow for non-strict diagonal dominance)
    Aaug1 = MakeDiagDom(Aaug1)

    # Check diagonal dominance by both row and column
    if is_diagonally_dominant(Aaug1):
        print("System 1 is diagonally dominant by both row and column.")
    else:
        print("System 1 is NOT diagonally dominant by both row and column.")

    x0_1 = [0, 0, 0]
    Niter = 15

    solution1 = GaussSeidel(Aaug1, x0_1, Niter)
    print("\nSolution for system 1:")
    for num in solution1:
        print(f"{num:.5f}")

    # System 2
    Aaug2 = [
        [1, -10, 2, 4, 2],
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]
    ]
    # Ensure diagonal dominance (but allow for non-strict diagonal dominance)
    Aaug2 = MakeDiagDom(Aaug2)

    # Check diagonal dominance by both row and column
    if is_diagonally_dominant(Aaug2):
        print("\nSystem 2 is diagonally dominant by both row or column.")
    else:
        print("\nSystem 2 is NOT diagonally dominant by both row or column.")

    x0_2 = [0, 0, 0, 0]
    Niter = 15

    solution2 = GaussSeidel(Aaug2, x0_2, Niter)
    print("\nSolution for system 2:")
    for num in solution2:
        print(f"{num:.5e}")

if __name__ == "__main__":
    main()