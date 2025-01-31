from copy import deepcopy
from NumericalMethods import GaussSeidel
import numpy as np


def main():
    # Augmented matrix [A|b] for the system
    Aaug_initial = np.array([[4, -1, 0, 0, 15],
                            [-1, 4, -1, 0, 10],
                            [0, -1, 4, -1, 10],
                            [0, 0, -1, 3, 10]], dtype=float)

    # Initial guess (could also be np.zeros, or some other guess)
    x0 = np.zeros(Aaug_initial.shape[0])

    # Solve the system using Gauss-Seidel method
    solution = GaussSeidel(Aaug_initial, x0)

    # Print the solution
    return "Solution:", solution

    pass

if __name__=="__main":
    main()
