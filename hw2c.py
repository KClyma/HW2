from copy import deepcopy
from NumericalMethods import GaussSeidel


def main():
    # Example system 1
    Aaug1 = [
        [3, 1, -1, 2],
        [1, 4, 1, 12],
        [2, 1, 2, 10]
    ]
    x0_1 = [0, 0, 0]
    Niter = 15

    solution1 = GaussSeidel(Aaug1, x0_1, Niter)
    print("Solution for system 1:", solution1)

    # Example system 2
    Aaug2 = [
        [1, -10, 2, 4, 2],
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]
    ]
    x0_2 = [0, 0, 0, 0]

    solution2 = GaussSeidel(Aaug2, x0_2, Niter)
    print("Solution for system 2:", solution2)

if __name__ == "__main__":
    main()