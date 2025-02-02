#region imports
import Gauss_Elim as GE  # this is the module from lecture 2 that has usefule matrix manipulation functions
from math import sqrt, pi, exp, cos
from Gauss_Elim import MakeDiagDom
#endregion

#region function definitions
def Probability(PDF, args, c, GT=True):
    """
    This is the function to calculate the probability that x is >c or <c depending
    on the GT boolean.
    Step 1:  unpack args into mu and stDev
    Step 2:  compute lhl and rhl for Simpson
    Step 3:  package new tuple args1=(mu, stDev, lhl, rhl) to be passed to Simpson
    Step 4:  call Simpson with GNPDF and args1
    Step 5:  return probability
    :param PDF: the probability density function to be integrated
    :param args: a tuple with (mean, standard deviation)
    :param c: value for which we ask the probability question
    :param GT: boolean deciding if we want probability x>c (True) or x<c (False)
    :return: probability value
    """
    mu, sig = args
    if GT:
        lhl, rhl = c, mu + 5 * sig #integrate from c to upper bound
    else:
        lhl,rhl = mu -5 * sig, c #Integrate from lower bound to c

    args1 = (mu, sig, lhl, rhl)
    probability = Simpson(PDF, args1)
    return probability

def GPDF(args):
    """
    Here is where I will define the Gaussian probability density function.
    This requires knowing the population mean and standard deviation.
    To compute the GPDF at any value of x, I just need to compute as stated
    in the homework assignment.
    Step 1:  unpack the args tuple into variables called: x, mu, stDev
    Step 2:  compute GPDF value at x
    Step 3:  return value
    :param args: (x, mean, standard deviation)  tuple in that order
    :return: value of GPDF at the desired x
    """
    # Step 1: unpack args
    x, mu, sig = args
    # step 2: compute GPDF at x
    fx = (1 / (sig * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sig) ** 2)
    # step 3: return value
    return fx

def Simpson(fn, args, N=100):
    """
    This executes the Simpson 1/3 rule for numerical integration (see page 832, Table 19.4).
    As I recall:
    1. divide the range from x=lhl to x=rhl into an even number of parts. Perhaps 20?
    2. compute fx at each x value between lhl and rhl
    3. sum the even and odd values of fx as prescribed
    4. return the area beneath the function fx
    :param fx: some function of x to integrate
    :param args: a tuple containing (mean, stDev, lhl, rhl)
    :return: the area beneath the function between lhl and rhl
    """
    mu, sig, lhl, rhl =args
    if N % 2 != 0:
        N += 1 #ensure N is even


    h = (rhl - lhl) / N
    x_values = [lhl + i * h for i in range(N + 1)]
    fx_values = [fn((x, mu, sig)) for x in x_values]

    area = fx_values[0] + fx_values[-1]
    area += 4 * sum(fx_values[1:N:2]) #odd terms
    area += 2 * sum(fx_values[2:N-1:2]) #even terms

    area *= h / 3
    return area

def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    This funciton implements th Secant method to find the root of an equation.  You should write your equation in a form
    fcn = 0 such that when the correct value of x is selected, the fcn actually equals zero (or very close to it).
    :param fcn: the function for which we want to find the root
    :param x0: x value in neighborhood of root (or guess 1)
    :param x1: another x value in neighborhood of root (or guess x0+1)
    :param maxiter: exit if the number of iterations (new x values) equals this number
    :param xtol:  exit if the |xnewest - xprevious| < xtol
    :return: tuple with: (the final estimate of the root (most recent value of x), number of iterations)
    """
    for i in range(maxiter):
        f_x0, f_x1 = fcn(x0), fcn(x1)

        # Check if the difference is too small to avoid division by zero
        if abs(f_x1 - f_x0) < 1e-12:  # Avoid division by an almost zero number
            raise ValueError(
                "Small difference between function values, leading to possible division by zero in Secant method.")

        # Secant method update
        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)

        # Check if the estimated root has converged
        if abs(x2 - x1) < xtol:
            return x2, i + 1  # Return root and iterations

        # Update for the next iteration
        x0, x1 = x1, x2

    return x2, maxiter  # Return after max iterations


def GaussSeidel(Aaug, x, Niter = 15):
    """
    This should implement the Gauss-Seidel method (see page 860, Tabl 20.2) for solving a system of equations.
    :param Aaug: The augmented matrix from Ax=b -> [A|b]
    :param x:  An initial guess for the x vector. if A is nxn, x is nx1
    :param Niter:  Number of iterations to run the GS method
    :return: the solution vector x
    """

    Aaug = MakeDiagDom(Aaug)

    n = len(Aaug)  # Number of equations
    for _ in range(Niter):
        for i in range(n):
            # Extract the current row and the corresponding value in b
            row = Aaug[i][:-1]  # All elements except the last
            b_i = Aaug[i][-1]  # The last element in the row (b[i])

            # Sum up the contributions of the other variables
            sum_ax = sum(Aaug[i][j] * x[j] for j in range(n) if j != i)

            # Update the i-th variable using the Gauss-Seidel formula
            x[i] = (b_i - sum_ax) / row[i]

    return x

def main():
    '''
    This is a function I created for testing the numerical methods locally.
    :return: None
    '''
    #region testing GPDF
    fx = GPDF((0,0,1))
    print("{:0.5f}".format(fx))  # Does this match the expected value?
    #edregion

    #region testing Simpson
    p=Simpson(GPDF,(0,1,-5,0)) # should return 0.5
    print("p={:0.5f}".format(p))  # Does this match the expected value?
    #endregion

    #region testing Probability
    p1 = Probability(GPDF, (0,1),0,True)
    print("p1={:0.5f}".format(p1))  # Does this match the expected value?
    #endregion
    pass

#endregion

#region function calls
if __name__ == '__main__':
    main()
#endregion