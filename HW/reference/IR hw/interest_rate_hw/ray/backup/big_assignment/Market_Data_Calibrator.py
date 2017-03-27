__author__ = 'Ray'

import numpy as np
import scipy.optimize as opt

class MarketRateCalibrator(object):
    """
    Class that performs optimization using scipy.optimize
    """
    def __init__(self, lambdas, tvec, rc):
        """
        Constructor of the class
        :param lambdas: The list of lambda used in Tikhonov regularizer
        :param tvec: the knots point array
        :param rc: The rate calculator
        """
        self._lambdas = lambdas

        #The knot points time array
        self._tvec = tvec

        #n is the number of knot points
        self._n = len(tvec)
        self._rc = rc
        self._x0 = np.repeat(0.01, 2 * self._n)


    def optimize(self, *fitobjs):
        """
        Perform optimization using scipy.optimize

        :param fitobjs: The array of market rate objects to be used
                        in the optimizer
        :return:
        """
        f = self.objfun(self._rc, self._n, *fitobjs)

        factr = 1e8

        bounds = None
        res = opt.lbfgsb.fmin_l_bfgs_b(f, self._x0, bounds=bounds, factr=factr, approx_grad=True)
        xopt = res[0]
        fmin = res[1]
        #print "fmin = ", fmin
        return xopt


    def square_error(self, x, y):
        """
        Function to return the square error
        :param x: input array of x
        :param y: input array of y
        :return: the square error
        """
        return sum((np.asarray(x) - np.asarray(y))**2)


    def objfun(self, rc, n, *mkt_rate_objs):
        """
        Objective function
        :param rc: the rate calculator object
        :param n: The number of knot points in the basis functions
        :param mkt_rate_objs: The list of market rate objects
        :return: the function that calculate the result of the objective function
        """
        def f(x):
            error = 0.0
            for item in mkt_rate_objs:
                error += self.square_error(item.calibrate(x[:n], x[n:]), item.act())
            pnty = rc.Tikhonov_regularizer(x[:n], x[n:], self._lambdas)
            return error + pnty
        return f
