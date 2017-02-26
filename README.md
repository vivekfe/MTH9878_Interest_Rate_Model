# MTH9878_Interest_Rate_Model
This repo contains projects for Spring2017 MTH9878 Interest Rate Models class

## HW1 - Project 1
* Kernel version: Python 2.7
* Packages: pandas, numpy, datetime, matplotlib, scipy
* Data: From given xls file: DataSheetCurve.xls
* Notes:
    * The discount factor function (dist) contains the option that the inserest rate is not constant, i.e., a function of time
    * The f_const, f_instant examples are arbitrary, and can be modified to the corresponding constant or function.
* Defined funcitons to compute discount factor, LIBOR/OIS forward rates, LIBOR/OIS instantaneous rate, and present value of a swap.
* Also defined numerical method (simpson's rule) to compute the integration

