import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

def error_abs(array1, array2):
    # returning the cumulated, over a time period normalized, mean absolute error in rad
    n_timesteps = np.size(array1)
    error_cum = 0

    for i in range(n_timesteps-1):
        error_cum = error_cum + np.abs((array1[i] - 1) - array2[i])

    error_abs = error_cum/n_timesteps
    # print('abs error: ' + str(round(np.rad2deg(error_abs(array1[:, 1], array2[:, 0].real)), 6)) + ' grad')
    return error_abs

def error_rel(array1, array2):
    # returning the cumulated, over a time period normalized, mean relative error
    n_timesteps = np.size(array1) # how many single errors; altought not working when size(array1) != size(array2)
    error_cum = 0

    # adding the singel errors up
    for i in range(n_timesteps-1):
        error_cum = error_cum + np.abs(1 - (array1[i] - 1)/array2[i])

    error_rel = error_cum/n_timesteps
    # print('rel error: ' + str(round(error_rel(array1[:, 1], array2[:, 0].real)*100, 2)) + ' %')
    return error_rel

def w_dataset():
    return

def r_datasets():
    return

def r_parameterset():
    return

def w_parameterset():
    return