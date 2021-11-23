# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:56:50 2021

@author: amaro
"""

#REFERENCES: Code inspired and adapted from Core Worksheets 1/2/3/4

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from scipy.optimize import curve_fit

dis_data = sp.loadtxt("Distance_Mpc.txt")
s_dis_data = dis_data[sp.argsort(dis_data[:, 0])] #sorts the data numerically by the 1st column
c_s_dis_data = s_dis_data[sp.where(s_dis_data[:,2]==1)] #produces an array where all values in the 3rd column are 1 (clean data)

spec_file = sp.loadtxt("Halpha_spectral_data.csv", skiprows=4, delimiter=',')
freq = spec_file[0,1:] #frequency data
spec_data = spec_file[1:,:] #spectral data
s_spec_data = spec_data[sp.argsort(spec_data[:,0])] #sorts the data numerically by the 1st column
c_s_spec_data = s_spec_data[sp.where(s_dis_data[:,2]==1)] #removes rows from the array corresponding to the "unclean data" in the other file

l_speed = 2.9979*10**8
le = 656.28e-9 #emitted wavelength
vel = [] #empty array that will be used to feed velocity data into
distance = c_s_dis_data[:,1] #distance data

def initial_guesses(x,y): #function used to calculate all the initial guesses for the parameters in the Gaussian equation
    a = (1/3)*(max(y)-min(y)) #amplitude of Gaussian is approximately one 3rd of the highest y-value minus the lowest y-value
    m,c = np.polyfit(x,y,1) #fits a 1st order polynomial to the data and extracts the gradient and y-intercept
    residual = 0
    for i in range(1000): #iterates through all 1000 x-values
        if y[i] - (m*x[i]+c) > residual: #determines the value of x furthest away from the fitted straight line (peak of Gaussian = mu)
            residual = y[i] - (m*x[i]+c) #updates residual variable if larger residual found
            mu = x[i] #updates mu variable if larger residual found
    return [a,mu,1e13,m,c] #no code to guess sig as approximately same for all curves (of order e13)

def fit_func(x,a,mu,sig,m,c): #given initial parameters, draws the sum of a Gaussian curve and straight line
    gaus = a*sp.exp(-(x-mu)**2/(2*sig**2))
    line = m*x+c
    return gaus + line

for j in range(0,25): #25 good values (5 values removed when cleaning data) - iterates the code through all
    ig = initial_guesses(freq,c_s_spec_data[j,1:]) #determines the initial guesses for a specific iteration of the for loop
    fit = curve_fit(fit_func, freq, c_s_spec_data[j,1:],ig) #determines the optimized parameters for a curve fitted
    data_fit = fit_func(freq, *fit[0]) #recreates the fitted curve using the optimized parameters

    po,po_cov = fit #stores the optimized parameters and covariance matrix in separate variables
    def neg_fit_func(x): #negative version of Gaussian curve so we can find the minimum value using fmin which corresponds to maximum
        gaus = po[0]*sp.exp(-(x-po[1])**2/(2*po[2]**2))
        line = po[3]*x+po[4]
        return -(gaus + line)
    max_f = float(op.fmin(neg_fit_func,po[1])) #stores the maxima of Gaussian curve as a float value
    w_length = l_speed/max_f #calculates wavelength
    v = l_speed * (w_length**2 - le**2)/(le**2 + w_length**2) #calculates red-shift velocity
    vel.append(v) #adds this to vel array that was defined earlier
    
    plt.plot(freq, c_s_spec_data[j,1:]) #plots frequency against spectra
    plt.plot(freq, data_fit) #plots the fitted curve
    plt.show()

velocity = np.array(vel) #creates an array based on the list of velocties (vel)
clean_distance = np.delete(distance, [10,18,19]) #removes data points which gave odd results (negative values for velocity)
clean_velocity = np.delete(velocity, [10,18,19]) #removes data points which gave odd results (negative values for velocity)
plt.plot(clean_distance, clean_velocity, 'o', mew=1, ms=5, color='red') #plots clean velocity data against clean distance data

fits,cov = sp.polyfit(clean_distance,clean_velocity,1,cov=True) #fits a 1st order polynomial to the data and extracts the gradient and y-intercept
fit_values = sp.poly1d(fits)
plt.plot(clean_distance,fit_values(clean_distance)) #plots the line of best fit
plt.xlabel("Distance (Mpc)") #plots x-label
plt.ylabel("Redshift (km/s)") #plots y-label
plt.savefig("Redshift_V-Distance.png", dpi = 500)
plt.show()

gradient = fits[0]/1000 #extracts gradient from fits array
uncertainty = np.sqrt(cov[0,0])/1000 #extracts uncertainty in gradient from cov array
print("The calculated value for H0 is:", "%.3f +/- %.3f" %(gradient, uncertainty), "km/s/Mpc") #prints Hubble's Constant = gradient/1000 and associated uncertainty