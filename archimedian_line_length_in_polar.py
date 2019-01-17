# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 15:44:48 2018

@author: kai.lu
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import simps
import sympy as sym
#import sys

# line integration should be between 53.6 and 63.7
Top_Curve_Ratio, Top_Curve_Ri, start_angle, end_angle = 4.464e-3, 0.5, 0, 1120

theta = np.arange(0, end_angle + 1, 1)
theta_len = len(theta)
theta_rad = theta * np.pi / 180
theta_delta = theta_rad[1] - theta_rad[0]
line_len = []

def super_archimedean(ang):
    return (Top_Curve_Ri + Top_Curve_Ratio * ang * 180 / np.pi)

def super_archimedean_from_sympy(ang):
    return (Top_Curve_Ri + Top_Curve_Ratio * ang * 180 / np.pi)

mode_num = 1
for m in range(1, mode_num + 1):
    r = (Top_Curve_Ri + Top_Curve_Ratio * theta)
    # as the abs function can not be deviated in a good way, so
    
    line_int = 0
    line_elements = []
    for i in range(theta_len - 1):
        r_delta = r[i + 1] - r[i]
        # https://www.math24.net/line-integrals-scalar-functions/
        line_element = theta_delta * (((r[i + 1] + r[i]) / 2)**2 + (r_delta / theta_delta)**2)**0.5
        line_elements.append(line_element)
        line_int += line_element
    line_int = float(format(line_int, '.1f'))
    print('The length of curve for m={} is: {}'.format(m, line_int))
    line_elements = np.asarray(line_elements)
    
    x = sym.symbols('x')
    r_symbol = super_archimedean_from_sympy(x)
    dr_symbol = sym.diff(super_archimedean_from_sympy(x), x)  # tuple(theta_rad)
    dl = (dr_symbol ** 2 + r_symbol ** 2)**0.5
    print(dl)
    
    dl_exp = sym.lambdify(x, dl, "numpy") #    dr_values = dr_symbol.evalf(subs={x: 90})
    dl_data = dl_exp(theta_rad)
    linelenth_symb = sym.integrate(dl, x)
    print(linelenth_symb)
    print("Now let's see the difference between three methods:")
    linelenth_simps = simps(dl_data, theta_rad)
    print(f'The line integral using simps method is: {linelenth_simps}')
    linelenth_quad = integrate.quad(dl_exp, theta_rad[0], theta_rad[-1])
    print(f'The line integral using quad method is: {linelenth_quad[0]}')
    print(f'The line integral using Kai method is: {line_int}')

    line_len.append(line_int)
    
    ax = plt.subplot(111, projection='polar')
    ax.plot(theta*np.pi/180, r)
    ax.set_rmax(20)
    ax.set_rticks(list(range(0, 20+1, 2)))  # Less radial ticks
    ax.set_rlabel_position(90)
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)
    
    ax.set_title("An Archimedean Curve", va='bottom')
    ax.text(-2, 9, 'Curve Length: {}'.format(line_int), fontsize=10)
    ax.text(-2, 6, r'$m={}$'.format(m), fontsize=10)
    # plt.figure(figsize=(20,10))
    plt.show()
    plt.tight_layout()
    plt.close()
print('Finished!')
