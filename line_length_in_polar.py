# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 10:23:37 2018

@author: kai.lu
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import simps
import sympy as sym
import sys


    
a, b, m, n1, n2, n3, Top_Curve_Ratio, Top_Curve_Ri, start_angle, end_angle = 1, 1, 10, 5, 5, 5, 0.069, 4.52, 0, 150

theta = np.arange(0, 151, 1)
theta_len = len(theta)
theta_rad = theta * np.pi / 180
theta_delta = theta_rad[1] - theta_rad[0]
line_len = []


def super_archimedean(ang):
    return (abs(np.cos(m * ang / 4) / a) ** n2 + abs(np.sin(m * ang / 4) / b) ** n3)**(-1 / n1)\
        * (Top_Curve_Ri + Top_Curve_Ratio * ang * 180 / np.pi)

def super_archimedean_from_sympy(ang):
    return ((sym.cos(m * ang / 4) / a) ** n2 + (sym.sin(m * ang / 4) / b) ** n3)**(-1 / n1)\
        * (Top_Curve_Ri + Top_Curve_Ratio * ang * 180 / np.pi)

mode_num = 10
for m in range(1, mode_num + 1):
#    print(f"value of m is: {m}")
    r = (abs(np.cos(m * theta_rad / 4) / a) ** n2 + abs(np.sin(m * theta_rad / 4) / b) ** n3)**(-1 / n1)\
        * (Top_Curve_Ri + Top_Curve_Ratio*theta)
    # as the abs function can not be deviated in a good way, so
    # to verify the line integration
#    r = 1 * theta_rad / theta_rad
#     plt.close()
    
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
    
#    x = sym.symbols('x')
#    dr_from_symbol = sym.diff(super_archimedean_from_sympy(x), x)  # tuple(theta_rad)
#    print(dr_from_symbol)
#
#    dr_values = dr_from_symbol.evalf(subs={x: 90})
#    dr_exp = sym.lambdify(x, dr_from_symbol, "numpy")
#    f = dr_exp(theta)
#    print(f)
#    dr_values = dr_exp(theta_rad)
#     dl_from_symbol = theta_delta * ((r)**2 + (dr_values)**2)**0.5
#    print(dl_from_symbol - line_elements)
#     sys.exit()

#    line_int_by_sy_integrate = float(format(integrate.quad(lambda x: super_archimedean(x), theta_rad[0], theta_rad[-1])[0], '.1f'))
#    print('The length of curve for m={} by_sym_integrate is: {}'.format(m, line_int_by_sy_integrate))
#    line_int_by_sy_simps = float(format(simps(r, theta_rad), '.1f'))
#    print('The length of curve for m={} by_sym_simps is: {}'.format(m, line_int_by_sy_simps))
#    print('The integral with the scipy methods need the infinitesimal to the curve lengths')

    line_len.append(line_int)
    
    ax = plt.subplot(111, projection='polar')
    ax.plot(theta*np.pi/180, r)
    ax.set_rmax(20)
    ax.set_rticks(list(range(0, 20+1, 2)))  # Less radial ticks
    ax.set_rlabel_position(90)
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)
    
    ax.set_title("A Super-Archimedean Curve", va='bottom')
    ax.text(-2, 9, 'Curve Length: {}'.format(line_int), fontsize=10)
    ax.text(-2, 6, r'$m={}$'.format(m), fontsize=10)
    # plt.figure(figsize=(20,10))
    plt.show()
    plt.close()
print('Finished!')
