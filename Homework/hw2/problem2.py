"""
Find the original signals f that correspond to them.
"""


import pywt
import math

cA_a = [2*math.sqrt(2), -math.sqrt(2), 3*math.sqrt(2), -math.sqrt(2)]
cD_a = [0, math.sqrt(2), 0, math.sqrt(2)]
y_a = pywt.idwt(cA_a, cD_a, 'haar')

cA_b = [4*math.sqrt(2), 3*math.sqrt(2), -math.sqrt(2), 2*math.sqrt(2)]
cD_b = [math.sqrt(2), -math.sqrt(2), 0, 2*math.sqrt(2)]
y_b = pywt.idwt(cA_b, cD_b, 'haar')

cA_c = [3*math.sqrt(2), 2*math.sqrt(2), 2*math.sqrt(2), 0]
cD_c = [2*math.sqrt(2), -math.sqrt(2), math.sqrt(2), 0]
y_c = pywt.idwt(cA_c, cD_c, 'haar')

cA_d = [4*math.sqrt(2), 5*math.sqrt(2), 7*math.sqrt(2), -4*math.sqrt(2)]
cD_d = [math.sqrt(2), 2*math.sqrt(2), -2*math.sqrt(2), math.sqrt(2)]
y_d = pywt.idwt(cA_d, cD_d, 'haar')

print('The original signals f that correspond to them:')
print('(a):  ' + str(y_a))
print('(b):  ' + str(y_b))
print('(c):  ' + str(y_c))
print('(d):  ' + str(y_d))
