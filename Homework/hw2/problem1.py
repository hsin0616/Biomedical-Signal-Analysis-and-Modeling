"""
Find the signal f = (2, 2, 2, 4, 4, 4), find its first level Haar transform.
"""


import pywt


x = [2, 2, 2, 4, 4, 4]
cA, cD = pywt.dwt(x, 'haar')

print('Approximation: ')
print(cA)
print('Detail: ')
print(cD)
print('The first level haar transform of the signal f: ' + str(cA) + ' | ' + str(cD))