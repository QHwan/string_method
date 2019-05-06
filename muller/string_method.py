import numpy as np
import sys
import random
import math

from muller import (f, fx, fy)
from lib import cubic_spline, runge_kutta



y_vec = np.arange(-0.2,1.5+0.1,0.1)
n = len(y_vec)
x_vec = np.zeros(n)
xi_vec = x_vec
yi_vec = y_vec

step = 1.
tol = 10.
tol_vec = np.zeros(n)
s_vec = np.zeros(n) # not n-1, s_vec[0] = 0.
step = 0

while (tol > max(10.**-4,1e-10)):
	step += 1

	dt = 0.001 * min(0.2,1./n)

	xn_vec = np.zeros(n)
	yn_vec = np.zeros(n)

	# Runge-Kutta method
	# Solve image_i+1 = image_i - dt*Del(V(image_i))
	xn_vec, yn_vec = runge_kutta(fx, fy, x_vec, y_vec, dt)
	

	# Parameterization by equal arc length
	# s is the arc length
	# s0 = 0, s_last = 1
	# s_i = s_i-1 + |image_i - image_i-1|
	s_vec = np.zeros(n)
	for i in range (n):
		if i==0:
			s_vec[i] = 0.
		else:
			s_vec[i] = s_vec[i-1] + ((xn_vec[i]-xn_vec[i-1])**2 + (yn_vec[i]-yn_vec[i-1])**2)**0.5

	# normalize s_vec
	for i in range (n):
		s_vec[i] /= s_vec[n-1]





	# Arc length is not equidistant.
	# So we relocate image by doing cubic spline along s
	# x, y = a + bs + cs2 + ds3
	ax_vec, bx_vec, cx_vec, dx_vec = cubic_spline(s_vec,xn_vec)
	ay_vec, by_vec, cy_vec, dy_vec = cubic_spline(s_vec,yn_vec)


	# re-locate the s_vec & x_vec, y_vec, with uniform distance

	dist = 1./(n-1)

	snew_vec = np.zeros(n)
	for i in range (n):
		snew_vec[i] = dist*i

	xnew_vec = np.zeros(n)
	ynew_vec = np.zeros(n)
	for i in range (n):
		if i==0:
			xnew_vec[0] = xn_vec[0]
			ynew_vec[0] = yn_vec[0]
		elif i==n-1:
			xnew_vec[n-1] = xn_vec[n-1]
			ynew_vec[n-1] = yn_vec[n-1]  # first and last points are always same
		else:
			s = snew_vec[i]
			for j in range (n-1):
				if s_vec[j]<s and s_vec[j+1]>s:
					ax = ax_vec[j+1]
					bx = bx_vec[j+1]
					cx = cx_vec[j+1]
					dx = dx_vec[j+1]
					ay = ay_vec[j+1]
					by = by_vec[j+1]
					cy = cy_vec[j+1]
					dy = dy_vec[j+1]
					break
			xnew = ax + bx*s + cx*s*s + dx*s*s*s
			ynew = ay + by*s + cy*s*s + dy*s*s*s

			xnew_vec[i] = xnew
			ynew_vec[i] = ynew


	# calculate the tol
	tol_vec = np.zeros(n)
	for i in range (n):
		if i==0 or i==n-1:
			tol_vec[i] = 0.
		else:
			tol_vec[i] = ((xnew_vec[i]-x_vec[i])**2 + (ynew_vec[i]-y_vec[i])**2)**0.5

	tol = abs(np.amax(tol_vec))

	dist_vec = np.zeros(n-1)
	for i in range (len(dist_vec)):
		dist_vec[i] = ((xnew_vec[i+1]-xnew_vec[i])**2 + (ynew_vec[i+1]-ynew_vec[i])**2)**0.5


	print(step, tol)

	# update density profile
	x_vec = xnew_vec
	y_vec = ynew_vec

	if step == 10:
		oMat = []
		for i in range (n):
			oMat.append([xnew_vec[i],ynew_vec[i]])
		np.savetxt(str(step)+'.xvg',oMat,fmt='%e')

	


					
oMat = []
for i in range (n):
	oMat.append([xi_vec[i],yi_vec[i],xnew_vec[i],ynew_vec[i]])

np.savetxt('1.xvg',oMat,fmt='%e')
		


