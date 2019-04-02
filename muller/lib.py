import numpy as np

def div_dif(xb,x,xn,yb,y,yn):
	return ((yn-y)/(xn-x) - (y-yb)/(x-xb))/(xn-xb)

def cubic_spline(x_vec,y_vec):
	n = len(x_vec)
	A = np.zeros((n,n))
	x = np.zeros(n)
	b = np.zeros(n)

	for i in range (n):
		if i==0:
			lam = 0.
			A[i][i] = 2.
			A[i][i+1] = lam
		elif i==n-1:
			mu = 0.
			A[i][i] = 2.
			A[i][i-1] = mu
		else:
			x = x_vec[i]
			xn = x_vec[i+1]
			xb = x_vec[i-1]
			h = x-xb
			hn = xn-x
			mu = h/(h+hn)
			lam = 1.-mu
			A[i][i] = 2.
			A[i][i-1] = mu
			A[i][i+1] = lam

	for i in range (n):
		if i==0:
			x = x_vec[i]
			xn = x_vec[i+1]
			y = y_vec[i]
			yn = y_vec[i+1]
			d = 0.
		elif i==n-1:
			x = x_vec[i]
			xb = x_vec[i-1]
			y = y_vec[i]
			yb = y_vec[i-1]
			d = 0.
		else:
			x = x_vec[i]
			xn = x_vec[i+1]
			xb = x_vec[i-1]
			y = y_vec[i]
			yn = y_vec[i+1]
			yb = y_vec[i-1]
			d = 6*div_dif(xb,x,xn,yb,y,yn)
		b[i] = d


	Ainv = np.linalg.inv(A)
	sol_vec = Ainv.dot(b)

	a = np.zeros(n)
	b = np.zeros(n)
	c = np.zeros(n)
	d = np.zeros(n)

	for i in range (n):
		if i==0:
			a[i] = 0.
			b[i] = 0.
			c[i] = 0.
			d[i] = 0.
		else:
			x = x_vec[i]
			xb = x_vec[i-1]
			y = y_vec[i]
			yb = y_vec[i-1]
			M = sol_vec[i]
			Mb = sol_vec[i-1]
			h = x-xb

			a[i] = (Mb*x*x*x/6./h) - (M*xb*xb*xb/6./h) + (yb-Mb*h*h/6.)*x/h - (y-M*h*h/6.)*xb/h
			b[i] = -(Mb*x*x/2./h) + (M*xb*xb/2./h) - (yb-Mb*h*h/6.)*1/h + (y-M*h*h/6.)*1/h
			c[i] = (Mb*x/2./h) - (M*xb/2./h)
			d[i] = (1./6./h)*(-Mb+M)

	return a, b, c ,d

def runge_kutta(fx, fy, x_vec, y_vec, dt):
	n = len(x_vec)
	xn_vec = np.zeros(n)
	yn_vec = np.zeros(n)

	for i in range (n):
		xi = x_vec[i]
		yi = y_vec[i]

		kx1 = dt*fx(xi,yi)
		ky1 = dt*fy(xi,yi)
		kx2 = dt*fx(xi+0.5*kx1,yi+0.5*ky1)
		ky2 = dt*fy(xi+0.5*kx1,yi+0.5*ky1)
		kx3 = dt*fx(xi+0.5*kx2,yi+0.5*ky2)
		ky3 = dt*fy(xi+0.5*kx2,yi+0.5*ky2)
		kx4 = dt*fx(xi+kx3,yi+ky3)
		ky4 = dt*fy(xi+kx3,yi+ky3)

		xin = xi - (1./6.)*kx1 - (1./3.)*kx2 - (1./3.)*kx3 - (1./6.)*kx4
		yin = yi - (1./6.)*ky1 - (1./3.)*ky2 - (1./3.)*ky3 - (1./6.)*ky4

		xn_vec[i] = xin
		yn_vec[i] = yin

	return xn_vec, yn_vec