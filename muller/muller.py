import math

## Define muller Potential
def f(x,y):
	A_vec = [-200.,-100.,-170.,15]
	a_vec = [-1.,-1.,-6.5,0.7]
	b_vec = [0.,0.,11.,0.6]
	c_vec = [-10.,-10.,-6.5,0.7]
	x_vec = [1.,0.,-0.5,-1.]
	y_vec = [0.,0.5,1.5,1.]

	n = len(A_vec)

	V = 0.
	x = 0.
	y = 1.5
	for i in range (n):
		V += A_vec[i]*math.exp(a_vec[i]*(x-x_vec[i])**2 + b_vec[i]*(x-x_vec[i])*(y-y_vec[i]) + c_vec[i]*(y-y_vec[i])**2)

	return V


def fx(x,y):
	A_vec = [-200.,-100.,-170.,15]
	a_vec = [-1.,-1.,-6.5,0.7]
	b_vec = [0.,0.,11.,0.6]
	c_vec = [-10.,-10.,-6.5,0.7]
	x_vec = [1.,0.,-0.5,-1.]
	y_vec = [0.,0.5,1.5,1.]

	n = len(A_vec)

	V = 0.
	for i in range (n):
		V += A_vec[i]*a_vec[i]*2*(x-x_vec[i])*math.exp(a_vec[i]*(x-x_vec[i])**2 + b_vec[i]*(x-x_vec[i])*(y-y_vec[i]) + c_vec[i]*(y-y_vec[i])**2)
		V += A_vec[i]*b_vec[i]*(y-y_vec[i])*math.exp(a_vec[i]*(x-x_vec[i])**2 + b_vec[i]*(x-x_vec[i])*(y-y_vec[i]) + c_vec[i]*(y-y_vec[i])**2)

	return V

def fy(x,y):
	A_vec = [-200.,-100.,-170.,15]
	a_vec = [-1.,-1.,-6.5,0.7]
	b_vec = [0.,0.,11.,0.6]
	c_vec = [-10.,-10.,-6.5,0.7]
	x_vec = [1.,0.,-0.5,-1.]
	y_vec = [0.,0.5,1.5,1.]

	n = len(A_vec)

	V = 0.
	for i in range (n):
		V += A_vec[i]*c_vec[i]*2*(y-y_vec[i])*math.exp(a_vec[i]*(x-x_vec[i])**2 + b_vec[i]*(x-x_vec[i])*(y-y_vec[i]) + c_vec[i]*(y-y_vec[i])**2)
		V += A_vec[i]*b_vec[i]*(x-x_vec[i])*math.exp(a_vec[i]*(x-x_vec[i])**2 + b_vec[i]*(x-x_vec[i])*(y-y_vec[i]) + c_vec[i]*(y-y_vec[i])**2)

	return V