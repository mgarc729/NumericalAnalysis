# coding=utf-8
import numpy as np


class NaturalSpline(object):
	"""Calculates the Natural Cubic Spline for a set of points


	:args
		x_i: set of points in the x coordinate
		a_i: set of f(x) points
	:param
		n: Maximum sub index of points
		x_i(np.array): set of points in the x coordinate
		a_i(np.array): set of f(x) points
		result(np.array): set of all the coefficients per spline(n-1 splines)

	"""
	def __init__(self, x_i, a_i):
		if len(x_i) is not len(a_i):
			raise ValueError("Points array(x_i) and f(x_i) array(a_i) has to have the same size")

		self.n = len(a_i) - 1
		self.a_i = np.array(a_i)
		self.x_i = np.array(x_i)
		self.result = None

	def generate(self):
		h_i = np.zeros(self.n)

		# Step 1
		"""
			For i = 0,1, ... ,n - 1 set h_i = x_i+1 - x_i
		"""
		for i in range(0, self.n):
			h_i[i] = self.x_i[i+1] - self.x_i[i]

		# Step 2
		"""
			For i = 1,2,...,n−1
			set alpha i = 3(a_i+1 −a_i)/h_i− 3(a_i −a_i−1)/h_i-1.
		"""
		alpha_i = np.zeros(self.n + 1)

		for i in range(1,self.n):
			alpha_i[i] = (3/h_i[i]) * (self.a_i[i+1] - self.a_i[i]) - (3/h_i[i-1]) * (self.a_i[i] - self.a_i[i-1])

		# Step 3
		"""
			Set l_0 = 1; 
			μ_0 = 0;
			z_0 = 0.
		"""

		l_i = np.zeros(self.n + 1)
		mu_i = np.zeros(self.n)
		z_i = np.zeros(self.n + 1)

		l_i[0] = 1
		mu_i[0] = 0
		z_i[0] = 0

		# Step 4
		"""
			For i=1,2,...,n−1
				set l_i = 2(x_i+1 − x_i−1) − h_i−1 * mu_i−1;
				mu_i = h_i/l_i;
				z_i = (alpha_i − h_i−1 * z_i−1)/l_i.
		"""

		for i in range(1, self.n):
			l_i[i] = 2 * (self.x_i[i + 1] - self.x_i[i-1]) - h_i[i-1] * mu_i[i-1]
			mu_i[i] = h_i[i] / l_i[i]
			z_i[i] = (alpha_i[i] - h_i[i-1] * z_i[i-1])/l_i[i]

		# Step 5
		"""
			Set ln = 1; 
			zn = 0; 
			cn = 0.
		"""

		c_i = np.zeros(self.n + 1)

		l_i[self.n] = 1
		z_i[self.n] = 0
		c_i[self.n] = 0

		# Step 7
		"""
			For j= n−1,n−2,...,0 
			set c_j =z_j −mu_j * c_j+1;
			b_j = (a_j+1 − a_j)/h_j − h_j(c_j+1 + 2c_j)/3; 
			d_j = (c_j+1 − c_j)/(3h_j).
		"""
		b_i = np.zeros(self.n)
		d_i = np.zeros(self.n)

		for i in range(self.n - 1, -1, -1):
			c_i[i] = z_i[i] - mu_i[i] * c_i[i+1]
			b_i[i] = (self.a_i[i+1] - self.a_i[i])/h_i[i] - h_i[i] * (c_i[i+1] + 2 * c_i[i]) / 3
			d_i[i] = (c_i[i+1] - c_i[i]) / (3 * h_i[i])

		# Step 8
		"""
			OUTPUT (aj,bj,cj,dj for j = 0,1,...,n − 1);
		"""

		result = np.zeros((self.n, 4))
		for i in range(0, self.n):
			result[i, 0] = self.a_i[i]
			result[i, 1] = b_i[i]
			result[i, 2] = c_i[i]
			result[i, 3] = d_i[i]

		self.result = result
		return result

	def _cubic(self, coefficients, x, x_i):
		"""

		:type coefficients: np.array with the coefficients of the cubic function
		"""
		return coefficients[0] + \
			coefficients[1] * (x - x_i) + \
			coefficients[2] * np.power((x - x_i), 2) + \
			coefficients[3] * np.power((x - x_i), 3)

	def evaluate(self, x):
		for i in range(0, self.n):
			if self.x_i[i] <= x <= self.x_i[i+1]:
				return self._cubic(self.result[i], x, self.x_i[i])




if __name__ == "__main__":
	x_i_1 = [1, 2, 5, 6, 7, 8, 10, 13, 17]
	a_i_1 = [3, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5]

	x_i_2 = [17, 20, 23, 24, 25, 27,27.7]
	a_i_2 = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1]

	x_i_3 = [27.7, 28, 29, 30]
	a_i_3 = [4.1, 4.3, 4.1, 3.0]

	print "Natural Spline"
	clamped_1 = NaturalSpline(x_i_1, a_i_1)
	clamped_2 = NaturalSpline(x_i_2, a_i_2)
	clamped_3 = NaturalSpline(x_i_3, a_i_3)

	print "Curve 1 a, b, c, d"
	clamped_1.generate()
	print clamped_1.result

	print "Curve 2 a, b, c, d"
	clamped_2.generate()
	print clamped_2.result

	print "Curve 3 a, b, c, d"
	clamped_3.generate()
	print clamped_3.result
