# coding=utf-8
import numpy as np


class NaturalSpline(object):
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
	# x_i = [1, 2, 5, 6, 7, 8, 10, 13, 17]
	# a_i = [3, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5]
	x_i = [0.9, 1.3, 1.9, 2.1, 2.6, 3.0, 3.9, 4.4, 4.7, 5.0, 6.0, 7.0, 8.0, 9.2, 10.5, 11.3, 11.6, 12.0, 12.6, 13.0, 13.3]
	a_i = [1.3, 1.5, 1.85, 2.1, 2.6, 2.7, 2.4, 2.15, 2.05, 2.1, 2.25, 2.3, 2.25, 1.95, 1.4, 0.9, 0.7, 0.6, 0.5, 0.4, 0.25]


	clamped = NaturalSpline(x_i, a_i)
	clamped.generate()


	print clamped.result
	# print clamped.evaluate(28.1)