# coding=utf-8
import numpy as np


class ClampedSpline(object):
	def __init__(self, x_i, a_i, fp0, fpn):
		if len(x_i) is not len(a_i):
			raise ValueError("Points array(x_i) and f(x_i) array(a_i) has to have the same size")

		self.n = len(a_i) - 1
		self.a_i = np.array(a_i)
		self.x_i = np.array(x_i)
		self.fp0 = fp0
		self.fpn = fpn
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
			Set alpha_0 = 3(a_1 - a_0)/h_0 - 3FP0
				alpha_n = 3FPN - 3(a_n - a_n-1)/h_n-1
		"""
		alpha_i = np.zeros(self.n + 1)

		alpha_i[0] = (3 * (self.a_i[1] - self.a_i[0])/h_i[0]) - 3 * self.fp0
		alpha_i[self.n] = 3 * self.fpn - (3 * (self.a_i[self.n] - self.a_i[self.n - 1]))/h_i[self.n - 1]

		# Step 3
		"""
			For =1,2,...,n−1
				set alpha_i = 3(alpha_i+1 − alpha_i)/h_i− 3 (alpha_i − alpha_i−1)/h_i−1
		"""

		for i in range(1, self.n):
			alpha_i[i] = (3/h_i[i]) * (alpha_i[i+1] - alpha_i[i]) - (3/h_i[i-1]) * (alpha_i[i] - alpha_i[i-1])

		# Step 4
		"""
		Set l_0 = 2h_0; 
			mu_0 = 0.5;
			z_0 = alpha_0/l_0.
		"""
		l_i = np.zeros(self.n + 1)
		mu_i = np.zeros(self.n)
		z_i = np.zeros(self.n + 1)

		l_i[0] = 2 * (h_i[0])
		mu_i[0] = 0.5
		z_i[0] = alpha_i[0] / l_i[0]

		# Step 5
		"""
			For i = 1,2,...,n−1
				set l_i = 2(x_i+1 − x_i−1) − h_i−1* mu_i−1
				mu_i = h_i/l_i;
				z_i = (alpha_i − h_i−1 * z_i−1)/l_i.
		"""

		for i in range(1, self.n):
			l_i[i] = 2 * (self.x_i[i + 1] - self.x_i[i-1]) - h_i[i-1] * mu_i[i-1]
			mu_i[i] = h_i[i] / l_i[i]
			z_i[i] = (alpha_i[i] - h_i[i-1] * z_i[i-1])/l_i[i]

		# Step 6
		"""
			Set l_n = h_n−1(2 − mu_n−1);
			z_n = (alpha_n − h_n−1 * z_n−1)/l_n;
			c_n = z_n.
		"""

		c_i = np.zeros(self.n + 1)

		l_i[self.n] = h_i[self.n - 1] * (2 - mu_i[self.n - 1])
		z_i[self.n] = (alpha_i[self.n] - h_i[self.n-1] * z_i[self.n-1]) / l_i[self.n]
		c_i[self.n] = z_i[self.n]

		# Step 7
		"""
			For j = n−1,n−2,...,0 
				set c_j =z_j − mu_j *c_j+1;
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
			print "i",i
			if self.x_i[i] <= x <= self.x_i[i+1]:
				print i
				return self._cubic(self.result[i], x, self.x_i[i])


if __name__ == "__main__":
	x_i_1 = [1, 2, 5, 6, 7, 8, 10, 13, 17]
	a_i_1 = [3, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5]

	x_i_2 = [17, 20, 23, 24, 25, 27, 27.7]
	a_i_2 = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1]

	x_i_3 = [27.7, 28, 29, 30]
	a_i_3 = [4.1, 4.3, 4.1, 3.0]

	fp0_1 = 1.0
	fp0_2 = 3.0
	fp0_3 = 0.33

	fpn_1 = -0.67
	fpn_2 = -4.0
	fpn_3 = -1.5

	print "Clamped Spline"
	clamped_1 = ClampedSpline(x_i_1, a_i_1, fp0_1, fpn_1)
	clamped_2 = ClampedSpline(x_i_2, a_i_2, fp0_2, fpn_2)
	clamped_3 = ClampedSpline(x_i_3, a_i_3, fp0_3, fpn_3)

	print "Curve 1 a, b, c, d"
	clamped_1.generate()
	print clamped_1.result

	print "Curve 2 a, b, c, d"
	clamped_2.generate()
	print clamped_2.result

	print "Curve 3 a, b, c, d"
	clamped_3.generate()
	print clamped_3.result
