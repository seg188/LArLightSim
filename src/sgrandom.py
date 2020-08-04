import numpy as np

	
def poisson(char_length):
	def pdf(t):
		return np.exp(-1.0*t/char_length)/char_length

	return pdf

class Distribution:

	_NPOINTS = 1000

	def __init__(self, function, fmin=0.0, fmax=1.0, fHighPrecission=False):
		self._function = function
		self._min = fmin
		self._max = fmax

		self._x = [] 
		self._y = []
		integral = 0.000
		if fHighPrecission:
			self._NPOINTS = self._NPOINTS*5
		step_size = (self._max - self._min) / float(self._NPOINTS)
		for k in range(self._NPOINTS):
			self._x.append(fmin + float(k)*step_size )
			self._y.append(self._function(self._x[k]))
			integral = integral + step_size*self._y[k]

		
		_sum = 0.0
		for k in range(self._NPOINTS):
			_sum = _sum + self._y[k]/integral
			self._y[k] = _sum

		self._ymax = max(self._y)
		self._ymin = min(self._y)

	def sample(self, n=1):
		vals = []
		rtn_vals = []
		rands = np.random.rand(n)

		for k, val in enumerate(rands):
			nval = self._ymin + val*(self._ymax - self._ymin)

			for i, y in enumerate(self._y):
				if y >= nval:
					side = np.random.rand(1)
					if side >= 0.50:
						rtn_vals.append(self._x[i])
					else:
						rtn_vals.append(self._x[i-1])


					break
		return rtn_vals







