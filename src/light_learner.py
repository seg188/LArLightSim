import numpy as np 
from sklearn.neural_network import MLPRegressor
class light_learner:

	def __init__(self, layer_sizes=(2000, 2000), activation_function='logistic', _max_iter=277):
		self._regressor = MLPRegressor(hidden_layer_sizes=layer_sizes, activation=activation_function,  max_iter=_max_iter, alpha=0.01)

	def train(self, x, y):
		self._regressor.fit(x, y)

class light_classifier:
	def __init__(self, layer_sizes=(500, 250, 100), activation_function='relu', _max_iter=277):
		self._class = MLPClassifier(hidden_layer_sizes=layer_sizes, activation=activation_function,  max_iter=_max_iter, alpha=0.1)

	def train(self, x, y):
		self._regressor.fit(x, y)
