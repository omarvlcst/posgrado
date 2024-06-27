import numpy as np
import matplotlib.pyplot as plt

class EstrellaUltracompacta:
	def __init__(self,masa,radio,distancia,radio_maximo,alpha,beta,gamma,temp_par,temp_perp)
		self.masa = masa
		self.radio = radio
		self.distancia = distancia
		self.radio_maximo = radio_maximo
		self.alpha = alpha
		self.beta = beta
		self.gamma = gamma
		self.temp_par = temp_par
		self.temp_perp = temp_perp

def chi_0(self):
	print(f"El cociente chi0 es: {self.temp_perp/self.temp_par}")
	
def compacidad_estrella(self):
	print(f"La compacidad de la estrella es: {self.radio/self.masa}")
	
def radio_maximo(self):
	if self.radio < 1.5:
		return 1.5*2*self.masa
	else
		return self.radio/np.sqrt(1-2*self.masa/self.radio)
	
estrella = EstrellaUltracompacta(2.99,10,2e16,18,30,60,25,4e5,3e5)
