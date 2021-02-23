import math

class Bhaskara:

	def delta(self, a, b, c):
		return b**2 - 4*a*c
		
	def main(self):
		a = float(input("Digite o valor de a: "))
		b = float(input("Digite o valor de b: "))
		c = float(input("Digite o valor de c: "))
		print(self.calcula_raizes(a, b, c))
	
	def calcula_raizes(self, a, b, c):
		d = self.delta(a, b, c)
		if d == 0:
			raiz = -b/(2*a)
			return 1, raiz
		elif d < 0:
			return 0
		else:
			raiz1 = (-b + math.sqrt(d))/(2*a)
			raiz2 = (-b - math.sqrt(d))/(2*a)
			return 2, raiz1, raiz2
		
if __name__ == '__main__':

	x = Bhaskara()
	x.main()
