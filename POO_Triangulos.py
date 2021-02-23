'''
Define uma Classe Triangulo cujo construtor recebe 3 valores inteiros correspondentes
aos lados "a", "b", e "c" de um triângulo. A classe triângulo possui: 1) O método 
perimetro(), que não recebe parâmetros e devolve um valor inteiro correspondente ao 
perímetro do triângulo; 2) O método tipo_lado(), que devolve uma string dizendo se o 
triângulo é: isósceles, equilátero ou escaleno; 3) O método retangulo(), que devolva 
True se o triângulo for retângulo, e False caso contrário; 4) O método 
semelhantes(triangulo), que recebe um objeto do tipo Triangulo como parâmetro e 
verifica se o triângulo atual é semelhante ao triângulo passado como parâmetro. Caso 
positivo, o método devolve True. Caso negativo, devolve False.
'''

class Triangulo:

	def __init__(self, a, b, c):
		self.a = a
		self.b = b
		self.c = c
	
	def perimetro(self):
		return self.a + self.b + self.c
		
	def tipo_lado(self):
		if self.a == self.b == self.c:
			return "equilátero"
		elif self.a != self.b and self.a != self.c:
			return "escaleno"
		else:
			return "isósceles"
	
	def retangulo(self):
		l = sorted([self.a, self.b, self.c])
		if l[2]**2 == l[0]**2 + l[1]**2:
			return True
		else:
			return False
			
	def semelhantes(self, triangulo):
		l = sorted([self.a, self.b, self.c])		
		l2 = sorted([triangulo.a, triangulo.b, triangulo.c])
		if l[0]/l2[0] == l[1]/l2[1] == l[2]/l2[2]:
			return True
		else:
			return False

if __name__ == "__main__":

	t1 = Triangulo(3, 4, 5)
	t2 = Triangulo(6, 8, 10)
	print(t1.a, t1.b, t1.c)	
	print(t1.perimetro())
	print(t1.tipo_lado())
	print(t1.retangulo())
	print(t1.semelhantes(t2))
