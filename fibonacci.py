# Fibonacci Numbers Module. Extraído da documentação oficial de Python 3

# Write Fibonacci series up to n
def fib(n): 
	a, b = 0, 1
	while b < n:
		print(b, end = " ")
		a, b = b, a + b
	print()
	
# Return Fibonacci series up to n
def fib2(n):	
	result = []
	a, b = 0, 1
	while b < n:
		result.append(b)
		a, b = b, a + b
	return result
	
print(__name__)

# A variável __name__ guarda o nome do módulo. A mudança no nome pode nos 
# ajudar a verificar se estamos executando como script ou importando como 
# módulo dentro de outro código. Se true, significa que estou executando 
# ele como script. Se false, está sendo importado como módulo em outro 
# lugar, para simplesmente usar as funções

if __name__ == "__main__":
	# Aqui vai o que queremos que seja executado quando for rodado 
	# como script
	
	# Módulo sys: Permite interagir com o sistema operacional
	import sys
	fib(int(sys.argv[1]))
	
	# Para executar: python fibonacci.py sys.argv[1], onde sys.argv[1]
	# é o parâmetro da função
