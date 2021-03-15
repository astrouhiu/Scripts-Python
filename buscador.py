'''
Algoritmos de Busca
===================


1) Algoritmo de Busca Sequencial
   =============================
   
Compara cada elemento da lista com o elemento buscado ("x").

2) Algoritmo da Busca Binária
   ==========================

Aqui levamos em conta o fato de que a lista está ordenada. Considera-se o
elemento "m" do meio da lista:

- Se x == m, o elemento buscado foi encontrado.
- Se x < m, procura apenas na 1a metade (esquerda).
- Se x > m, procura apenas na 2a metade (direita).

O procedimento é repetido até que "x" seja encontrado ou até que a sub-lista
em questão seja vazia.
'''

class Buscador:

	def busca_sequencial(self, lista, x):
	
		for i in range(len(lista)):
			if lista[i] == x:
				return "O elemento está no índice " + str(i)
				
		return "O elemento não está na lista"
		
	def busca_binaria(self, lista, x):
	
		primeiro = 0
		ultimo = len(lista) - 1
		
		while primeiro <= ultimo:
			m = (primeiro + ultimo)//2
			if lista[m] == x:
				return "O elemento está no índice " + str(m)
			elif x < lista[m]:
				ultimo = m - 1
			else:
				primeiro = m + 1
		
		return "O elemento não está na lista"
		
if __name__ == "__main__":

	o = Buscador()
	print(o.busca_sequencial([-40, 9, 2, 45, 9, -3], 4))
	print(o.busca_sequencial(['a', 'u', 'e', 'i'], 'e'))
	print(o.busca_binaria([-3, 6, 7, 18, 25, 32, 70], 18))
	print(o.busca_binaria(['a', 'e', 'i'], 'e'))	
