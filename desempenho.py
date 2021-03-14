import ordenador
import random
import time

class Desempenho:

	def lista_aleatoria(self, n):
	
		# Geramos uma lista de "n" números inteiros aleatórios entre 0 e 999
		lista = [random.randrange(1000) for i in range(n)]
			
		return lista
	
	def lista_quase_ordenada(self, n):
	
		# Lista ordenada
		lista = [i for i in range(n)]
		# Colocamos um número fora de ordem no primeiro décimo
		lista[n//10] = -500
		
		return lista
			
	def comparacao(self, n):
	
		lista_1 = self.lista_aleatoria(n)
		# Clonamos a lista_1 para que os outros métodos de ordenação não recebam
		# a lista_1 já ordenada, após executar o primeiro método 
		lista_2 = lista_1[:]
		lista_3 = lista_1[:]
		lista_4 = lista_1[:]
		lista_5 = lista_1[:]

		o = ordenador.Ordenador()
		
		print("Comparando os desempenhos dos algoritmos com uma lista aleatória")
		inicio = time.time()
		o.selection_sort(lista_1)
		fim = time.time()
		print("A execução do Algoritmo Selection Sort demorou:", fim - inicio, "segundos")	
		inicio = time.time()
		o.bubble_sort(lista_2)
		fim = time.time()
		print("A execução do Algoritmo Bubble Sort demorou:", fim - inicio, "segundos")
		inicio = time.time()
		o.short_bubble_sort(lista_3)
		fim = time.time()
		print("A execução do Algoritmo Bubble Sort curto demorou:", fim - inicio, "segundos")
		inicio = time.time()
		o.insertion_sort(lista_4)
		fim = time.time()
		print("A execução do Algoritmo Insertion Sort demorou:", fim - inicio, "segundos")
		inicio = time.time()
		o.merge_sort(lista_5)
		fim = time.time()
		print("A execução do Algoritmo Merge Sort demorou:", fim - inicio, "segundos")
		
		lista_1 = self.lista_quase_ordenada(n)
		lista_2 = lista_1[:]
		lista_3 = lista_1[:]
		lista_4 = lista_1[:]
		lista_5 = lista_1[:]
		
		print("\nComparando os desempenhos dos algoritmos com uma lista quase ordenada")
		inicio = time.time()
		o.selection_sort(lista_1)
		fim = time.time()
		print("A execução do Algoritmo Selection Sort demorou:", fim - inicio, "segundos")		
		inicio = time.time()
		o.bubble_sort(lista_2)
		fim = time.time()
		print("A execução do Algoritmo Bubble Sort demorou:", fim - inicio, "segundos")		
		inicio = time.time()
		o.short_bubble_sort(lista_3)
		fim = time.time()
		print("A execução do Algoritmo Bubble Sort curto demorou:", fim - inicio, "segundos")
		inicio = time.time()
		o.insertion_sort(lista_4)
		fim = time.time()
		print("A execução do Algoritmo Insertion Sort demorou:", fim - inicio, "segundos")	
		inicio = time.time()
		o.merge_sort(lista_5)
		fim = time.time()
		print("A execução do Algoritmo Merge Sort demorou:", fim - inicio, "segundos")	
	
if __name__ == "__main__":
		
	o = Desempenho()
	o.comparacao(900)		
