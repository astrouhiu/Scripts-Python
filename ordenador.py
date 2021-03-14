'''
Alguns Algoritmos de Ordenação
==============================


1) Algoritmo da Seleção Direta ou Selection Sort
   =============================================

Ordena a lista selecionando a cada i-ésima iteração o menor elemento no pedaço
ainda não ordenado da lista e o coloca na i-ésima posição da lista.

2) Algoritmo de Ordenação da Bolha ou Bubble Sort
   ==============================================

Pense na lista como um tubo de ensaio vertical onde os elementos mais leves
sobem à superfície como uma bolha e os mais pesados afundam.

O algoritmo consiste em percorrer a lista múltiplas vezes; a cada passagem, 
compara todos os elementos adjacentes e troca de lugar os que estiverem fora de
ordem. A cada passagem, o elemento mais pesado se deposita no fundo do tubo de 
ensaio, i.e., o maior elemento da lista vai para o final dela.

3) Melhoria no Algoritmo da Bolha
   ==============================

A melhoria no Algoritmo de Ordenação da Bolha consiste em que, se em uma das
iterações nenhuma troca é realizada, o que significa que a lista já está 
ordenada, finalizamos o algoritmo.

4) Algoritmo de Ordenação por Inserção ou Insertion Sort
   =====================================================

A cada iteração "i" do algoritmo, troca de posição o elemento lista[i] com cada 
elemento maior que ele que esteja à sua esquerda.

5) Algoritmo de Ordenação por Intercalação ou Merge Sort
   =====================================================

Divide a lista na metade recursivamente, até que cada sublista contém apenas
1 elemento. Repetidamente, intercala as sublistas para produzir novas listas 
ordenadas. Isto é repetido até que tenhamos apenas 1 lista no final (que 
estará ordenada).
'''	

class Ordenador:

	def selection_sort(self, lista):
	
		fim = len(lista)
	
		for i in range(fim - 1):
			pos_menor = i
			for j in range(i + 1, fim):
				if lista[j] < lista[pos_menor]:
					pos_menor = j
			lista[i], lista[pos_menor] = lista[pos_menor], lista[i]

		return lista

	def bubble_sort(self, lista):
	
		fim = len(lista)
		
		for i in range(fim - 1, 0, -1):
			for j in range(i):
				if lista[j] > lista[j + 1]:
					lista[j], lista[j + 1] = lista[j + 1], lista[j]
					
		return lista
		
	def short_bubble_sort(self, lista):
	
		fim = len(lista)
		
		if fim <= 1:
			return lista
		
		for i in range(fim - 1, 0, -1):
			trocou = False
			for j in range(i):
				if lista[j] > lista[j + 1]:
					lista[j], lista[j + 1] = lista[j + 1], lista[j]
					trocou = True
			if not trocou:
				return lista
				
	def insertion_sort(self, lista):

		for i in range(1, len(lista)):
			x = lista[i]
			j = i - 1
			while j >= 0 and x < lista[j]:
				lista[j + 1] = lista[j]
				j -= 1
			lista[j + 1] = x
		
		return lista
		
	def merge_sort(self, lista):

		if len(lista) <= 1:
			return lista
	
		meio = len(lista)//2
	
		lado_esquerdo = self.merge_sort(lista[:meio])
		lado_direito = self.merge_sort(lista[meio:])

		return self.merge(lado_esquerdo, lado_direito)

	def merge(self, lado_esquerdo, lado_direito):

		if not lado_esquerdo:
			return lado_direito	
		if not lado_direito:
			return lado_esquerdo

		if lado_esquerdo[0] < lado_direito[0]:
			return [lado_esquerdo[0]] + self.merge(lado_esquerdo[1:], lado_direito)
		return [lado_direito[0]] + self.merge(lado_esquerdo, lado_direito[1:])
	
if __name__ == "__main__":

	lista_1 = [0, -7, 9, 5, 25, 30]
	lista_2 = lista_1[:]
	lista_3 = lista_1[:]
	lista_4 = lista_1[:]
	lista_5 = lista_1[:]
	
	o = Ordenador()
	print(o.selection_sort(lista_1))
	print(o.bubble_sort(lista_2))
	print(o.short_bubble_sort(lista_3))
	print(o.insertion_sort(lista_4))
	print(o.merge_sort(lista_5))
