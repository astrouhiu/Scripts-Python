'''
Busca uma música dentro de uma playlist
'''

class Musica:

	def __init__(self, titulo, interprete, ano):
	
		self.titulo = titulo
		self.interprete = interprete
		self.ano = ano
		
class Buscador:

	def busca_por_titulo(self, playlist, musica):
	
		for i in range(len(playlist)):
			if playlist[i].titulo.lower() == musica.lower():
				return i
		return -1
		
	def consulta(self, musica):
	
		playlist = [Musica("Hay un Universo de Pequeñas Cosas", "Alejandro Sanz", 2000),
			    Musica("Io Sì", "Laura Pausini", 2020),
			    Musica("Siempre es de Noche", "Alejandro Sanz", 1997),
			    Musica("Un Charquito de Estrellas", "Alejandro Sanz", 1997)]
		
		indice = self.busca_por_titulo(playlist, musica)
		
		if indice == -1:
			print("Essa música não está na playlist")
		else:
			resultado = playlist[indice]
			print(resultado.titulo, resultado.interprete, resultado.ano, sep = ', ')

if __name__ == "__main__":
			
	b = Buscador()
	b.consulta("Siempre es de noche")	
