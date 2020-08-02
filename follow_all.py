'''
Este codigo produz o arquivo follow.out com freq. 1, usando o executavel follow_swift-all.x em cada uma das subpastas de uma dada corrida e muda o nome segundo a enumeracao da subpasta (follow.out -> follow_simxyz.out). Todos esses arquivos sao armazenados numa nova pasta para poder fazer a descarga dela posteriormente.
'''

import os

# ********************************** Dados de Entrada **********************************

# nproc:           num. de processadores usados ou num. de subpastas de uma dada corrida
# pasta_corrida:   nome da pasta da corrida que contem as subpastas
# pasta_follow:    nome da pasta que sera criada e contera os arquivos follow_simxyz

nproc = 4
pasta_corrida = '/home/jessica/caceres/Dispersos/q9_58_au_parte_2/'
pasta_follow = '/home/jessica/caceres/Dispersos/Follow_q9_58_au_parte_2/' 
rota_executavel = '/home/jessica/programs/swift/tools/follow_swift-all.x'

# **************************************************************************************

comando1 = 'mkdir ' + pasta_follow
os.system(comando1)

for i in range(1, nproc + 1):

        num = 'sim%3.3d'%i
	path = pasta_corrida + num 
	os.chdir(path)

	comando2 = rota_executavel + '<../../../frecuencia.dat'
	os.system(comando2)
	
	comando3 = 'mv follow.out ' + pasta_follow + 'follow_' + num + '.out'
	os.system(comando3)
