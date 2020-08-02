'''
Este codigo produz o arquivo follow.out com freq. 1, usando o executavel follow_symba.x (porque o arquivo tp.in nao existe neste tipo de corridas) em cada uma das subpastas de uma dada corrida e muda o nome segundo a enumeracao do planeta e da subpasta (follow.out -> followi_simxyz.out). Todos esses arquivos sao armazenados numa nova pasta para poder fazer a descarga dela posteriormente.
'''

import os

# ************************************************* Dados de Entrada ****************************************************

# nproc:           num. de processadores usados ou num. de subpastas de uma dada corrida
# pasta_corrida:   rota da pasta da corrida que contem as subpastas
# pasta_follow:    rota da pasta que sera criada e contera os arquivos followi_simxyz.out, onde "i" e o indice do planeta

nproc = 361
pasta_corrida = '/home/jessica/caceres/Classicos/lamb_q9_62_au/'
pasta_follow = '/home/jessica/caceres/Classicos/Follow_lamb_q9_62_au/' 
rota_executavel = '/home/jessica/programs/swift/tools/follow_symba.x'

# ***********************************************************************************************************************

comando1 = 'mkdir ' + pasta_follow
os.system(comando1)

for i in range(1, nproc + 1):

        num = 'sim%3.3d'%i
	path = pasta_corrida + num 
	os.chdir(path)

        for j in range(2, 7):

		comando2 = rota_executavel + '<../../../dados' + str(j) + '.dat'
		os.system(comando2)
	
		comando3 = 'mv follow.out ' + pasta_follow + 'follow' + str(j) + '_' + num + '.out'
		os.system(comando3)
