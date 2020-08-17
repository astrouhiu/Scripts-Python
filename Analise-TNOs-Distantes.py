'''
Esse script realiza a analise para os TNOs distantes e cria os seguintes arquivos:

- Na pasta_compilacao: entrada-hist.mp4
	      	       estatistica.txt
		       tabela1.txt
		       s.txt
		       entrada-conf-vs-t.pdf
	               entrada-N1-N2-vs-t.pdf
		       entrada-Dangles-vs-t.pdf
		       entrada-N1-Nrandom-vs-t.pdf
		       entrada-s-vs-t.pdf
		    
- Nas respectivas pastas dos arquivos de entrada: particles.txt
		 	                           arquivo_10Myr_iv4_i.txt
'''

from matplotlib.ticker import AutoMinorLocator
from multiprocessing import Pool
from itertools import combinations
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math
import random
import os
import time
# Librerias criadas com codigos Fortran do Swift
import xv2el
import el2xv

rc('text', usetex = True)
rc('font', size = 14)
rc('legend', fontsize = 10)
rc('font', **{'weight': 'light','serif':['Times']})

def positive(angle):

	# angle é um elemento e está em radianos
	
	if (angle < 0):

		angle = angle + 2*math.pi

	return angle

def detialpha(e):

	if e >= 0 and e < 1:

		# elliptic orbit
		
		ialpha = -1

	else: 

		if e > 1:

			# hyperbolic orbit
			
			ialpha = 1
	
		else:

			# parabolic orbit (e == 1)
			
			ialpha = 0

	return ialpha

def coord_h2b(nbod, mpl, xh, yh, zh, vxh, vyh, vzh):

	xb = []
	yb = []
	zb = []
	vxb = []
	vyb = []
	vzb = []

	msys = gm
	xtmp = 0.0
	ytmp = 0.0
	ztmp = 0.0
	vxtmp = 0.0
	vytmp = 0.0
	vztmp = 0.0

	for n in range(0, 4):

		msys = msys + mpl[n]
		xtmp = xtmp + mpl[n]*xh[n]
		ytmp = ytmp + mpl[n]*yh[n]
		ztmp = ztmp + mpl[n]*zh[n]
		vxtmp = vxtmp + mpl[n]*vxh[n]
		vytmp = vytmp + mpl[n]*vyh[n]
		vztmp = vztmp + mpl[n]*vzh[n]

	xbsun = -xtmp/msys
	ybsun = -ytmp/msys
	zbsun = -ztmp/msys
	vxbsun = -vxtmp/msys
	vybsun = -vytmp/msys
	vzbsun = -vztmp/msys

	for n in range(0, nbod):

		xb.append(xh[n] + xbsun)
		yb.append(yh[n] + ybsun)
		zb.append(zh[n] + zbsun)
		vxb.append(vxh[n] + vxbsun)
		vyb.append(vyh[n] + vybsun)
		vzb.append(vzh[n] + vzbsun)

	return xb, yb, zb, vxb, vyb, vzb

def rotation(inc, alpha, x, y, z, vx, vy, vz):

	# matriz de rotacao
	
	D11 = math.cos(alpha)
	D12 = math.sin(alpha)
	D13 = 0.0
	D21 = -math.cos(inc)*math.sin(alpha)
	D22 = math.cos(inc)*math.cos(alpha)
	D23 = math.sin(inc)
	D31 = math.sin(inc)*math.sin(alpha)
	D32 = -math.sin(inc)*math.cos(alpha)
	D33 = math.cos(inc)

	# equacoes de transformacao
	
	xi = D11*x + D12*y + D13*z
	yi = D21*x + D22*y + D23*z
	zi = D31*x + D32*y + D33*z
	vxi = D11*vx + D12*vy + D13*vz
	vyi = D21*vx + D22*vy + D23*vz
	vzi = D31*vx + D32*vy + D33*vz

	return xi, yi, zi, vxi, vyi, vzi	

def zero360(angle):

	if np.isnan(angle):

		return float('nan')

	nper = int(angle/360)	
	angle = angle - nper*360

	if (angle < 0):

		angle = angle + 360

	return angle

def function_media(angle):

	# Evitamos elementos float('nan') para usar funções min e max (para tabela1.txt)
	
	angle = [i for i in angle if str(i) != 'nan']

	# Esta funcao encontra a media, desvio e confinamento de angle (um vetor consistente de angulos em graus, pelo que temos que ser
	# cuidadosos ao momento de obter a media ou desvio, usar funcoes 'mean' ou 'std' dariam valores errados). 

	angle = sorted(angle)
	num = len(angle)

	# Adicionamos mais um elemento  

	angle.append(angle[0] + 360)

	# Encontrarmos a max. diferenca entre os valores contiguos de angle

	dif = []

	for i in range(0, num):

		dif.append(angle[i+1] - angle[i])

	# Tamanho do confinamento em angle

	confangle = 360 - max(dif)

	# Obtemos o indice correspondente ao valor max. de dif

	ind = np.argmax(dif)

	# Fazemos uma varredura em angle e ordenamos de forma que todos os valores (exceto o primeiro elemento, 
	# porque ja estamos considerando-lo no elemento que foi adicionado) estejam por acima do valor de angle[ind].

	a = angle[ind]

	for i in range(1, len(angle)):

		if (angle[i] <= a):

			angle[i] = angle[i] + 360 

	media = sum(angle[1:])/num

	desvio = math.sqrt(sum([(x - media)**2 for x in angle[1:]])/num)

	# Reduzimos a media no intervalo de 0 a 360 graus         

	media = zero360(media)		

	return media, desvio, confangle	

def function_hist2(a, b, c, ini, final, j, entrada, t, capom, NDcapompeak1, NDcapompeak2, NDomegapeak1, NDomegapeak2, NDobarpeak1, NDobarpeak2, texto):

	pasta = rota_comum + entrada + '/'

	if len(a) == 0:

		f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
		ax1.plot([])
		ax2.plot([])
		ax3.plot([])
		f.subplots_adjust(hspace = 0.)
		plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
		ax3.set_xlabel('angle (deg)')
		ax1.set_ylabel('$N$ in $\Delta\Omega$')
		ax2.set_ylabel('$N$ in $\Delta\omega$')
		ax3.set_ylabel('$N$ in $\Delta\\varpi$')
		ax1.set_ylim(0, 1)
		ax1.set_xlim(0, 360)
		ax1.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
		ax2.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
		ax3.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
		ax3.set_xticks(np.linspace(0, 360, 7))
		ax1.text(180, 1.12, texto, fontsize = 14, ha = 'center', va = 'center')
		ax1.text(15, 0.82, 't = %4.2E yr' %t, size = 10)
		ax1.text(15, 0.62, 'ntp = %d' %(len(capom) - 5), size = 10)
		ax1.text(15, 0.42, 'ntp$_{range}$= %d' %len(a), size = 10)
		ax1.text(300, 0.81, '$N_1$= %.3f' %NDcapompeak1, size = 10)
		ax1.text(300, 0.61, '$N_2$= %.3f' %NDcapompeak2, size = 10)
		ax2.text(300, 0.81, '$N_1$= %.3f' %NDomegapeak1, size = 10)
		ax2.text(300, 0.61, '$N_2$= %.3f' %NDomegapeak2, size = 10)
		ax3.text(300, 0.81, '$N_1$= %.3f' %NDobarpeak1, size = 10)
		ax3.text(300, 0.61, '$N_2$= %.3f' %NDobarpeak2, size = 10)
		ax1.grid(which = 'both', axis = 'x', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax1.grid(which = 'major', axis = 'y', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax2.grid(which = 'both', axis = 'x', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax2.grid(which = 'major', axis = 'y', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax3.grid(which = 'both', axis = 'x', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax3.grid(which = 'major', axis = 'y', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax1.minorticks_on()
		ax2.minorticks_on()
		ax3.minorticks_on()
		ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'both')	
		ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
		ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
		ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
		f.savefig(pasta + 'hist' + str(j) + '.jpg')
		plt.close()

	else:

		# 12 barras de bins
		
		bins = np.linspace(ini, final, 13)

		# para normalizar a 1
		
		weightsa = np.ones_like(a)/float(len(a))

		f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
		ax1.hist(a, bins, weights = weightsa, color = 'cornflowerblue')
		ax2.hist(b, bins, weights = weightsa, color = 'cornflowerblue')
		ax3.hist(c, bins, weights = weightsa, color = 'cornflowerblue')
		f.subplots_adjust(hspace = 0.)
		plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
		ax3.set_xlabel('angle (deg)')
		ax1.set_ylabel('$N$ in $\Delta\Omega$')
		ax2.set_ylabel('$N$ in $\Delta\omega$')
		ax3.set_ylabel('$N$ in $\Delta\\varpi$')
		ax1.set_ylim(0, 1)
		ax1.set_xlim(0, 360)
		ax1.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
		ax2.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
		ax3.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
		ax3.set_xticks(np.linspace(0, 360, 7))
		ax1.text(180, 1.12, texto, fontsize = 14, ha = 'center', va = 'center')
		ax1.text(15, 0.82, 't = %4.2E yr' %t, size = 10)
		ax1.text(15, 0.62, 'ntp = %d' %(len(capom) - 5), size = 10)
		ax1.text(15, 0.42, 'ntp$_{range}$= %d' %len(a), size = 10)
		ax1.text(300, 0.81, '$N_1$= %.3f' %NDcapompeak1, size = 10)
		ax1.text(300, 0.61, '$N_2$= %.3f' %NDcapompeak2, size = 10)
		ax2.text(300, 0.81, '$N_1$= %.3f' %NDomegapeak1, size = 10)
		ax2.text(300, 0.61, '$N_2$= %.3f' %NDomegapeak2, size = 10)
		ax3.text(300, 0.81, '$N_1$= %.3f' %NDobarpeak1, size = 10)
		ax3.text(300, 0.61, '$N_2$= %.3f' %NDobarpeak2, size = 10)
		ax1.grid(which = 'both', axis = 'x', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax1.grid(which = 'major', axis = 'y', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax2.grid(which = 'both', axis = 'x', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax2.grid(which = 'major', axis = 'y', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax3.grid(which = 'both', axis = 'x', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax3.grid(which = 'major', axis = 'y', color = 'darkgrey', linestyle = '--', linewidth = 0.5)
		ax1.minorticks_on()
		ax2.minorticks_on()
		ax3.minorticks_on()
		ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'both')	
		ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
		ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
		ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
		f.savefig(pasta + 'hist' + str(j) + '.jpg')
		plt.close()
		
def function_hist(a, ini, final):

	# A diferenca entre os primeiro e segundo (se existir) picos máximos deve ser maior igual do que 60 graus. 

	# 12 barras de bins
	
	bins = np.linspace(ini, final, 13)

	# para normalizar a 1
	
	weightsa = np.ones_like(a)/float(len(a))

	# x: posições para construir os bins, 13 elementos
	# y: freq. em 'a', 12 elementos 
	
	y, x = np.histogram(a, bins, weights = weightsa)

	# Criamos um novo vetor 'y2', onde adicionamos a 'y' seus valores extremos, nos extremos opuestos
	
	yini = y[-1]
	yfinal = y[0]
	y2 = [yini] + y.tolist() + [yfinal]

	# Comparamos y2[i] com seus adjacentes só para os elementos de 'y', se ele é maior do que seus adjacentes o anexamos no vetor 
	# 'maxy' e anexamos o índice correspondente de 'x'. Na maioria das vezes teremos pelo menos um elemento de maxy (é muito 
	# improvável que todos os bins tenham o mesmo peso).
	
	maxy = []
	indx = []

	for j in range(0, len(y2) - 2):

		yizq = y2[j]
		ydir = y2[j + 2]

		if y2[j + 1] > yizq and y2[j + 1] > ydir:

			maxy.append(y2[j + 1])
			indx.append(j)

	if maxy != []:

		# Primeiro máximo. Indice da última ocorrencia de max(maxy) em maxy, ymax1	
		
		ind = [i for i, n in enumerate(maxy) if n == max(maxy)][-1]
		ymax1 = maxy[ind]
		arg_max1 = indx[ind]
		xmax1 = (x[arg_max1 + 1] - x[arg_max1])/2 + x[arg_max1]

		# Segundo máximo
		
		if len(maxy) > 1:

			# Se o valor máximo em maxy tem duplicidade
			
			if sorted(maxy)[-1] == sorted(maxy)[-2]:

				# Indice da penultima ocorrencia de max(maxy) em maxy, ymax2
				
				ind = [i for i, n in enumerate(maxy) if n == max(maxy)][-2]
				ymax2 = maxy[ind]
				arg_max2 = indx[ind]
				xmax2 = (x[arg_max2 + 1] - x[arg_max2])/2 + x[arg_max2]

			else:

				# Dá a primeira ocorrencia em caso de repeticao do segundo máximo
				
				ymax2 = sorted(maxy)[-2]
				arg_max2 = indx[maxy.index(ymax2)]
				xmax2 = (x[arg_max2 + 1] - x[arg_max2])/2 + x[arg_max2]

		else:

			ymax2 = float('nan')
			xmax2 = float('nan')

	else:

		xmax1 = float('nan')
		ymax1 = float('nan')
		xmax2 = float('nan')
		ymax2 = float('nan')

	return xmax1, xmax2, ymax1, ymax2

def function_est(lista):

	# Evitamos elementos float('nan') para usar funções min e max
	
	lista = [i for i in lista if str(i) != 'nan']
	fmaxlista = "%.3f"%max(lista)
	fmedlista = "%.3f"%((max(lista) - min(lista))/2 + min(lista))

	return fmaxlista, fmedlista

def AnaliseProjeto(argumentos):

	entrada, texto = argumentos

	# Rota à pasta que contém os arquivos de entrada e onde se almacenarão os arquivos criados no iv4
	
	pasta = rota_comum + entrada + '/'

	# Primeira Parte
	# ===============================================================================================================================
	'''
	Essa parte do codigo pega arquivos no formado do follow_swift-all.x (id a e inc capom obar lamb) gerados anteriormente em bash e 
	retorna os arquivos no iv4 no formato: id a e inc capom omega capm, para cada certo tempo contendo dados de planetas e partículas
	(reiteradamente ou nao, dependendo da corrida)
	'''

	# Mudamos dados dos arquivos de entrada que estão em relação ao sistema iv4 inicial heliocêntrico para o sistema iv4 instantáneo
	# heliocêntrico, definido pelo vetor momento angular dos planetas gigantes conhecidos mediante duas rotações (ângulos 'incplano' 
	# e 'alpha') 

	for j in range(1, narq + 1):

		print('j =', j)

		if entrada == 'a700-10-0-100-3':

			data = open(pasta + arquivo_entrada + str(j) + '.txt', 'r')
			lines1 = data.read().splitlines()
			data.close()

			data = open(rota_comum + 'a700-10-0-100-3-II/' + arquivo_entrada + str(j) + '.txt', 'r')
			lines2 = data.read().splitlines()
			data.close()

			lines = np.concatenate((lines1, lines2[5:]))

		else:
	
			data = open(pasta + arquivo_entrada + str(j) + '.txt', 'r')
			lines = data.read().splitlines()
			data.close()

		outputfile = open(pasta + arquivo_invariante + str(j) + '.txt', 'w')

		counterlines = 1

		for line in lines:

			text = line.split()

			if text[0] == '2' or text[0] == '3' or text[0] == '4' or text[0] == '5':

				if text[0] == '2':

					hx = 0.0
					hy = 0.0
					hz = 0.0
			
					n = 0 
	
				# Calculamos o momento angular resultante dos 4 planetas gigantes
					
				ind = text[0]
				a = float(text[1])
				e = float(text[2])
				inc = float(text[3])*math.pi/180
				capom = float(text[4])*math.pi/180
				obar = float(text[5])*math.pi/180
				lamb = float(text[6])*math.pi/180
				omega = obar - capom
				capm = lamb - obar
				omega = positive(omega)
				capm = positive(capm)
				ialpha = detialpha(e)

				# el2xv. Essa função precisa que os ângulos estem em radianos
					
				x, y, z, vx, vy, vz = el2xv.orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm)

				hx_m = np.float32(np.float32(mpl[n])*(np.float32(y)*np.float32(vz) - np.float32(z)*np.float32(vy)))
				hy_m = np.float32(np.float32(mpl[n])*(np.float32(z)*np.float32(vx) - np.float32(x)*np.float32(vz)))
				hz_m = np.float32(np.float32(mpl[n])*(np.float32(x)*np.float32(vy) - np.float32(y)*np.float32(vx)))

				hx = np.float32(hx + hx_m)
				hy = np.float32(hy + hy_m)
				hz = np.float32(hz + hz_m)

				n = n + 1

				if text[0] == '5':

					h2 = np.float32(hx*hx + hy*hy + hz*hz)
					h = np.float32(math.sqrt(h2))

					if (hz > h):
		
						hz = h	
						hx = 0.0
						hy = 0.0

					# math.acos(x): Return the arc cosine of x, in radians
						
					incplano = np.float32(math.acos(hz/h))
					fac = np.float32(math.sqrt(hx*hx + hy*hy)/h)
	
					if (fac < 2**-127 or incplano == 0):

						alpha = 0.0

					else:

						# math.atan2(y, x): Return atan(y/x), in radians. The result is between -pi and pi. The 
						# vector in the plane from the origin to point (x, y) makes this angle with the positive 
						# X axis. The point of atan2() is that the signs of both inputs are known to it, so it 
						# can compute the correct quadrant for the angle. For example, atan(1) and atan2(1, 1) 
						# are both pi/4, but atan2(-1, -1) is -3*pi/4
							
						alpha = np.float32(math.atan2(hx, -hy))
	
					if (alpha < 0):

						alpha = np.float32(alpha + 2*math.pi)
					
					flag = 1
			else:
				
				if text[0] == '6' or (entrada == 'mig-sempl9' and flag == 1):			

					nn = 0

					for x in lines[counterlines - 5: counterlines]:

						text = x.split()
						ind = text[0]
						a = float(text[1])
						e = float(text[2])
						inc = float(text[3])*math.pi/180
						capom = float(text[4])*math.pi/180
						obar = float(text[5])*math.pi/180
						lamb = float(text[6])*math.pi/180
						omega = obar - capom
						capm = lamb - obar
						omega = positive(omega)
						capm = positive(capm)
						ialpha = detialpha(e)
	
						# el2xv. Essa função precisa que os ângulos estem em radianos
	
						x, y, z, vx, vy, vz = el2xv.orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm)

						xi, yi, zi, vxi, vyi, vzi = rotation(np.float32(incplano), np.float32(alpha), np.float32(x), np.float32(y), np.float32(z), np.float32(vx), np.float32(vy), np.float32(vz))	

						# xv2el. Essa função retorna os ângulos em radianos

						gmsum = gm + mpl[nn]

						_, ai, ei, inci, capomi, omegai, capmi = xv2el.orbel_xv2el(xi, yi, zi, vxi, vyi, vzi, gmsum)
				
						outputfile.write(ind.ljust(8) + "  " + str("%6E" %ai) + "  " + str("%6E" %ei) + "  " + str("%6E" %(inci*180/math.pi)) + "  " + str("%6E" %(capomi*180/math.pi)) + "  " + str("%6E" %(omegai*180/math.pi)) + "  " + str("%6E" %(capmi*180/math.pi)) + "\n")	
		
						nn = nn + 1
							
					flag = 0

				else:

					if (line[10] != '*'):

						ind = text[0]
						a = float(text[1])
						e = float(text[2])
						inc = float(text[3])*math.pi/180
						capom = float(text[4])*math.pi/180
						obar = float(text[5])*math.pi/180
						lamb = float(text[6])*math.pi/180
						omega = obar - capom
						capm = lamb - obar
						omega = positive(omega)
						capm = positive(capm)
						ialpha = detialpha(e)
	
						# Não consideramos órbitas hiperbólicas
							
						if ialpha != 1:
	
							# el2xv. Essa função precisa que os ângulos estem em radianos
								
							x, y, z, vx, vy, vz = el2xv.orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm)

							xi, yi, zi, vxi, vyi, vzi = rotation(np.float32(incplano), np.float32(alpha), np.float32(x), np.float32(y), np.float32(z), np.float32(vx), np.float32(vy), np.float32(vz))	

							# xv2el. Essa função retorna os ângulos em radianos

							gmsum = gm
								
							_, ai, ei, inci, capomi, omegai, capmi = xv2el.orbel_xv2el(xi, yi, zi, vxi, vyi, vzi, gmsum)
				
							outputfile.write(ind.ljust(8) + "  " + str("%6E" %ai) + "  " + str("%6E" %ei) + "  " + str("%6E" %(inci*180/math.pi)) + "  " + str("%6E" %(capomi*180/math.pi)) + "  " + str("%6E" %(omegai*180/math.pi)) + "  " + str("%6E" %(capmi*180/math.pi)) + "\n")	

			counterlines = counterlines + 1

		outputfile.close()
	# ===============================================================================================================================

	# Segunda Parte
	# ===============================================================================================================================
	'''
	Nessa parte do código se analisam os dados gerados anteriormente (arquivos no iv4) e os resultados são salvos em gráficos
	'''

	if executeHist or executeMean:

		t = 0
		time = []

		Dcapompeak1 = []
		Dcapompeak2 = []
		Domegapeak1 = []
		Domegapeak2 = []
		Dobarpeak1 = []
		Dobarpeak2 = []
		NDcapompeak1 = []
		NDcapompeak2 = []
		NDomegapeak1 = []
		NDomegapeak2 = []
		NDobarpeak1 = []
		NDobarpeak2 = []

		flag = 0
		num = []
		N = []
		s_capom = []
		s_omega = []
		s_obar = []

		# Arquivo 'nND': número de pontos n, N determinado a partir de uma distribuição randômica com n pontos, desvio padrão 
		# referente ao N. Analise: Ver em quantos sigmas o valor de N1 difere da média em nND. Se o valor de N1 for menor que a 
		# média, considere zero

		n_random, N_random, D_random = np.loadtxt(rota_nND, usecols = (0, 1, 2), unpack = True)
		n_random = n_random.tolist()	
		N_random = N_random.tolist()	
		D_random = D_random.tolist()

		plt.errorbar(n_random, N_random, yerr = D_random, alpha = 0.5, color = 'lightcoral', label = 'errorbar')
		plt.plot(n_random, N_random, marker = 'o', markersize = 1.5, linestyle = '', color = 'black', alpha = 0.5, label = 'N')
		plt.plot(n_random, [x + y for x, y in zip(N_random, D_random)], marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4, alpha = 0.5, label = 'D')
		plt.plot(n_random, [x - y for x, y in zip(N_random, D_random)], marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4, alpha = 0.5)
		plt.xlabel('n')
		plt.ylabel('N')
		plt.legend()
		plt.show()
		
		mediamediacapom = []
		mediamediaomega = []
		mediamediaobar = []
		desviomediacapom = []
		desviomediaomega = []
		desviomediaobar = []
		mediaconfcapom = []
		mediaconfomega = []
		mediaconfobar = []
		desvioconfcapom = []
		desvioconfomega = []
		desvioconfobar = []

		for j in range(1, narq + 1):

			print('j =', j)

			ind = np.genfromtxt(pasta + arquivo_invariante + str(j) + '.txt', dtype = str, usecols = (0), unpack = True)
			a, e, inc, capom, omega = np.loadtxt(pasta + arquivo_invariante + str(j) + '.txt', usecols = (1, 2, 3, 4, 5), unpack = True)
			q = np.multiply(a, [1 - x for x in e])
			obar = [sum(x) for x in zip(omega, capom)]			
			obar = list(map(lambda x: zero360(x), obar))

			if j == 1:

				a_CI = [i for i, j in zip(a, ind) if j[0] == 'p']
				q_CI = [i for i, j in zip(q, ind) if j[0] == 'p']
				inc_CI = [i for i, j in zip(inc, ind) if j[0] == 'p']
				capom_CI = [i for i, j in zip(capom, ind) if j[0] == 'p']
				omega_CI = [i for i, j in zip(omega, ind) if j[0] == 'p']
				obar_CI = [i for i, j in zip(obar, ind) if j[0] == 'p']

				print('*****************************************************************************************')
				print('Entrada:', entrada)
				print('C.I. para as partículas:')
				print('min(a), max(a) em au =', min(a_CI), max(a_CI))
				print('min(q), max(q) em au =', min(q_CI), max(q_CI))
				print('min(inc), max(inc) em graus =', min(inc_CI), max(inc_CI))
				print('min(capom), max(capom) em graus =', min(capom_CI), max(capom_CI))
				print('min(omega), max(omega) em graus =', min(omega_CI), max(omega_CI))
				print('min(obar), max(obar) em graus =', min(obar_CI), max(obar_CI))
				print('*****************************************************************************************')

			time.append(t)

			# Desde aqui se tem em conta só as partículas pertencentes a intervalos em 'a', 'q' e 'inc' similar à dos TNOs 
			# distantes observados
			
			ind_tp = [l for i, j, k, l in zip(a, q, inc, ind) if i >= aminobs and i <= amaxobs and j >= qminobs and j <= qmaxobs and k >= incminobs and k <= incmaxobs and l != '6']
			a_tp = [i for i, j, k, l in zip(a, q, inc, ind) if i >= aminobs and i <= amaxobs and j >= qminobs and j <= qmaxobs and k >= incminobs and k <= incmaxobs and l != '6']
			capom_tp = [m for i, j, k, l, m in zip(a, q, inc, ind, capom) if i >= aminobs and i <= amaxobs and j >= qminobs and j <= qmaxobs and k >= incminobs and k <= incmaxobs and l != '6']
			omega_tp = [m for i, j, k, l, m in zip(a, q, inc, ind, omega) if i >= aminobs and i <= amaxobs and j >= qminobs and j <= qmaxobs and k >= incminobs and k <= incmaxobs and l != '6']
			obar_tp = [m for i, j, k, l, m in zip(a, q, inc, ind, obar) if i >= aminobs and i <= amaxobs and j >= qminobs and j <= qmaxobs and k >= incminobs and k <= incmaxobs and l != '6']
			
			capom9 = [i for i, j in zip(capom, ind) if j == '6']
			omega9 = [i for i, j in zip(omega, ind) if j == '6']
			obar9 = [i for i, j in zip(obar, ind) if j == '6']
			
			Dcapom_tp = []
			Domega_tp = []
			Dobar_tp = []

			# arquivo_10Myr_1.txt da corrida com GT 'a3000-e0.7-i30-w0-O150' contem os ids: 2, 3, 4, 5, 6, p1, p2, ..., 
			# p9200, 2, 3, 4, 5, 6, p9201, ..., p12400

			if entrada == 'a3000-e0.7-i30-w0-O150':
	
				Dcapom_tp.extend([capom9[0] - i if float(j[1:]) >= 1 and float(j[1:]) <= 9200 else capom9[1] - i for i, j in zip(capom_tp, ind_tp)])
				Domega_tp.extend([omega9[0] - i if float(j[1:]) >= 1 and float(j[1:]) <= 9200 else omega9[1] - i for i, j in zip(omega_tp, ind_tp)])
				Dobar_tp.extend([obar9[0] - i if float(j[1:]) >= 1 and float(j[1:]) <= 9200 else obar9[1] - i for i, j in zip(obar_tp, ind_tp)])

			# arquivo_10Myr_1.txt das corridas com GT contem os ids: 2, 3, 4, 5, 6, p1, p2, ..., p100000; excetuando a 
			# anterior e 'GT-a2000-30-350-400-5', essa última contem: 2, 3, 4, 5, 6, p1, p2, ..., p10000
			# arquivo_10Myr_1.txt das corridas do artigo contem os ids: 2, 3, 4, 5, 6, p1, p2, ..., p10000

			else:

				Dcapom_tp.extend([capom9[0] - i for i, j in zip(capom_tp, ind_tp)])
				Domega_tp.extend([omega9[0] - i for i, j in zip(omega_tp, ind_tp)])
				Dobar_tp.extend([obar9[0] - i for i, j in zip(obar_tp, ind_tp)])

			Dcapom_tp = list(map(lambda x: zero360(x), Dcapom_tp))
			Domega_tp = list(map(lambda x: zero360(x), Domega_tp))
			Dobar_tp = list(map(lambda x: zero360(x), Dobar_tp))

			# ---------------------------------------------------------------------------------------------------------------
			# Em cada tempo geramos histogramas em capom, omega, obar das partículas que pertencem ao intervalo similar ao 
			# observado, em relacao aos do planeta 9. Podemos ter plots vazios e dados com float('nan') se não temos 
			# partículas pertencentes a essas faixas ou se não satisfazem o critério para calcular o primeiro e/ou segundo 
			# pico nos histogramas
	
			if executeHist: 

				num.append(len(capom_tp))

				Dcapomhist1 = []
				Dcapomhist2 = []
				NDcapomhist1 = []
				NDcapomhist2 = []
				Domegahist1 = []
				Domegahist2 = []
				NDomegahist1 = []
				NDomegahist2 = []
				Dobarhist1 = []
				Dobarhist2 = []
				NDobarhist1 = []
				NDobarhist2 = [] 			

				# .......................................................................................................
				# Amostra
				
				if entrada == 'a1500-30-0-60-5' and j == 219:

					colors = cm.coolwarm(np.linspace(0, 1, 30))
					f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
					b = 7					
				# .......................................................................................................

				for z in range(0, 30):

					valueini = z
	
					Dcapomhist = [x + 360 if x < valueini else x for x in Dcapom_tp] 
					Domegahist = [x + 360 if x < valueini else x for x in Domega_tp]
					Dobarhist = [x + 360 if x < valueini else x for x in Dobar_tp]

					# ...............................................................................................
					# Amostra
					
					if entrada == 'a1500-30-0-60-5' and j == 219 and (z == 0 or z == 12 or z == 22):

						# 12 barras de bins
						
						bins = np.linspace(z, z + 360, 13)

						# normalizar a 1
						
						weightsa = np.ones_like(Dcapomhist)/float(len(Dcapomhist))

						ax1.hist(Dcapomhist, bins, weights = weightsa, color = colors[b], histtype = 'step', stacked = True, fill = False, label = entrada)
						ax2.hist(Domegahist, bins, weights = weightsa, color = colors[b], histtype = 'step', stacked = True, fill = False)
						ax3.hist(Dobarhist, bins, weights = weightsa, color = colors[b], histtype = 'step', stacked = True, fill = False)
						f.subplots_adjust(hspace = 0.15)
						plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
						ax3.set_xlabel('$\Delta$angle (deg)')
						ax1.set_ylabel('$N$ in $\Delta\Omega$')
						ax2.set_ylabel('$N$ in $\Delta\omega$')
						ax3.set_ylabel('$N$ in $\Delta\\varpi$')
						ax1.set_xlim(0, 390)
						b = b + 8
						ax1.minorticks_on()
						ax2.minorticks_on()
						ax3.minorticks_on()
						ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')	
						ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'major')
						ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'major')
						ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')	
						ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True,  which = 'minor')
						ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True,  which = 'minor')

						if z == 22:

							ax1.arrow(165, 0.32, 0, 0.15, head_width = 3, head_length = 0.05, fc = 'gray', ec = 'gray', linewidth = 0.5)
							ax1.arrow(315, 0.2, 0, 0.15, head_width = 3, head_length = 0.05, fc = 'gray', ec = 'gray', linewidth = 0.5)
							ax2.arrow(285, 0.32, 0, 0.15, head_width = 3, head_length = 0.05, fc = 'gray', ec = 'gray', linewidth = 0.5)
							ax2.arrow(15, 0.13, 0, 0.15, head_width = 3, head_length = 0.05, fc = 'gray', ec = 'gray', linewidth = 0.5)
							ax3.arrow(87, 0.53, 0, 0.15, head_width = 3, head_length = 0.05, fc = 'gray', ec = 'gray', linewidth = 0.5)
							ax3.arrow(267, 0.32, 0, 0.15, head_width = 3, head_length = 0.05, fc = 'gray', ec = 'gray', linewidth = 0.5)
							ax1.text(0.5, 1.2, texto, fontsize = 14, ha = 'center', va = 'center', transform = ax1.transAxes)
							ax1.text(165, 0.56, '$N_1$', size = 12, ha = 'center')
							ax1.text(315, 0.44, '$N_2$', size = 12, ha = 'center')
							ax2.text(285, 0.56, '$N_1$', size = 12, ha = 'center')
							ax2.text(15, 0.37, '$N_2$', size = 12, ha = 'center')
							ax3.text(87, 0.77, '$N_1$', size = 12, ha = 'center')
							ax3.text(267, 0.56, '$N_2$', size = 12, ha = 'center')
							ax1.set_xticks([0, 50, 100, 150, 200, 250, 300, 350])
							ax1.set_yticks([0.0, 0.5, 1.0])
							ax1.set_xlim(0, 390)
							ax1.set_ylim(0, 1)
							f.savefig(pasta_compilacao + 'amostra_hist.pdf', format = 'pdf', bbox_inches = 'tight')
							plt.close()
					# ...............................................................................................

					# Encontramos as frequências principais e secundárias e os correspondentes ângulos relativos, nos 
					# histogramas	
					
					xmax1, xmax2, ymax1, ymax2 = function_hist(Dcapomhist, z, 360 + z) 
					Dcapomhist1.append(xmax1)
					Dcapomhist2.append(xmax2)
					NDcapomhist1.append(ymax1)
					NDcapomhist2.append(ymax2)

					xmax1, xmax2, ymax1, ymax2 = function_hist(Domegahist, z, 360 + z) 
					Domegahist1.append(xmax1)
					Domegahist2.append(xmax2)
					NDomegahist1.append(ymax1)
					NDomegahist2.append(ymax2)

					xmax1, xmax2, ymax1, ymax2 = function_hist(Dobarhist, z, 360 + z) 
					Dobarhist1.append(xmax1)
					Dobarhist2.append(xmax2)
					NDobarhist1.append(ymax1)
					NDobarhist2.append(ymax2)

				# Se não existe a frequência principal ou secundária, ela e o respectivo elemento angular relativo são
				# preenchidos com float('nan'). Isto permitiria, ao final, plotar sem problemas os ângulos relativos e
				# frequências vs. o tempo; ao manterem a mesma longitude que a do vetor tempo (python omitiria esse 
				# ponto)

				# Pegamos o máximo valor do conjunto de frequências principais, as respectivas! frequências secundárias 
				# e os respectivos ângulos relativos

				Nmax = max(NDcapomhist1)
				arg = NDcapomhist1.index(Nmax)
				Dcapompeak1.append(Dcapomhist1[arg])
				Dcapompeak2.append(Dcapomhist2[arg])
				NDcapompeak1.append(Nmax)
				NDcapompeak2.append(NDcapomhist2[arg])

				Nmax = max(NDomegahist1)
				arg = NDomegahist1.index(Nmax)
				Domegapeak1.append(Domegahist1[arg])
				Domegapeak2.append(Domegahist2[arg])
				NDomegapeak1.append(Nmax)				
				NDomegapeak2.append(NDomegahist2[arg])

				Nmax = max(NDobarhist1)
				arg = NDobarhist1.index(Nmax)
				Dobarpeak1.append(Dobarhist1[arg])
				Dobarpeak2.append(Dobarhist2[arg])	
				NDobarpeak1.append(Nmax)
				NDobarpeak2.append(NDobarhist2[arg])
			
				# Pegamos ou determinamos por extrapolação o valor da frequência N correspondente a um certo núm. de 
				# pontos (igual ao número das partículas simuladas pertencendo ao intervalo similar ao observado) 
				# aleatórios desde o arquivo nND e determinamos em quantos sigmas o valor de N1 (calculado acima) difere 
				# da média em nND (N).

				# Se a longitude é um dos elementos da lista n_random
				
				if len(capom_tp) in n_random:

					k = n_random.index(len(capom_tp))
					N.append(N_random[k])
		
					if NDcapompeak1[-1] <= N[-1]:
						s_capom.append(0)
					else:
						s_capom.append((NDcapompeak1[-1] - N[-1])/D_random[k])

					if NDomegapeak1[-1] <= N[-1]:
						s_omega.append(0)
					else:
						s_omega.append((NDomegapeak1[-1] - N[-1])/D_random[k])

					if NDobarpeak1[-1] <= N[-1]:
						s_obar.append(0)
					else:
						s_obar.append((NDobarpeak1[-1] - N[-1])/D_random[k])

				else:

					# Obtemos o número mais próximo no nND ao nosso número e um número vizinho para extrapolar
					
					xfit1 = min(n_random, key = lambda x: abs(x - len(capom_tp)))

					# Minimo valor em n_random é ninf_nND
					
					if len(capom_tp) < ninf_nND:

						xfit2 = n_random[n_random.index(xfit1) + 1]

					# Máximo valor em n_random é nsup_nND
					
					elif len(capom_tp) > nsup_nND:

						xfit2 = n_random[n_random.index(xfit1) - 1]

					else:

						# Pegar vizinho da esquerda ou direita
						
						if len(capom_tp) < xfit1:
						
							xfit2 = n_random[n_random.index(xfit1) - 1]
							
						else:
							xfit2 = n_random[n_random.index(xfit1) + 1]

					# n_random
					
					xfit = [xfit1, xfit2]

					# N_random
					
					yfit = [N_random[n_random.index(xfit1)], N_random[n_random.index(xfit2)]]
					fit = np.polyfit(xfit, yfit, 1)
					fit_fn = np.poly1d(fit)

					N.append(fit_fn(len(capom_tp)))

					# D_random
					
					yfit = [D_random[n_random.index(xfit1)], D_random[n_random.index(xfit2)]]
					fit = np.polyfit(xfit, yfit, 1)
					fit_fn = np.poly1d(fit)
					
					D = fit_fn(len(capom_tp))

					# s
					
					if NDcapompeak1[-1] <= N[-1]:
					
						s_capom.append(0)
						
					else:
					
						s_capom.append((NDcapompeak1[-1] - N[-1])/D)

					if NDomegapeak1[-1] <= N[-1]:
					
						s_omega.append(0)
						
					else:
					
						s_omega.append((NDomegapeak1[-1] - N[-1])/D)

					if NDobarpeak1[-1] <= N[-1]:
					
						s_obar.append(0)
						
					else:
					
						s_obar.append((NDobarpeak1[-1] - N[-1])/D)
						
				# Realizamos histogramas a cada tempo só no primeiro loop (início do primeiro bin em 0 graus), mas as
				# lendas correspondem às máximas frequências principais e as respectivas frequências secundárias 
				# considerando todo o loop

				function_hist2(Dcapom_tp, Domega_tp, Dobar_tp, 0, 360, j, entrada, t, capom, NDcapompeak1[-1], NDcapompeak2[-1], NDomegapeak1[-1], NDomegapeak2[-1], NDobarpeak1[-1], NDobarpeak2[-1], texto) 

				# Se em algum tempo se tem poucas partículas (< 50) pertencendo ao intervalo similar ao observado 
				# salvamos: o número de partículas dentro do intervalo similar ao observado, o número total de partículas 
				# na simulacao, as frequências máximas considerando os 30 loops e o tempo

				if len(capom_tp) < 50:

					if flag == 0:

						particles = open(rota_comum + entrada + 'particles.txt', 'w')
						particles.write("%18s"%entrada + "\n")
						particles.write("%6s"%"ntp<50" + " " + "%6s"%"ntp" + " " + "%11s"%"Nmax_Dcapom" + " " + "%11s"%"Nmax_Domega" + " " + "%10s"%"Nmax_Dobar" + " " + "%9s"%"time (yr)" + "\n")
						particles.close()

						flag = 1

					else:

						particles = open(rota_comum + entrada + 'particles.txt', 'a+')
						particles.write(str("%6d"%len(Dcapom_tp)) + " " + str("%6d"%(len(capom) - 5)) + " " + str("%11.4f"%NDcapompeak1[-1]) + " " + str("%11.4f"%NDomegapeak1[-1]) + " " + str("%10.4f"%NDobarpeak1[-1]) + " " + str("%9.2E"%t) + "\n")
						particles.close()
			# ---------------------------------------------------------------------------------------------------------------
		
			# ---------------------------------------------------------------------------------------------------------------
			# Em cada tempo se obtém a média e o confinamento nos ângulos das partículas. Se pelo menos 'numobs' particulas 
			# existem dentro do intervalo similar ao observado em 'a', 'q' e 'inc', criamos outros novos vetores (capomtp, 
			# omegatp, obartp) de 'numobs' elementos a partir da escolha aleatória de ('numobs') índices dentre a longitude 
			# dos vetores capom_tp, omega_tp, obar_tp

			if executeMean:

				if (len(a_tp) >= numobs):

					mediacapom = []
					mediaomega = []
					mediaobar = []
					confcapom = []
					confomega = []
					confobar = []

					# Realizamos esse procedimento 'numite' vezes para obter a média e std
					
					for z in range(numite):
			
						# random.sample takes a population and a sample size k and returns k random members of      
						# the population. Essa variacao de random.sample(list, m) faz com que as listas 
						# resultantes se correspodam em índices
						
						capomtp, omegatp, obartp = zip(*random.sample(list(zip(capom_tp, omega_tp, obar_tp)), numobs))
					
						# Em cada escolha aleatória vamos determinar o menor confinamento possível nos ângulos 
						# dentre todas as combinações possíveis para um certo núm. de objetos dentre a amostra
						# (no caso de 'obar', o min. confinamento dentre as combinações possíveis de 9 objetos 
						# dentre o total (11 objetos))
						
						m1 = []
						m2 = []
						m3 = []
						c1 = []
						c2 = []
						c3 = []
	
						for c in combinations(capomtp, numobsconfcapom):	
							x1, _, x3 = function_media(c)
							m1.append(x1)
							c1.append(x3)
						
						for c in combinations(omegatp, numobsconfomega):	
							x1, _, x3 = function_media(c)
							m2.append(x1)
							c2.append(x3)
	
						for c in combinations(obartp, numobsconfobar):	
							x1, _, x3 = function_media(c)
							m3.append(x1)
							c3.append(x3)

						# Pegamos os valores que fornecem o menor confinamento
						
						mediacapom.append(m1[c1.index(min(c1))])
						mediaomega.append(m2[c2.index(min(c2))])
						mediaobar.append(m3[c3.index(min(c3))])
						confcapom.append(min(c1))
						confomega.append(min(c2))
						confobar.append(min(c3))

					# Para um certo tempo temos uma coleção de médias e confinamentos nos ângulos que fornecem o 
					# menor confinamento. Calculamos a média e desvio padrão das médias e confinamentos
					
					x1, x2, _ = function_media(mediacapom)
					mediamediacapom.append(x1) 
					desviomediacapom.append(x2)
					x1, x2, _ = function_media(mediaomega)
					mediamediaomega.append(x1)
					desviomediaomega.append(x2)
					x1, x2, _ = function_media(mediaobar)
					mediamediaobar.append(x1)
					desviomediaobar.append(x2)

					x1, x2, _ = function_media(confcapom)
					mediaconfcapom.append(x1)
					desvioconfcapom.append(x2)
					x1, x2, _ = function_media(confomega)
					mediaconfomega.append(x1)
					desvioconfomega.append(x2)
					x1, x2, _ = function_media(confobar)
					mediaconfobar.append(x1)
					desvioconfobar.append(x2)
				
				else:

					mediamediacapom.append(float('nan')) 
					desviomediacapom.append(float('nan'))
					mediamediaomega.append(float('nan'))
					desviomediaomega.append(float('nan'))
					mediamediaobar.append(float('nan'))
					desviomediaobar.append(float('nan'))
					mediaconfcapom.append(float('nan'))
					desvioconfcapom.append(float('nan'))
					mediaconfomega.append(float('nan'))
					desvioconfomega.append(float('nan'))
					mediaconfobar.append(float('nan'))
					desvioconfobar.append(float('nan'))
					
			# ---------------------------------------------------------------------------------------------------------------
	
			t = t + deltat

		# -----------------------------------------------------------------------------------------------------------------------
		
		if entrada == 'init-close-net/a3000-e0.7-i30-w0-O150':
		
			entrada = 'init-close-net-a3000-e0.7-i30-w0-O150'
		
		imag1 = pasta_compilacao + entrada + '-conf-vs-t.pdf'
		imag2 = pasta_compilacao + entrada + '-N1-N2-vs-t.pdf'
		imag3 = pasta_compilacao + entrada + '-Dangles-vs-t.pdf'
		imag4 = pasta_compilacao + entrada + '-N1-Nrandom-vs-t.pdf'
		imag5 = pasta_compilacao + entrada + '-s-vs-t.pdf'
		
		if executeHist:	

			os.chdir(pasta)
			os.system('ffmpeg -r 3 -i hist%d.jpg -vcodec libx264 -y -an -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" x2.mp4')
			os.system('mv x2.mp4 ' + pasta_compilacao + entrada + '-hist.mp4')
			os.system('rm *.jpg')

			# ---------------------------------------------------------------------------------------------------------------
			# Dados do analise estatístico (fmax e fmed), para tudo o tempo de integração e para os últimos 500 Myr 
					
			fmaxcapom1, fmedcapom1 = function_est(NDcapompeak1)	
			fmaxomega1, fmedomega1 = function_est(NDomegapeak1)
			fmaxobar1, fmedobar1 = function_est(NDobarpeak1)
			d = [i for i, n in enumerate(time) if n >= 4*10**9][0]
			fmaxcapom1d, fmedcapom1d = function_est(NDcapompeak1[d:])	
			fmaxomega1d, fmedomega1d = function_est(NDomegapeak1[d:])
			fmaxobar1d, fmedobar1d = function_est(NDobarpeak1[d:])

			est = open(pasta_compilacao + 'estatistica.txt', 'a+')
			est.write("%18s"%entrada + " & " + "%6s"%str(fmaxcapom1) + " & " + "%6s"%str(fmedcapom1) + " & " + "%6s"%str(fmaxcapom1d) + " & " + "%6s"%str(fmedcapom1d) + " & " + "%6s"%str(fmaxomega1) + " & " + "%6s"%str(fmedomega1) + " & " + "%6s"%str(fmaxomega1d) + " & " + "%6s"%str(fmedomega1d) + " & " + "%6s"%str(fmaxobar1) + " & " + "%6s"%str(fmedobar1) + " & " + "%6s"%str(fmaxobar1d) + " & " + "%6s"%str(fmedobar1d) + "\\\\" + "\n")
			est.close()		
			
			# ---------------------------------------------------------------------------------------------------------------
			# Plot de ângulos relativos correspondentes aos dos picos nas frequências nos histogramas vs. o tempo

			Dcapom = list(map(lambda x: zero360(x), Dcapompeak1))
			Dcapom2 = list(map(lambda x: zero360(x), Dcapompeak2))
			Domega = list(map(lambda x: zero360(x), Domegapeak1))
			Domega2 = list(map(lambda x: zero360(x), Domegapeak2))
			Dobar = list(map(lambda x: zero360(x), Dobarpeak1))
			Dobar2 = list(map(lambda x: zero360(x), Dobarpeak2))

			f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
			ax1.plot(time, Dcapom, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'steelblue', markersize = 1.5)
			ax1.plot(time, Dcapom2, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'lightcoral', markersize = 1.5)
			ax2.plot(time, Domega, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'steelblue', markersize = 1.5)
			ax2.plot(time, Domega2, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'lightcoral', markersize = 1.5)
			ax3.plot(time, Dobar, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'steelblue', markersize = 1.5)
			ax3.plot(time, Dobar2, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'lightcoral', markersize = 1.5)
			f.subplots_adjust(hspace = 0.15)
			plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
			ax3.set_xlabel('time (Gyr)')
			ax1.set_ylabel('$\Delta\Omega$ (deg)')
			ax2.set_ylabel('$\Delta\omega$ (deg)')
			ax3.set_ylabel('$\Delta\\varpi$ (deg)')
			ax1.set_ylim(0, 360)
			ax1.set_xlim(min(time), max(time))
			ax1.text(0.5, 1.2, texto, fontsize = 14, ha = 'center', va = 'center', transform = ax1.transAxes)
			ax1.set_yticks([0, 100, 200, 300])
			ax1.set_xticks([0*10**9, 0.5*10**9, 1.0*10**9, 1.5*10**9, 2.0*10**9, 2.5*10**9, 3.0*10**9, 3.5*10**9, 4.0*10**9, 4.5*10**9])
			ax1.set_xticklabels(['$0$', '$0.5$', '$1.0$', '$1.5$', '$2.0$', '$2.5$', '$3.0$', '$3.5$', '$4.0$', '$4.5$'])
			ax1.minorticks_on()
			ax2.minorticks_on()
			ax3.minorticks_on()
			ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')	
			ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')	
			ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			f.savefig(imag3, format = 'pdf')
			
			# ---------------------------------------------------------------------------------------------------------------
			# Dados da média e desvio padrão dos ângulos relativos associados a N1 desde 4 Gyr para elaborar a Tabela 1
			
			d = [i for i, n in enumerate(time) if n >= 4*10**9][0]
			mDcapom500, dDcapom500, _ = function_media(Dcapom[d:])
			mDomega500, dDomega500, _ = function_media(Domega[d:])
			mDobar500, dDobar500, _ = function_media(Dobar[d:])

			outputfile = open(pasta_compilacao + 'tabela1.txt', 'a+')
			outputfile.write("%18s"%entrada + " & " + str("%5.1f"%mDcapom500) + " & " + str("%5.1f"%dDcapom500) + " & " + str("%5.1f"%mDomega500) + " & " + str("%5.1f"%dDomega500) + " & " + str("%5.1f"%mDobar500) + " & " + str("%5.1f"%dDobar500) + " \\\\" + "\n")
			outputfile.close()
			
			# ---------------------------------------------------------------------------------------------------------------
			# Plot de 'N1' e 'N2' vs. 't'

			f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
			ax1.plot(time, NDcapompeak1, marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4, label = '$N_1$')
			ax1.plot(time, NDcapompeak2, marker = '', linestyle = '-', color = 'lightcoral', linewidth = 0.4, label = '$N_2$')
			ax2.plot(time, NDomegapeak1, marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4)
			ax2.plot(time, NDomegapeak2, marker = '', linestyle = '-', color = 'lightcoral', linewidth = 0.4)
			ax3.plot(time, NDobarpeak1, marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4)
			ax3.plot(time, NDobarpeak2, marker = '', linestyle = '-', color = 'lightcoral', linewidth = 0.4)
			f.subplots_adjust(hspace = 0.15)
			plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
			ax3.set_xlabel('time (Gyr)')
			ax1.set_ylabel('$N$ in $\Delta\Omega$')
			ax2.set_ylabel('$N$ in $\Delta\omega$')
			ax3.set_ylabel('$N$ in $\Delta\\varpi$')
			ax3.set_ylim(0, 1)
			ax3.set_xlim(min(time), max(time))
			ax1.legend(loc = 2)
			ax1.text(0.5, 1.2, texto, fontsize = 14, ha = 'center', va = 'center', transform = ax1.transAxes)
			ax1.set_ylim(0, 1)
			ax1.set_yticks([0.0, 0.5, 1.0])
			ax1.set_xticks([0*10**9, 0.5*10**9, 1.0*10**9, 1.5*10**9, 2.0*10**9, 2.5*10**9, 3.0*10**9, 3.5*10**9, 4.0*10**9, 4.5*10**9])
			ax1.set_xticklabels(['$0$', '$0.5$', '$1.0$', '$1.5$', '$2.0$', '$2.5$', '$3.0$', '$3.5$', '$4.0$', '$4.5$'])
			ax1.minorticks_on()
			ax2.minorticks_on()
			ax3.minorticks_on()
			ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')	
			ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')	
			ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			f.savefig(imag2, format = 'pdf')
			
			# ---------------------------------------------------------------------------------------------------------------
			# Plot de 'N1' e 'N' da distribuição aleatória vs. 't'

			f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
			ax1.plot(time, NDcapompeak1, marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4, label = '$N_1$')
			ax1.plot(time, N, marker = '', linestyle = '-', color = 'darkorange', linewidth = 0.4, label = '$N_{random}$')
			ax2.plot(time, NDomegapeak1, marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4)
			ax2.plot(time, N, marker = '', linestyle = '-', color = 'darkorange', linewidth = 0.4)
			ax3.plot(time, NDobarpeak1, marker = '', linestyle = '-', color = 'steelblue', linewidth = 0.4)
			ax3.plot(time, N, marker = '', linestyle = '-', color = 'darkorange', linewidth = 0.4)
			f.subplots_adjust(hspace = 0.15)
			plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
			ax3.set_xlabel('time (Gyr)')
			ax1.set_ylabel('$N$ in $\Delta\Omega$')
			ax2.set_ylabel('$N$ in $\Delta\omega$')
			ax3.set_ylabel('$N$ in $\Delta\\varpi$')
			ax3.set_ylim(0, 1)
			ax3.set_xlim(min(time), max(time))
			ax1.legend(loc = 2)
			ax1.text(0.5, 1.2, texto, fontsize = 14, ha = 'center', va = 'center', transform = ax1.transAxes)
			ax1.set_ylim(0, 1)
			ax1.set_yticks([0.0, 0.5, 1.0])
			ax1.set_xticks([0*10**9, 0.5*10**9, 1.0*10**9, 1.5*10**9, 2.0*10**9, 2.5*10**9, 3.0*10**9, 3.5*10**9, 4.0*10**9, 4.5*10**9])
			ax1.set_xticklabels(['$0$', '$0.5$', '$1.0$', '$1.5$', '$2.0$', '$2.5$', '$3.0$', '$3.5$', '$4.0$', '$4.5$'])
			ax1.minorticks_on()
			ax2.minorticks_on()
			ax3.minorticks_on()
			ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')	
			ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')	
			ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			f.savefig(imag4, format = 'pdf')

			# ---------------------------------------------------------------------------------------------------------------
			# Plot de 's' vs. 't'
			
			f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
			ax1.semilogy(time, s_capom, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'black', markersize = 1.5)
			ax2.semilogy(time, s_omega, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'black', markersize = 1.5)
			ax3.semilogy(time, s_obar, marker = 'o', markeredgewidth = 0.0, linestyle = '', color = 'black', markersize = 1.5)
			f.subplots_adjust(hspace = 0.15)
			plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
			ax3.set_xlabel('time (Gyr)')
			ax1.set_ylabel('$s$ in $\Omega$')
			ax2.set_ylabel('$s$ in $\omega$')
			ax3.set_ylabel('$s$ in $\\varpi$')
			ax3.set_xlim(min(time), max(time))
			ax1.text(0.5, 1.2, texto, fontsize = 14, ha = 'center', va = 'center', transform = ax1.transAxes)
			ax1.set_xticks([0*10**9, 0.5*10**9, 1.0*10**9, 1.5*10**9, 2.0*10**9, 2.5*10**9, 3.0*10**9, 3.5*10**9, 4.0*10**9, 4.5*10**9])
			ax1.set_xticklabels(['$0$', '$0.5$', '$1.0$', '$1.5$', '$2.0$', '$2.5$', '$3.0$', '$3.5$', '$4.0$', '$4.5$'])
			ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')	
			ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'major')
			ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'major')
			ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')	
			ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True,  which = 'minor')
			ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True,  which = 'minor')
			locmaj = matplotlib.ticker.LogLocator(base = 10, numticks = 100) 
			ax1.yaxis.set_major_locator(locmaj)
			locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.2,0.4,0.6,0.8), numticks = 100)
			ax1.yaxis.set_minor_locator(locmin)
			ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
			f.savefig(imag5, format = 'pdf')

			# ---------------------------------------------------------------------------------------------------------------
			# Dados de s, para tudo o tempo de integração e para os últimos 500 Myr 
					
			smax_capom, smed_capom = function_est(s_capom)	
			smax_omega, smed_omega = function_est(s_omega)
			smax_obar, smed_obar = function_est(s_obar)
			d = [i for i, n in enumerate(time) if n >= 4*10**9][0]
			smax_capomd, smed_capomd = function_est(s_capom[d:])	
			smax_omegad, smed_omegad = function_est(s_omega[d:])
			smax_obard, smed_obard = function_est(s_obar[d:])

			est2 = open(pasta_compilacao + 's.txt', 'a+')
			est2.write("%18s"%entrada + " & " + "%6s"%str(smax_capom) + " & " + "%6s"%str(smed_capom) + " & " + "%6s"%str(smax_capomd) + " & " + "%6s"%str(smed_capomd) + " & " + "%6s"%str(smax_omega) + " & " + "%6s"%str(smed_omega) + " & " + "%6s"%str(smax_omegad) + " & " + "%6s"%str(smed_omegad) + " & " + "%6s"%str(smax_obar) + " & " + "%6s"%str(smed_obar) + " & " + "%6s"%str(smax_obard) + " & " + "%6s"%str(smed_obard) + "\\\\" + "\n")
			est2.close()		
	
		# -----------------------------------------------------------------------------------------------------------------------

		# -----------------------------------------------------------------------------------------------------------------------
		if executeMean:
		
			# Plot da evolução dos confinamentos
			
			mediaconfcapomdpmais = [x + y for x, y in zip(mediaconfcapom, desvioconfcapom)]
			mediaconfcapomdpmenos = [x - y for x, y in zip(mediaconfcapom, desvioconfcapom)]
			mediaconfomegadpmais = [x + y for x, y in zip(mediaconfomega, desvioconfomega)]
			mediaconfomegadpmenos = [x - y for x, y in zip(mediaconfomega, desvioconfomega)]
			mediaconfobardpmais = [x + y for x, y in zip(mediaconfobar, desvioconfobar)]
			mediaconfobardpmenos = [x - y for x, y in zip(mediaconfobar, desvioconfobar)]
			
			f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
			ax1.fill_between(time, mediaconfcapomdpmenos, mediaconfcapomdpmais, color = 'lightgray')
			ax1.plot(time, mediaconfcapom, marker = '', linestyle = '-', color = 'brown', linewidth = 0.4)
			ax2.fill_between(time, mediaconfomegadpmenos, mediaconfomegadpmais, color = 'lightgray')
			ax2.plot(time, mediaconfomega, marker = '', linestyle = '-', color = 'brown', linewidth = 0.4)
			ax3.fill_between(time, mediaconfobardpmenos, mediaconfobardpmais, color = 'lightgray')
			ax3.plot(time, mediaconfobar, marker = '', linestyle = '-', color = 'brown', linewidth = 0.4)
			f.subplots_adjust(hspace = 0.15)
			plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
			ax3.set_xlabel('time (Gyr)')
			ax1.set_ylabel('Conf. $\Omega$ (deg)')
			ax2.set_ylabel('Conf. $\omega$ (deg)')
			ax3.set_ylabel('Conf. $\\varpi$ (deg)')
			ax1.set_ylim(0, 360)
			ax1.set_xlim(min(time), max(time))
			ax1.text(0.5, 1.2, texto, fontsize = 14, ha = 'center', va = 'center', transform = ax1.transAxes)
			ax1.set_yticks([0, 100, 200, 300])
			ax1.set_xticks([0*10**9, 0.5*10**9, 1.0*10**9, 1.5*10**9, 2.0*10**9, 2.5*10**9, 3.0*10**9, 3.5*10**9, 4.0*10**9,
 4.5*10**9])
			ax1.set_xticklabels(['$0$', '$0.5$', '$1.0$', '$1.5$', '$2.0$', '$2.5$', '$3.0$', '$3.5$', '$4.0$', '$4.5$'])
			ax1.minorticks_on()
			ax2.minorticks_on()
			ax3.minorticks_on()
			ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')	
			ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
			ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')	
			ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
			f.savefig(imag1, format = 'pdf')
			
		# -----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

	t0 = time.time()

	# Dados de Entrada
	# ===============================================================================================================================
	#
	# numobs:			num. de objetos observados com a > 250 au e q > 40 au
	# numobsconfcapom:		num. de objetos observados confinados em long. do nodo
	# numobsconfomega:		num. de objetos observados confinados em argumento do periélio
	# numobsconfobar:		num. de objetos observados confinados em long. do periélio
	# numite: 	 		num. de iterações para calcular as médias e confinamentos em um certo tempo
	# aminobs, amaxobs:		limites min. e max. similar a dos objetos reais confinados em semieixo maior
	# qminobs, qmaxobs:		limites min. e max. similar a dos objetos reais confinados em distância ao periélio
	# incminobs, incmaxobs:	limites min. e max. similar a dos objetos reais confinados em inclinação
	# 				Ter em conta que os minimos valores fornecidos não sejam maiores a dos objetos observados!
	# executeHist:			é igual a True se deseja realizar analise dos histogramas
	# executeMean:			é igual a True se deseja realizar analise das média e confinamento
	# arquivo_entrada:		nome genérico dos arquivos de entrada (ex.: arquivo_i)
	# arquivo_invariante:		nome genérico dos arquivos invariantes que serão criados (ex.: arquivo_iv4_i)
	# narq:			num. total de arquivos de entrada (arquivo_i)
	# deltat:			intervalo de tempo entre os tempos considerados na análise, em anos
	# ninf_nND:			mínimo 'n' em nND
	# nsup_nND:			máximo 'n' em nND
	# rota_comum:			rota anterior à pasta de entrada (que contém os arquivos para a análise)
	# rota_nND:			rota ao arquivo nND
	# pasta_compilacao:		rota à pasta que será criada e contera os diversos plots
	# gm:				G vezes a massa central, no nosso caso o Sol. Aqui G = 1 e 1 massa solar é aprox. 2.96e-4
	#				unidades de massa
	# gmsum:			G*(M1+M2)
	# mpl:		 		lista que contem a massa dos 4 planetas gigantes e a do pl9 em massas solares. Nota: 
	#				para a9 = 700 au, massa9 = 3e-5 massas solares
	#				para a9 = 1500 au, massa9 = 5e-5 massas solares
	
	numobs = 11
	numobsconfcapom = 11
	numobsconfomega = 11
	numobsconfobar = 9
	numite = 100
	aminobs = 250
	amaxobs = 1300
	qminobs = 40
	qmaxobs = 90
	incminobs = 0
	incmaxobs = 30
	executeHist = True
	executeMean = True
	arquivo_entrada = 'arquivo_10Myr_' 	
	arquivo_invariante = 'arquivo_10Myr_iv4_'
	narq = 450
	deltat = 10**7	
	ninf_nND = 5 
	nsup_nND = 87541
	rota_comum = '/home/jessi/Documentos/Projeto-II/'
	rota_nND = '/home/jessi/Dropbox/Codigos-Projeto-GT/nND2.txt'
	pasta_compilacao = '/home/jessi/Dropbox/Compilacao-previa/'
	gm = 0.01720209895**2
	gmsum = gm
#	mpl = [0.000954786104, 0.000285836787, 4.37273165e-05, 5.17759138e-05, 3e-5]
	mpl = [0.000954786104, 0.000285836787, 4.37273165e-05, 5.17759138e-05, 5e-5]
	mpl = np.float32(np.multiply(gm, mpl))

	# -------------------------------------------------------------------------------------------------------------------------------
	# Nomes das corridas 
	
	entrada1 = 'a1500-30-0-60-5'
	entrada2 = 'a1500-10-0-100-5'
	entrada3 = 'a1500-30-0-100-5'
	entrada4 = 'a1500-60-0-100-5'
	entrada5 = 'a1500-30-0-200-5'
	entrada6 = 'a1500-30-0-300-5'
	entrada7 = 'a1500-30-90-60-5'
	entrada8 = 'a700-30-0-60-3'
	entrada9 = 'a700-10-0-100-3' 
	entrada10 = 'a700-30-0-100-3'
	entrada11 = 'a700-60-0-100-3'
	entrada12 = 'a700-30-0-200-3'
	entrada13 = 'a700-30-0-300-3'
	entrada14 = 'a700-30-90-60-3' 
	entrada15 = 'a2000-e0.8-i30-w100-O50'
	entrada16 = 'a2500-e0.8-i30-w30-O140'
	entrada17 = 'a3000-e0.7-i30-w0-O150'
	entrada18 = 'a3000-e0.7-i30-w90-O30'
	entrada19 = 'a3000-e0.9-i30-w310-O20'
	entrada20 = 'init-close-net/a3000-e0.7-i30-w0-O150'
	entrada21 = 'GT-a2000-30-350-400-5'
	entrada22 = 'a3000-e0.8-i30-w130-O20' 
	entrada23 = 'a3000-e0.9-i30-w130-O190'
	entrada24 = 'mig-a3000-e0.9-i30-w310-O20'
	entrada25 = 'mig-sempl9'
	entrada26 = 'mig-a3000-e0.7-i30-w0-O150'

	# -------------------------------------------------------------------------------------------------------------------------------
	# Texto segundo elementos orbitais iniciais

	texto1 = r"$a_9=1500$ au$,$ $q_9=60$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"	
	texto2 = r"$a_9=1500$ au$,$ $q_9=100$ au$,$ $i_9=10^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"	
	texto3 = r"$a_9=1500$ au$,$ $q_9=100$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto4 = r"$a_9=1500$ au$,$ $q_9=100$ au$,$ $i_9=60^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto5 = r"$a_9=1500$ au$,$ $q_9=200$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto6 = r"$a_9=1500$ au$,$ $q_9=300$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto7 = r"$a_9=1500$ au$,$ $q_9=60$ au$,$ $i_9=30^\circ,$ $\omega_9=90^\circ,$ $\Omega_9=150^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto8 = r"$a_9=700$ au$,$ $q_9=60$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=10$ $M_{\oplus}$"
	texto9 = r"$a_9=700$ au$,$ $q_9=100$ au$,$ $i_9=10^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=10$ $M_{\oplus}$"
	texto10 = r"$a_9=700$ au$,$ $q_9=100$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=10$ $M_{\oplus}$"
	texto11 = r"$a_9=700$ au$,$ $q_9=100$ au$,$ $i_9=60^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=10$ $M_{\oplus}$"
	texto12 = r"$a_9=700$ au$,$ $q_9=200$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=10$ $M_{\oplus}$"
	texto13 = r"$a_9=700$ au$,$ $q_9=300$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=0^\circ,$" "\n" r"$m_9=10$ $M_{\oplus}$"
	texto14 = r"$a_9=700$ au$,$ $q_9=60$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=90^\circ,$" "\n" r"$m_9=10$ $M_{\oplus}$"
	texto15 = r"$a_9=2000$ au$,$ $q_9=400$ au$,$ $i_9=30^\circ,$ $\omega_9=100^\circ,$ $\Omega_9=50^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto16 = r"$a_9=2500$ au$,$ $q_9=500$ au$,$ $i_9=30^\circ,$ $\omega_9=30^\circ,$ $\Omega_9=140^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto17 = r"$a_9=3000$ au$,$ $q_9=900$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=150^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto18 = r"$a_9=3000$ au$,$ $q_9=900$ au$,$ $i_9=30^\circ,$ $\omega_9=90^\circ,$ $\Omega_9=300^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto19 = r"$a_9=3000$ au$,$ $q_9=300$ au$,$ $i_9=30^\circ,$ $\omega_9=310^\circ,$ $\Omega_9=20^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto20 = r"$a_9=3000$ au$,$ $q_9=900$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=150^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto21 = r"$a_9=2000$ au$,$ $q_9=400$ au$,$ $i_9=30^\circ,$ $\omega_9=350^\circ,$ $\Omega_9=160^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto22 = r"$a_9=3000$ au$,$ $q_9=600$ au$,$ $i_9=30^\circ,$ $\omega_9=130^\circ,$ $\Omega_9=20^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto23 = r"$a_9=3000$ au$,$ $q_9=300$ au$,$ $i_9=30^\circ,$ $\omega_9=130^\circ,$ $\Omega_9=190^\circ,$" "\n" r"$m_9=16.5$ $M_{\oplus}$"
	texto24 = r"$a_9=3000$ au$,$ $q_9=600$ au$,$ $i_9=30^\circ,$ $\omega_9=130^\circ,$ $\Omega_9=20^\circ,$" "\n" "$m_9=16.5$ $M_{\oplus}$"
	texto25 = "Migration without planet 9"
	texto26 = r"$a_9=3000$ au$,$ $q_9=900$ au$,$ $i_9=30^\circ,$ $\omega_9=0^\circ,$ $\Omega_9=150^\circ,$" "\n" "$m_9=16.5$ $M_{\oplus}$"
	
	lista_entrada = [entrada1, entrada2, entrada3, entrada4, entrada5, entrada6, entrada7, entrada8, entrada9, entrada10, entrada11, entrada12, entrada13, entrada14, entrada15, entrada16, entrada17, entrada18, entrada19, entrada20, entrada21, entrada22, entrada23, entrada24, entrada25, entrada26]
	
	lista_texto = [texto1, texto2, texto3, texto4, texto5, texto6, texto7, texto8, texto9, texto10, texto11, texto12, texto13, texto14, texto15, texto16, texto17, texto18, texto19, texto20, texto21, texto22, texto23, texto24, texto25, texto26]
	
	p = Pool()
	p.map(AnaliseProjeto, zip(lista_entrada[3:5], lista_texto[3:5]))

	t1 = time.time()
	seconds = t1 - t0
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print ("h:min:s =", "%d:%02d:%02d" % (h, m, s))
