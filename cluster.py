'''
- Script que gera condições iniciais de estrelas do cluster.
- Analisa os resultados a partir do arquivo de saída e cria novas condições iniciais 
a partir dele, para as próximas simulações que adicionam J, S e partículas testes em
torno a uma estrela similar ao Sol.
'''

import math
import random
import matplotlib.pyplot as plt
import numpy as np
import time
import os
from matplotlib import rc
# Criando módulos a partir de subrotinas Fortran do Swift
os.system('f2py -c -m orbel_xv2el orbel_xv2el.f')
os.system('f2py -c -m orbel_el2xv orbel_el2xv.f orbel_scget.f orbel_ehybrid.f orbel_fhybrid.f orbel_schget.f orbel_zget.f orbel_flon.f orbel_fget.f orbel_esolmd.f orbel_eget.f orbel_ehie.f')
import orbel_xv2el
import orbel_el2xv

# Função que randomiza o parâmetro "r" com probabilidade uniforme
def randomize(r):

	X2 = random.random()
	X3 = random.random()
	z = (1 - 2*X2)*r
	x = math.sqrt(r**2 - z**2)*math.cos(2*math.pi*X3)
	y = math.sqrt(r**2 - z**2)*math.sin(2*math.pi*X3)
	
	return x, y, z
	
def g(q):

	return q**2*(1 - q**2)**(7/2)

# Função que converte posições e velocidades em relação ao baricentro
def coord_h2b(nbod, mpl, xh, yh, zh, vxh, vyh, vzh):

	xb = []
	yb = []
	zb = []
	vxb = []
	vyb = []
	vzb = []

	msys = 1e-20
	xtmp = 0.0
	ytmp = 0.0
	ztmp = 0.0
	vxtmp = 0.0
	vytmp = 0.0
	vztmp = 0.0

	for n in range(nbod):

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

	for n in range(nbod):

		xb.append(xh[n] + xbsun)
		yb.append(yh[n] + ybsun)
		zb.append(zh[n] + zbsun)
		vxb.append(vxh[n] + vxbsun)
		vyb.append(vyh[n] + vybsun)
		vzb.append(vzh[n] + vzbsun)

	return xb, yb, zb, vxb, vyb, vzb

def CI_cluster(modelo, rho_0, c):

	os.system('mkdir ' + modelo)

	# Massa do cluster (gás mais estrelas), em massas solares
	# -----------------------------------------------------------------------------------------------------------

	# O fator 210000 é de conversão de au à parsec
	M = 4*math.pi/3*rho_0*(c/210000)**3 
	print('Massa do cluster, em massas solares:', M)

	# Coleção de massas de estrelas
	# -----------------------------------------------------------------------------------------------------------

	# Gera um float aleatório uniformemente na faixa semiaberta [0.0, 1.0)
	epsilon = random.random() 
	# M_j: Massa da estrela j, em massas solares
	M_j = [0.01 + (0.19*epsilon**1.55 + 0.05*epsilon**0.6)/(1 - epsilon)**0.58]
	# M_estrelas: Massa acumulada de estrelas dentro de um raio, em massas solares
	M_estrelas = [M_j[-1]]
	# fracao: Fração de massa do cluster contida nas estrelas
	fracao = 0.66*M

	while M_estrelas[-1] <= fracao:
	
		epsilon = random.random()
		M_j.append(0.01 + (0.19*epsilon**1.55 + 0.05*epsilon**0.6)/(1 - epsilon)**0.58)
		M_estrelas.append(M_estrelas[-1] + M_j[-1])
		
	M_j.pop()
	M_estrelas.pop()

	print('Massa contida nas estrelas, em massas solares:', M_estrelas[-1])
	print('Fração de massa no cluster contida nas estrelas:', M_estrelas[-1]/M)
	print('Massas estelares minima, média e máxima, em massas solares:', min(M_j), ',', np.mean(M_j), ',', max(M_j))
	print('Número de estrelas:', len(M_j))

	# Energia do cluster em unidades tal que G = 1 (longitudes em au, velocidades em au/dia e uma massa solar é 
	# aproximadamente 2.96e-4 unidades de massa)
	# -----------------------------------------------------------------------------------------------------------

	G = 1
	energy = -3*math.pi*G*(M*2.96*10**(-4))**2/(64*c)

	# Escolhemos um referencial onde a distância da origem ao centro de massa das estrelas seja pequena (< 100 au 
	# já é satisfatório). Essas posições e velocidades as definimos em relação ao baricentro.
	# -----------------------------------------------------------------------------------------------------------

	r_cm = 101
	while r_cm > 50:

		# Definindo matriz que salva dados
		m = []
		x_cm = 0
		y_cm = 0
		z_cm = 0

		for i in range(len(M_j)):

			# Tomando G = 1, M = 1, c = 1 por conveniência

			# X1: massa acumulada do cluster (estrelas mais gás) até um certo raio dado em função da massa 
			# acumulada de estrelas, assumindo uma razão de estrelas à gás constante através do cluster, 
			# entre a massa do cluster
			X1 = M_estrelas[i]*(1 + (1 - 0.66)/0.66)/M
			# Distância da estrela ao centro do cluster
			r = (X1**(-2/3) - 1)**(-1/2)
			x, y, z = randomize(r)
		
			# Geramos dois números aleatórios X4 e X5; se 0.1*X5 < g(X4), adotamos q = X4; caso contrário, 
			# um novo par de números aleatórios é tentado
			X4 = random.random()
			X5 = random.random()
			while 0.1*X5 >= g(X4):
				X4 = random.random()
				X5 = random.random()
			q = X4
			
			# Velocidade de escape
			ve = 2**(1/2)*(1 + r**2)**(-1/4)
			# Módulo da velocidade
			v = q*ve
			vx, vy, vz = randomize(v)

			# Ajustando as escalas, enquanto mantém G = 1. Longitudes devem ser multiplicadas por 
			# 3*math.pi*M**2/(64*abs(energy)), e velocidades por 64*abs(energy)**(1/2)/(3*math.pi*M**(1/2))
			x = x*3*math.pi*(M*2.96*10**(-4))**2/(64*abs(energy))
			y = y*3*math.pi*(M*2.96*10**(-4))**2/(64*abs(energy))
			z = z*3*math.pi*(M*2.96*10**(-4))**2/(64*abs(energy))
			vx = vx*64*abs(energy)**(1/2)/(3*math.pi*(M*2.96*10**(-4))**(1/2))
			vy = vy*64*abs(energy)**(1/2)/(3*math.pi*(M*2.96*10**(-4))**(1/2))
			vz = vz*64*abs(energy)**(1/2)/(3*math.pi*(M*2.96*10**(-4))**(1/2))
			r = r*3*math.pi*(M*2.96*10**(-4))**2/(64*abs(energy))
		
			m.append([X1*M, r, x, y, z, vx, vy, vz])
			
			x_cm = x_cm + M_j[i]*x
			y_cm = y_cm + M_j[i]*y
			z_cm = z_cm + M_j[i]*z
			
		r_cm = math.sqrt( (x_cm/M_estrelas[-1])**2 + (y_cm/M_estrelas[-1])**2 + (z_cm/M_estrelas[-1])**2 )
		
	print('Distância da origem ao centro de massa, em au:', r_cm)

	m = np.asarray(m)
	xb, yb, zb, vxb, vyb, vzb = coord_h2b(len(M_j), M_j, m[:, 2], m[:, 3], m[:, 4], m[:, 5], m[:, 6], m[:, 7])

	# Escrevendo arquivo 'big.in' em formatos Cartesian e Asteroidal
	# -----------------------------------------------------------------------------------------------------------

	data = open(modelo + '/big.in_cartesian', 'w')
	data2 = open(modelo + '/big.in_asteroidal', 'w')
	string1 = ")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n"
	string2 = ") Lines beginning with `)' are ignored.\n"
	string3 = ")---------------------------------------------------------------------\n"
	string4 = " style (Cartesian, Asteroidal, Cometary) = Cartesian\n"
	string5 = " style (Cartesian, Asteroidal, Cometary) = Ast\n"
	string6 = " epoch (in days) = 0\n"
	string7 = ")---------------------------------------------------------------------\n"
	data.write(string1 + string2 + string3 + string4 + string6 + string7)
	data2.write(string1 + string2 + string3 + string5 + string6 + string7)

	for i in range(len(m)):

		# m: massa do corpo em massas solares
		# r: distância máxima que constitui um encontro próximo em raios de Hill
		# d: densidade em g/cm**3
		string = "{:5s}  {:25s}  r=1  d=1.408\n".format('s' + str(i + 1), 'm=' + str(M_j[i]))
		data.write(string)
		data2.write(string)
		# posições e velocidades em au e au/dia
		data.write("{:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  0  0  0\n".format(xb[i], yb[i], zb[i], vxb[i], vyb[i], vzb[i]))
		# gmsum em unidades tal que G = 1
		gmsum = (M_estrelas[-1] + M_j[i])*2.96*10**-4
		# Este método retorna os ângulos em radianos
		_, a, e, inc, capom, omega, capm = orbel_xv2el.orbel_xv2el(xb[i], yb[i], zb[i], vxb[i], vyb[i], vzb[i], gmsum)
		data2.write("{:25.4f}  {:2.9f}  {:8.4f}  {:8.4f}  {:8.4f}  {:8.4f}  0  0  0\n".format(a, e, inc*180/math.pi, omega*180/math.pi, capom*180/math.pi, capm*180/math.pi))

	data.close()
	data2.close()

	# Visualização de dados
	# ----------------------------------------------------------------------------------------------------------

	rc('font', **{'weight': 'light', 'family': 'serif', 'serif': ['Times'], 'size': 24})
	rc('text', usetex = True)
	rc('axes', titlesize = 24)

	fig = plt.figure(figsize = (24, 16))
	fig.subplots_adjust()
	ax1 = plt.subplot2grid(shape = (2, 6), loc = (0, 0), colspan = 2)
	ax2 = plt.subplot2grid((2, 6), (0, 2), colspan = 2)
	ax3 = plt.subplot2grid((2, 6), (0, 4), colspan = 2)
	ax4 = plt.subplot2grid((2, 6), (1, 1), colspan = 2, projection = '3d')
	ax5 = plt.subplot2grid((2, 6), (1, 3), colspan = 2, projection = '3d')
	ax1.hist(M_j, 50)
	ax1.set_xlabel('Mass (solar mass)')
	ax1.set_ylabel('N')
	ax1.set_title('Mass distribution')
	ax2.hist(m[:, 1], 50)
	ax2.set_xlabel('r (au)')
	ax2.set_ylabel('N')
	ax2.set_title('Radial distance distribution')
	ax3.plot(m[:, 1], m[:, 0], 'o', markersize = 2)
	ax3.set_xlabel('r (au)')
	ax3.set_ylabel('M(r) (solar mass)')
	ax3.set_title('Accumulated mass vs. radial distance')
	# Distâncias radiais no espaço em relação ao baricentro	
	ax4.plot3D(xb, yb, zb, 'o', markersize = 3, alpha = 0.5)
	ax4.set_xlabel('x (au)', labelpad = 20)
	ax4.set_ylabel('y (au)', labelpad = 30)
	ax4.set_zlabel('z (au)', labelpad = 40)
	ax4.set_title('Radial distances')
	# Módulos das velocidades no espaço em relação ao baricentro
	ax5.plot3D(vxb, vyb, vzb, 'o', markersize = 2.5, alpha = 0.5)
	ax5.set_xlabel('vx (au)', labelpad = 20)
	ax5.set_ylabel('vy (au)', labelpad = 30)
	ax5.set_zlabel('vz (au)', labelpad = 30)
	ax5.set_title('Velocity modules')
	ax1.minorticks_on()
	ax2.minorticks_on()
	ax3.minorticks_on()
	ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
	ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
	ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
	fig.tight_layout()
	plt.savefig(modelo + '/' + modelo + '.pdf')
	plt.close()
	
def Analise_cluster(modelo):

	# Extraindo as massas e identificadores das estrelas alvos
	# ---------------------------------------------------------------------------------------------------

	data = open(modelo + '/big.in_cartesian')
	lines = data.read().splitlines()
	data.close()

	stars = []
	masses = []
	for i in range(6, len(lines), 2):
	
		m = float(lines[i].split()[1].replace('m=', ''))
		
		if m >= 0.9 and m <= 1.1:
		
			stars.append(lines[i].split()[0])
			masses.append(m)
			
	print('Número de estrelas com massas entre 0.9 - 1.1 massas solares:', len(stars))
	
	# Visualizando dados das estrelas alvos
	# ---------------------------------------------------------------------------------------------------
	
	n = np.genfromtxt(modelo + '/planet.out', dtype = str, usecols = (0), unpack = True)
	xb, yb, zb = np.loadtxt(modelo + '/planet.out', usecols = (2, 3, 4), unpack = True)
	rb = (xb**2 + yb**2 + zb**2)**(1/2)
	rcmin = min(rb)
	rcmean = np.mean(rb)
	rcmax = max(rb)

	t = np.arange(0, 10**7 + 10**3, 10**3)	
	
	rini = []
	for k in stars:

		x = [j for i, j in zip(n, xb) if i == k]
		y = [j for i, j in zip(n, yb) if i == k]
		z = [j for i, j in zip(n, zb) if i == k]
		r = [j for i, j in zip(n, rb) if i == k]
		
		rini.append(r[-1])
		rmin = min(r)
		rmean = np.mean(r)
		rmax = max(r)
		
		fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (12, 4))	
		ax1.plot(x, y, marker = 'o', markersize = 0.25)
		ax1.set_xlabel('x (au)')
		ax1.set_ylabel('y (au)')
		ax2.plot(x, z, marker = 'o', markersize = 0.25)
		ax2.set_xlabel('x (au)')
		ax2.set_ylabel('z (au)')
		ax3.plot(y, z, marker = 'o', markersize = 0.25)
		ax3.set_xlabel('y (au)')
		ax3.set_ylabel('z (au)')
		fig.subplots_adjust(top = 0.8, wspace = 0.6)
		ax1.minorticks_on()
		ax2.minorticks_on()
		ax3.minorticks_on()
		ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
		ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
		ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True,  which = 'both')
		fig.suptitle('mass = ' + str('%.2f'%masses[stars.index(k)]) + ' solar masses\n$r_{{ini}}=${:.2e} au, $r_{{min}}=${:.2e} au, $r_{{max}}=${:.2e} au, $r_{{mean}}=${:.2e} au'.format(rini[-1], rmin, rmean, rmax))
		plt.savefig(modelo + '/' + str(k) + '.png')
		plt.close()
	
		fig = plt.figure()
		ax = plt.axes(projection = '3d')
		ax.plot3D(x, y, z, 'gray', markersize = 2)
		ax.scatter3D(x, y, z, c = t, norm = plt.Normalize(vmin = -max(t), vmax = max(t)), cmap = 'Greens', s = 0.05)
		ax.set_xlabel('x (au)') 
		ax.set_ylabel('y (au)')
		ax.set_zlabel('z (au)')
		ax.set_title('mass = ' + str('%.2f'%masses[stars.index(k)]) + ' solar masses\n$r_{{ini}}=${:.2e} au\n$r_{{min}}=${:.2e} au, $r_{{max}}=${:.2e} au, $r_{{mean}}=${:.2e} au'.format(rini[-1], rmin, rmean, rmax))
		plt.savefig(modelo + '/' + k + '_3d.png')
		plt.close()
		
	# Escolhemos 3 estrelas para testes com massa similar a do Sol, com distâncias iniciais próximas ao
	# centro, meio e borda externa do cluster	
	# ---------------------------------------------------------------------------------------------------
	
	dc = [rcmin, rcmean, rcmax]
	rini_testes = []
	stars_testes = []
	
	for i in range(3):
	
		rini_testes.append(min(rini, key = lambda j: abs(j - dc[i])))
		stars_testes.append(stars[rini.index(rini_testes[-1])])

	print('\nEstrelas testes:', stars_testes)
	print('Distâncias radiais iniciais das estrelas testes:', rini_testes)
	print('Extensão radial mínima, meia e máxima do cluster:', dc)
	print('Distâncias radiais iniciais de todas as estrelas alvos:', rini)
		
	# Criando novas condições iniciais: cluster e J, S e partículas testes em torno a uma das 
	# stars_testes. Isto para cada uma das 3 stars_testes
	# ---------------------------------------------------------------------------------------------------
	
	# CI de Júpiter e Saturno em RMM 3:2, em relação à estrela (em au e au/dia)
	jup = [3.614955269300000E+00, 4.211299632400000E+00, 1.814977848300000E-03, -5.671044120300000E-03, 4.657968368600000E-03, -2.410766223100000E-06]
	sat = [-7.379921071000000E+00, -2.766075235900000E+00, 5.700863805400000E-03, 2.214674421500000E-03, -5.456490811600000E-03, -5.740283918800000E-06]

	string1 = ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n"
	string2 = ") Lines beginning with `)' are ignored.\n"
	string3 = ")---------------------------------------------------------------------\n"
	string4 = " style (Cartesian, Asteroidal, Cometary) = Cartesian\n"
	string5 = ")---------------------------------------------------------------------\n"

	for i in range(3):

		os.system('cp ' + modelo + '/big.in_cartesian ' + modelo + '/big.in_cartesian_' + str(stars_testes[i]))		
		data = open(modelo + '/big.in_cartesian', 'r')
		flag = 0
		for linha in data:
			if flag == 1:	
				xa = float(linha.split()[0])
				ya = float(linha.split()[1])
				za = float(linha.split()[2])
				vxa = float(linha.split()[3])
				vya = float(linha.split()[4])
				vza = float(linha.split()[5])
				break
			if linha.split()[0] == stars_testes[i]:
				flag = 1
		data.close()
				
		data = open(modelo + '/big.in_cartesian_' + str(stars_testes[i]), 'a')		
		data.write("jup    m=0.000954786104   d=1.3\n")
		data.write("{:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  0  0  0\n".format(jup[0] + xa, jup[1] + ya, jup[2] + za, jup[3] + vxa, jup[4] + vya, jup[5] + vza))
		data.write("sat    m=0.000285836787   d=0.7\n")
		data.write("{:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  0  0  0\n".format(sat[0] + xa, sat[1] + ya, sat[2] + za, sat[3] + vxa, sat[4] + vya, sat[5] + vza))	
		data.close()	
		
		data = open(modelo + '/small.in', 'w')
		data.write(string1 + string2 + string3 + string4 + string5)
		step = (11 - 4)/(10000 - 1)
		a = 4
		k = 1
		gm = masses[stars.index(stars_testes[i])]*2.96*10**-4
		print('massa = ',masses[stars.index(stars_testes[i])])
		while a <= 11:
			e = 0
			inc = 0
			capom = random.uniform(0, 360)*math.pi/180
			omega = random.uniform(0, 360)*math.pi/180
			capm = random.uniform(0, 360)*math.pi/180
			x, y, z, vx, vy, vz = orbel_el2xv.orbel_el2xv(gm, -1, a, e, inc, capom, omega, capm)
			data.write('p' + str(k) + '\n')
			data.write("{:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  0  0  0\n".format(x + xa, y + ya, z + za, vx + vxa, vy + vya, vz + vza))
			a = a + step
			k = k + 1
		data.close()

# Script
# -----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

	t0 = time.time() 

	modelo = 'Modelo_1'
	# rho_0: Densidade central do cluster, em massas solares/parsec**3
	rho_0 = 3*10**4
	# c: Raio de Plummer, em au
	c = 30000

	#CI_cluster(modelo, rho_0, c)
	Analise_cluster(modelo)

	t1 = time.time()
	seconds = t1 - t0
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print ("h:min:s =", "%d:%02d:%02d" % (h, m, s))
