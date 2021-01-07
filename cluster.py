'''
- Script que gera condições iniciais de estrelas do cluster.
- Visualiza resultados do arquivo de saída da simulação e cria novos conjuntos de condições iniciais, cada um
adicionando J, S e partículas testes em torno a uma das estrelas com massas similares a do Sol.
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

	# Extraindo os dados iniciais e particularmente os dados de estrelas com massas similares a do Sol
	# ---------------------------------------------------------------------------------------------------

	data = open(modelo + '/big.in_cartesian')
	lines = data.read().splitlines()
	data.close()

	s = []
	xb0 = []
	yb0 = []
	zb0 = []
	vxb0 = []
	vyb0 = []
	vzb0 = []
	ss = []
	ms = []	
	
	for i in range(6, len(lines)):
		if i % 2 == 0:
			s.append(lines[i].split()[0])
			m = float(lines[i].split()[1].replace('m=', ''))
			if m >= 0.9 and m <= 1.1:
				ss.append(s[-1])
				ms.append(m)
		else:
			xb0.append(float(lines[i].split()[0]))
			yb0.append(float(lines[i].split()[1]))
			zb0.append(float(lines[i].split()[2]))
			vxb0.append(float(lines[i].split()[3]))
			vyb0.append(float(lines[i].split()[4]))
			vzb0.append(float(lines[i].split()[5]))		
			
	rb0 = [math.sqrt(i*i + j*j + k*k) for i, j, k in zip(xb0, yb0, zb0)]	
	
	rs = []
	for i in ss:
		rs.append(math.sqrt(xb0[s.index(i)]**2 + yb0[s.index(i)]**2 + zb0[s.index(i)]**2))

	# Visualizando trajetórias de estrelas com massas similares a do Sol, durante o tempo de integração
	# ---------------------------------------------------------------------------------------------------
	
	n = np.genfromtxt(modelo + '/planet.out', dtype = str, usecols = (0), unpack = True)
	xb, yb, zb = np.loadtxt(modelo + '/planet.out', usecols = (2, 3, 4), unpack = True)
	rb = (xb**2 + yb**2 + zb**2)**(1/2)
	
	t = np.arange(0, 10**7 + 10**3, 10**3)	
	
	for k in ss:

		x = [j for i, j in zip(n, xb) if i == k]
		y = [j for i, j in zip(n, yb) if i == k]
		z = [j for i, j in zip(n, zb) if i == k]
		r = [j for i, j in zip(n, rb) if i == k]
		
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
		fig.suptitle('mass = ' + str('%.2f'%ms[ss.index(k)]) + ' solar masses\n$r_{{ini}}=${:.2e} au, $r_{{min}}=${:.2e} au, $r_{{max}}=${:.2e} au, $r_{{mean}}=${:.2e} au'.format(r[0], min(r), np.mean(r), max(r)))
		plt.savefig(modelo + '/' + str(k) + '.png')
		plt.close()
	
		fig = plt.figure()
		ax = plt.axes(projection = '3d')
		ax.plot3D(x, y, z, 'gray', markersize = 2)
		ax.scatter3D(x, y, z, c = t, norm = plt.Normalize(vmin = -max(t), vmax = max(t)), cmap = 'Greens', s = 0.05)
		ax.set_xlabel('x (au)') 
		ax.set_ylabel('y (au)')
		ax.set_zlabel('z (au)')
		ax.set_title('mass = ' + str('%.2f'%ms[ss.index(k)]) + ' solar masses\n$r_{{ini}}=${:.2e} au\n$r_{{min}}=${:.2e} au, $r_{{max}}=${:.2e} au, $r_{{mean}}=${:.2e} au'.format(r[0], min(r), np.mean(r), max(r)))
		plt.savefig(modelo + '/' + k + '_3d.png')
		plt.close()
	
	# Pegamos três estrelas dentre aquelas com massas similares a do Sol que tenham distâncias iniciais 
	# próximas ao centro, meio e borda externa do cluster, param serem nossos alvos	
	# ---------------------------------------------------------------------------------------------------	
	dc = [min(rb0), np.mean(rb0), max(rb0)]
	
	ra = []
	sa = []
	
	for i in range(3):
		ra.append(min(rs, key = lambda j: abs(j - dc[i])))
		sa.append(ss[rs.index(ra[-1])])

	# Adicionando algumas estrelas com massas similares a do Sol cujas trajetórias tem formas 
	# distorsionadas de rosetas
	if modelo == "Modelo_1":
		sa.append('s117')
		ra.append(rs[ss.index(sa[-1])])
	if modelo == "Modelo_2":
		sa.append('s545')
		ra.append(rs[ss.index(sa[-1])])

	print('\nNúmero de estrelas no cluster:', len(s))	
	print('\nNúmero de estrelas com massas entre 0.9 - 1.1 massas solares:', len(ss))
	print('--Nomes:', ss)
	print('--Massas:', ms)
	print('--Distâncias radiais iniciais:', rs)
	print('\nDistâncias radiais iniciais mínima, meia e máxima dentre as estrelas do cluster:', dc)	
	print('\nEstrelas alvos para testes:')		
	print('--Nomes:', sa)
	print('--Distâncias radiais iniciais:', ra)	
		
	# Criando novas CI: cluster e J, S e partículas testes em torno a uma das estrelas alvos. Isto para 
	# cada uma das estrelas alvos
	# ---------------------------------------------------------------------------------------------------
	
	# Posições e velocidades de J e S em RMM 3:2, em relação à estrela alvo (em au e au/dia)
	jup = [3.614955269300000E+00, 4.211299632400000E+00, 1.814977848300000E-03, -5.671044120300000E-03, 4.657968368600000E-03, -2.410766223100000E-06]
	sat = [-7.379921071000000E+00, -2.766075235900000E+00, 5.700863805400000E-03, 2.214674421500000E-03, -5.456490811600000E-03, -5.740283918800000E-06]
	print('Jupiter =', orbel_xv2el.orbel_xv2el(jup[0], jup[1], jup[2], jup[3], jup[4], jup[5], 2.96*10**-4))
	print('Saturn =', orbel_xv2el.orbel_xv2el(sat[0], sat[1], sat[2], sat[3], sat[4], sat[5], 2.96*10**-4))
	
	string1 = ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n"
	string2 = ") Lines beginning with `)' are ignored.\n"
	string3 = ")---------------------------------------------------------------------\n"
	string4 = " style (Cartesian, Asteroidal, Cometary) = Cartesian\n"
	string5 = ")---------------------------------------------------------------------\n"

	for i in range(len(sa)):

		os.system('cp ' + modelo + '/big.in_cartesian ' + modelo + '/big.in_' + str(sa[i]))	
		
		# Posições e velocidades da estrela alvo em relação ao baricentro
		xba = xb0[s.index(sa[i])]
		yba = yb0[s.index(sa[i])]
		zba = zb0[s.index(sa[i])]
		vxba = vxb0[s.index(sa[i])]
		vyba = vyb0[s.index(sa[i])]
		vzba = vzb0[s.index(sa[i])]
				
		data = open(modelo + '/big.in_' + str(sa[i]), 'a')		
		data.write("jup    m=0.000954786104   d=1.3\n")
		# Posições e velocidades de J em relação ao baricentro
		xbJ = jup[0] + xba
		ybJ = jup[1] + yba
		zbJ = jup[2] + zba
		vxbJ = jup[3] + vxba
		vybJ = jup[4] + vyba
		vzbJ = jup[5] + vzba
		data.write("{:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  0  0  0\n".format(xbJ, ybJ, zbJ, vxbJ, vybJ, vzbJ))
		data.write("sat    m=0.000285836787   d=0.7\n")
		# Posições e velocidades de S em relação ao baricentro
		xbS = sat[0] + xba
		ybS = sat[1] + yba
		zbS = sat[2] + zba
		vxbS = sat[3] + vxba
		vybS = sat[4] + vyba
		vzbS = sat[5] + vzba
		data.write("{:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  0  0  0\n".format(xbS, ybS, zbS, vxbS, vybS, vzbS))	
		data.close()	
		
		data = open(modelo + '/small.in_' + str(sa[i]), 'w')
		data.write(string1 + string2 + string3 + string4 + string5)
		step = (11 - 4)/(10000 - 1)
		a = 4
		k = 1
		gm = ms[ss.index(sa[i])]*2.96*10**-4
		xbp = []
		ybp = []
		zbp = []
		vxbp = []
		vybp = []
		vzbp = []
		while a < 11 + step:
			# Elementos orbitais das partículas testes em relação à estrela alvo
			e = 0
			inc = 0
			capom = random.uniform(0, 360)*math.pi/180
			omega = random.uniform(0, 360)*math.pi/180
			capm = random.uniform(0, 360)*math.pi/180
			xsp, ysp, zsp, vxsp, vysp, vzsp = orbel_el2xv.orbel_el2xv(gm, -1, a, e, inc, capom, omega, capm)
			# Posições e velocidades das partículas testes em relação ao baricentro
			xbp.append(xsp + xba)
			ybp.append(ysp + yba)
			zbp.append(zsp + zba)
			vxbp.append(vxsp + vxba)
			vybp.append(vysp + vyba)
			vzbp.append(vzsp + vzba)			
			data.write('p' + str(k) + '\n')
			data.write("{:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  {:25.20f}  0  0  0\n".format(xbp[-1], ybp[-1], zbp[-1], vxbp[-1], vybp[-1], vzbp[-1]))
			a = a + step
			k = k + 1
		data.close()
		
		fig = plt.figure()
		ax = plt.axes(projection = '3d')
		ax.scatter3D(xb0, yb0, zb0, c = 'lightblue', s = 20, alpha = 0.5)
		ax.scatter3D(xba, yba, zba, c = 'brown', s = 20)
		ax.scatter3D([xbJ, xbS], [ybJ, ybS], [zbJ, zbS], c = 'lightcoral', s = 20)
		ax.scatter3D(xbp, ybp, zbp, c = 'gray', s = 0.5, alpha = 0.5)
		ax.set_xlabel('x (au)') 
		ax.set_ylabel('y (au)')
		ax.set_zlabel('z (au)')
		ax.set_xlim(xba - 15, xba + 15)
		ax.set_ylim(yba - 15, yba + 15)
		ax.set_zlim(zba - 15, zba + 15)
		ax.set_title(str(sa[i]))
		plt.show()
	
# Script
# -----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

	t0 = time.time() 

	# modelo: Nome da pasta que conterá os arquivos criados e o arquivo planet.out
	modelo = 'Modelo_2'
	# rho_0: Densidade central do cluster, em massas solares/parsec**3
	rho_0 = 3*10**4
	# c: Raio de Plummer, em au
	c = 30000

	# Cria CI de somente estrelas do cluster
	#CI_cluster(modelo, rho_0, c)
	
	# Visualiza resultados da simulação e cria conjuntos novos de CI
	Analise_cluster(modelo)

	t1 = time.time()
	seconds = t1 - t0
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print ("h:min:s =", "%d:%02d:%02d" % (h, m, s))
