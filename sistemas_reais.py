'''
- Script que cria os três arquivos de entrada do Swift ('start.in', 'param.in', 'pl.in') 
para o integrador SyMBA. SyMBA não pede o arquivo 'tp.in', pois todos os corpos do 
sistema são considerados massivos.
- Simula o modelo.
- Plota resultados.
- Cria grade de condições iniciais para Troianos.
- Gera diretórios com arquivos de entrada para as próximas integrações com rmvs.
'''

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import time
# Criando módulo a partir de subrotina Fortran do Swift
os.system('f2py -c -m orbel_xv2el orbel_xv2el.f')
import orbel_el2xv

def detialpha(e):

	if e >= 0 and e < 1:
		# Órbita elíptica
		ialpha = -1
	elif e > 1:
		# Órbita hiperbólica
		ialpha = 1	
	else:
		# Órbita parabólica (e == 1)
		ialpha = 0

	return ialpha

def criar_param(nome, t0, tstop, dt, dtout, dtdump, L1, L2, L3, L4, L5, L6,
	rmin, rmax, rmaxu, qmin, lclose, binary_outputfile, status_flag_for_open_statements):

	param = open(nome, 'w')
	param.write('{} {} {}\n'.format(t0, tstop, dt))
	param.write('{} {}\n'.format(dtout, dtdump))
	param.write('{} {} {} {} {} {}\n'.format(L1, L2, L3, L4, L5, L6))
	param.write('{} {} {} {} {}\n'.format(rmin, rmax, rmaxu, qmin, lclose))
	param.write(binary_outputfile + '\n')
	param.write(status_flag_for_open_statements)
	param.close()
	
def grade_Troianos(pos, apmax, apmin, apmed, incp, capomp, omegap, capmp):

	sigma_a = apmax - apmin
	delta_a = 5*sigma_a
	
	inc = incp
	capom = capomp
	omega = 60 + omegap
	capm = capmp

	# Redução de omega no intervalo 0 - 360 graus
	if omega > 360:
		nper = int(omega/360)
		omega = omega - nper*360

	# Passos no semieixo maior e na excentricidade da grade
	diva = 2*delta_a/99
	dive = 0.5/99
	
	# Semieixo maior inicial da grade
	a = apmed - delta_a
	n = 1

	data = open('Troianos_' + pos + '.txt', 'w') 
	for i in range(100):
	
		# Excentricidade inicial da grade
		e = 0
		for j in range(100):
		
			data.write("{:5d}  {:17.15f}  {:7.5f}  {:5.1f}  {:5.1f}  {:5.1f}  {:5.1f}\n".format(n, a, e, inc, capom, omega, capm))
			n = n + 1
			e = e + dive
		
		a = a + diva
 
	data.close()

def sistemas_reais(sistema):
	
	data = open(sistema + '.txt', 'r')
	lines = data.read().splitlines()
	data.close()

	# Criando arquivo param.in pro SyMBA e rmvs
	# -------------------------------------------------------------------------------------------------------------------------------

	# Periodo do planeta mais interno em dias
	P = float(lines[4].split()[3])**(3/2)*365.25

	# Radio da estrela em au
	R = float(lines[3].split()[1])*695508/149597870.700 

	t0 = "0.0d0"; tstop = 20*10**6*365.25; dt = "%.2f"%(P/20)
	dtout = "2.d4"; dtdump = "2.d4"
	L1 = "F"; L2 = "T"; L3 = "F"; L4 = "F"; L5 = "T"; L6 = "F"
	rmin = R; rmax = "10."; rmaxu = "-1."; qmin = R; lclose = "T"
	binary_outputfile = "bin.dat"
	status_flag_for_open_statements = "unknown"

	criar_param('param.in_SyMBA', t0, tstop, dt, dtout, dtdump, L1, L2, L3, L4, L5, L6,
	rmin, rmax, rmaxu, qmin, lclose, binary_outputfile, status_flag_for_open_statements)	
	
	if sistema == 'Kepler-9':
		tstop = '1.d7'
		dtout = '1.d3'
		dtdump = '1.d3'
	else:
		tstop = '1.d8'
		dtout = '1.d4'
		dtdump = '1.d4'
	
	criar_param('param.in_rmvs', t0, tstop, dt, dtout, dtdump, L1, L2, L3, L4, L5, L6,
	rmin, rmax, rmaxu, qmin, lclose, binary_outputfile, status_flag_for_open_statements)	
	
	# Criando arquivo pl.in pro SyMBA e rmvs
	# -------------------------------------------------------------------------------------------------------------------------------
	
	# Massa da estrela em unidades tal que G = 1 (uma massa solar equivale a 2.96*10**-4 unidades de massa)	
	gm = float(lines[3].split()[0])*2.959139768995959*10**-4

	pl_SyMBA = open('pl.in_SyMBA', 'w')
	pl_SyMBA.write(str(len(lines) - 3) + '\n')
	pl_SyMBA.write(str(gm) + '\n')
	pl_SyMBA.write('{} {} {}\n'.format(.0, .0, .0))
	pl_SyMBA.write('{} {} {}\n'.format(.0, .0, .0))
	
	pl_rmvs = open('pl.in_rmvs', 'w')
	pl_rmvs.write(str(len(lines) - 3) + '\n')
	pl_rmvs.write(str(gm) + '\n')
	pl_rmvs.write('{} {} {}\n'.format(.0, .0, .0))
	pl_rmvs.write('{} {} {}\n'.format(.0, .0, .0))

	# Salvando algumas CI dos planetas para usá-las na construção de grades de Troianos
	incp = []
	capomp = []
	omegap = []
	capmp = []

	for line in lines[4:]:
	
		text = line.split()
		# Massa do planeta em unidades de massa certa
		m = float(text[1])*1898.19*10**24*2.959139768995959*10**-4/(1988500*10**24)
		# Raio do planeta em au
		raio = float(text[2])*71492/149597870.700
		a = float(text[3])
		e = float(text[4])
		inc = float(text[5])*math.pi/180
		capom = float(text[6])*math.pi/180
		omega = float(text[7])*math.pi/180
		capm = float(text[8])*math.pi/180
		rhill = a*(1 - e)*pow(m/(3*gm), 1/3)	
		ialpha = detialpha(e)

		# orbel_el2xv. Essa função precisa que os ângulos estem em radianos
		x, y, z, vx, vy, vz = orbel_el2xv.orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm)

		pl_SyMBA.write('{:<30} {:<30} {}\n'.format(m, rhill, raio))
		pl_SyMBA.write('{:<30} {:<30} {}\n'.format(x, y, z))
		pl_SyMBA.write('{:<30} {:<30} {}\n'.format(vx, vy, vz))
		
		pl_rmvs.write('{:<30} {}\n'.format(m, raio))
		pl_rmvs.write('{:<30} {:<30} {}\n'.format(x, y, z))
		pl_rmvs.write('{:<30} {:<30} {}\n'.format(vx, vy, vz))
	
		incp.append(float(text[5]))
		capomp.append(float(text[6]))
		omegap.append(float(text[7]))
		capmp.append(float(text[8]))
	
	pl_SyMBA.close()
	pl_rmvs.close()

	# Criando arquivo start.in pro rmvs e SyMBA
	# -------------------------------------------------------------------------------------------------------------------------------

	start = open('start.in_rmvs', 'w')
	start.write('param.in\n')
	start.write('pl.in\n')
	start.write('tp.in\n')
	start.write('0.0')
	start.close()

	os.system("sed -e '3d' start.in_rmvs > start.in_SyMBA")

	# Executando a corrida do sistema 
	# -------------------------------------------------------------------------------------------------------------------------------
	
	os.chdir('/home/planeta9/Dropbox/Artigo_mestrado/Corridas/' + sistema)
	os.system('mv ../Scripts/*SyMBA .')
	os.system('mv start.in_SyMBA start.in')	
	os.system('mv param.in_SyMBA param.in')
	os.system('mv pl.in_SyMBA pl.in')
	os.system('/home/planeta9/Documents/swift/main/swift_symba7.x < start.in')	
	print('Corrida do Swift do sistema ' + sistema + ' finalizada')	
	os.system('cp ../Scripts/frequencia.dat .')
	os.system('/home/planeta9/Documents/swift/tools/follow_symba-all.x < frequencia.dat')
	
	# Plot 
	# -------------------------------------------------------------------------------------------------------------------------------
	
	data = open('follow.out', 'r')
	lines = data.read().splitlines()
	data.close()
	
	t = []
	a2 = []
	e2 = []
	inc2 = []
	a3 = []
	e3 = []
	inc3 = []
	a4 = []
	e4 = []
	inc4 = []	
	
	for line in lines:
	
		if len(line) < 30:
		
			t.append(float(line.split()[0]))
			
		elif line.split()[0] == "-2":
		
			a2.append(float(line.split()[1]))
			e2.append(float(line.split()[2]))
			inc2.append(float(line.split()[3]))
		
		elif line.split()[0] == "-3":
		
			a3.append(float(line.split()[1]))
			e3.append(float(line.split()[2]))
			inc3.append(float(line.split()[3]))
			
		else:
		
			a4.append(float(line.split()[1]))
			e4.append(float(line.split()[2]))
			inc4.append(float(line.split()[3]))

	print("Todos os planetas sobrevivem no tempo de integração?", len(t) == len(a2) == len(a3) == len(a4))
		
	# Convertimos dias a anos
	t = list(map(lambda i: i/365.25, t))
	
	fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True)
	fig.suptitle(sistema)
	ax1.plot(t, a2, marker = '', linestyle = '-', color = 'cornflowerblue', linewidth = 0.4)
	ax1.plot(t, a3, marker = '', linestyle = '-', color = 'darkred', linewidth = 0.4)
	ax1.plot(t, a4, marker = '', linestyle = '-', color = 'gray', linewidth = 0.4)
	ax1.set_ylabel('$a$ (au)')
	ax2.plot(t, e2, marker = '', linestyle = '-', color = 'cornflowerblue', linewidth = 0.4)
	ax2.plot(t, e3, marker = '', linestyle = '-', color = 'darkred', linewidth = 0.4)
	ax2.plot(t, e4, marker = '', linestyle = '-', color = 'gray', linewidth = 0.4)
	ax2.set_ylabel('$e$')
	ax3.plot(t, inc2, marker = '', linestyle = '-', color = 'cornflowerblue', linewidth = 0.4, alpha = 0.4)
	ax3.plot(t, inc3, marker = '', linestyle = '-', color = 'darkred', linewidth = 0.4, alpha = 0.4)
	ax3.plot(t, inc4, marker = '', linestyle = '-', color = 'gray', linewidth = 0.4, alpha = 0.4)
	ax3.set_ylabel('$i$ (deg)')
	ax3.set_xlabel('time (yr)')
	ax1.minorticks_on()
	ax2.minorticks_on()
	ax3.minorticks_on()
	ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')	
	ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
	ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
	ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')	
	ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
	ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
	plt.savefig(sistema + '.pdf')
	
	# Construindo uma grade de Troianos na órbita de cada planeta
	# -------------------------------------------------------------------------------------------------------------------------------

	pos = ['i', 'm', 'o']	
	amax = [max(a2), max(a3), max(a4)]
	amin = [min(a2), min(a3), min(a4)]
	amed = [np.mean(a2), np.mean(a3), np.mean(a4)]
		
	for i in range(3):
	
		grade_Troianos(pos[i], amax[i], amin[i], amed[i], incp[i], capomp[i], omegap[i], capmp[i])
	
	# Construindo diretórios de corridas com Troianos
	# -------------------------------------------------------------------------------------------------------------------------------

	os.system('mv /home/jessi/Dropbox/Artigo_mestrado/Scripts/*rmvs .')

	corridas = ['KeplerXi_i', 'KeplerXim_i', 'KeplerXio_i', 'KeplerXimo_i',
	'KeplerXm_m', 'KeplerXim_m', 'KeplerXmo_m', 'KeplerXimo_m',
	'KeplerXo_o', 'KeplerXmo_o', 'KeplerXio_o', 'KeplerXimo_o']

	[os.system('mkdir ' + corrida) for corrida in corridas]
	
	# Criando diferentes versões do arquivo original 'pl.in_rmvs' considerando todas as combinações possíveis entre planetas
	os.system("sed -e '8,$d' pl.in_rmvs > pl.in_i")
	os.system("sed -e '5,7d;11,$d' pl.in_rmvs > pl.in_m")
	os.system("sed -e '5,10d' pl.in_rmvs > pl.in_o")
	os.system("sed -e '11,$d' pl.in_rmvs > pl.in_im")
	os.system("sed -e '8,10d' pl.in_rmvs > pl.in_io")
	os.system("sed -e '5,7d' pl.in_rmvs > pl.in_mo")	
	
	for corrida in corridas:

		os.system('cp start.in_rmvs ' + corrida + '/start.in')
		os.system('cp param.in_rmvs ' + corrida + '/param.in')
		
		if corrida[7:10] == 'imo':
			os.system('cp pl.in_rmvs ' + corrida + '/pl.in')
		elif corrida[7:9] == 'im':
			os.system('cp pl.in_im ' + corrida + '/pl.in')
		elif corrida[7:9] == 'io':
			os.system('cp pl.in_io ' + corrida + '/pl.in')
		elif corrida[7:9] == 'mo':
			os.system('cp pl.in_mo ' + corrida + '/pl.in')
		elif corrida[7] == 'i':
			os.system('cp pl.in_i ' + corrida + '/pl.in')	
		elif corrida[7] == 'm':
			os.system('cp pl.in_m ' + corrida + '/pl.in')
		else:
			os.system('cp pl.in_o ' + corrida + '/pl.in')
					
		if corrida[-1] == 'i':
			os.system('cp Troianos_i.txt ' + corrida + '/Troianos.txt')
		elif corrida[-1] == 'm':
			os.system('cp Troianos_m.txt ' + corrida + '/Troianos.txt')
		else:
			os.system('cp Troianos_o.txt ' + corrida + '/Troianos.txt')
	
	os.system('rm *rmvs')
	
	# Modificando os nomes dos diretórios criados
	if sistema == 'Kepler-9':
		for corrida in corridas:
			corrida2 = corrida.replace('X', '9').replace('i', 'd').replace('m', 'b').replace('o', 'c')
			os.system('mv ' + corrida + ' ' + corrida2)	
	else:
		for corrida in corridas:
			corrida2 = corrida.replace('X', '56').replace('i', 'b').replace('m', 'c').replace('o', 'd')
			os.system('mv ' + corrida + ' ' + corrida2)	
	
# Script
# ---------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

	t0 = time.time()

	sistema = 'Kepler-9'
	sistemas_reais(sistema)

	t1 = time.time()
	seconds = t1 - t0
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print("h:min:s =", "%d:%02d:%02d" %(h, m, s))
