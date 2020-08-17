'''
Esse script considera dados dos planetas e TNOs observados para uma mesma data e cria os seguintes arquivos

- Na rota path: names.sh
		arquivos dos planetas e TNOs com extensao .txt
		J2000.txt, compilando elementos orbitais dos planetas e TNOs no formato: a e inc g n M (a e inc omega capom capm)
		big.in, big.ina700, big.ina1500 e 'n' arquivos small.in no iv4. Todos no formato asteroidal: a e inc g n M (a e inc omega
		capom capm). Arquivos big.ina700 e big.ina1500 incluim dados do planeta 9
'''

import numpy as np
import math
import os
# Librerias criadas com codigos Fortran do Swift
import xv2el
import el2xv

# Dados de Entrada
# =======================================================================================================================================
#
# x:		So objetos observados em mais do que uma opossicao? (True ou False)
# path:	Rota para os arquivos de entrada: TNOs.txt, osc_tbl (como executable) e osc_tbl.inp; e contera os arquivos de saida
# gm:		G vezes a massa central. Aqui G = 1 e 1 massa solar é aprox. 2.96e-4 unidades de massa.
# mpl:		Lista que contem a massa dos 4 planetas gigantes em solar masses. 

x = True
path = '/home/jessi/Dropbox/Codigos_Projeto/CI-TNOs/'
gm = 0.01720209895**2
mpl = [0.000954786104, 0.000285836787, 4.37273165e-05, 5.17759138e-05]
mpl = np.float32(np.multiply(gm, mpl))
# =======================================================================================================================================

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

# ---------------------------------------------------------------------------------------------------------------------------------------
# Criamos um arquivo em Bash em um formato específico para executar o script do JPL para os planetas e TNOs, e o executamos. Os nomes dos
# TNOs são obtidos do arquivo TNOs.txt do MPC

inputfile = open(path + 'TNOs.txt', 'r')
lines = inputfile.read().splitlines()
inputfile.close()

outputfile = open(path + 'names.sh', 'w')
outputfile.write('#!/bin/bash' + '\n')
outputfile.write('\n')
outputfile.write("./osc_tbl " + '"[599]" ' + 'Jup.txt' + '\n')
outputfile.write("./osc_tbl " + '"[699]" ' + 'Sat.txt' + '\n')
outputfile.write("./osc_tbl " + '"[799]" ' + 'Ura.txt' + '\n')
outputfile.write("./osc_tbl " + '"[899]" ' + 'Net.txt' + '\n')

name2 = ['Jup', 'Sat', 'Ura', 'Net']

if x == True:

	for line in lines[2:]:

		if (line[112] != "("):
	
			if line[141:146] == 'Pluto':

				name = '[999]'
				name2.append('Pluto')

			else:

				name = line[27:37]
				name2.append(name.replace(" ", ""))

			outputfile.write("./osc_tbl " + '"' + name + '" ' + name2[-1] + ".txt" + '\n')

else:

	for line in lines[2:]:

		if line[141:146] == 'Pluto':

			name = '[999]'
			name2.append('Pluto')

		else:

			name = line[27:37]
			name2.append(name.replace(" ", ""))

		outputfile.write("./osc_tbl " + '"' + name + '" ' + name2[-1] + ".txt" + '\n')

outputfile.close()

os.chdir(path)
os.system('chmod +x names.sh')
os.system('./names.sh')

# ---------------------------------------------------------------------------------------------------------------------------------------
# Criamos o arquivo 'J2000.txt' contendo os elementos orbitais para planetas e TNOs para uma mesma data, considerando o plano J2000.     
# Corpo central é o Sol e o plano é a ecliptica. Formato: name a e inc g n M (name a e inc omega capom capm)

outputfile = open(path + 'J2000.txt', 'w')

# Indice para Plutão

p = name2.index('Pluto')

for k in range(len(name2)):

	print(k)

	fp = open(path + name2[k] + '.txt')

	if k < 4 or k == p:

		for i, line in enumerate(fp):

			if i == 29:

				e = line.split()[1]
				inc = line.split()[5]
		
			if i == 30:

				capom = line.split()[1]
				omega = line.split()[4]

			if i == 31:

				capm = line.split()[4]

			if i == 32:

				a = line.split()[2]
				outputfile.write(k + "  " + a + "  " + e + "  " + inc + "  " + omega + "  " + capom + "  " + capm + '\n')
				break

	else:

		for i, line in enumerate(fp):

			if i == 42:

				e = line.split()[1]
				inc = line.split()[5]
		
			if i == 43:

				capom = line.split()[1]
				omega = line.split()[4]

			if i == 44:

				capm = line.split()[4]

			if i == 45:

				a = line.split()[2]
				outputfile.write(k + "  " + a + "  " + e + "  " + inc + "  " + omega + "  " + capom + "  " + capm + '\n')
				break

	fp.close()

outputfile.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
# Criamos n-arquivos small.in e 3 arquivos big.in, a partir do arquivo J2000.txt gerado anteriormente, em realacao ao iv4 heliocêntrico

inputfile = open(path + 'J2000.txt', 'r')
lines = inputfile.read().splitlines()
inputfile.close()

outputbig = open(path + 'big.in', 'w')
outputbig.write(')O+_06 Big-body initial data  (WARNING: Do not delete this line!!)' + '\n')
outputbig.write(") Lines beginning with `)' are ignored." + '\n')
outputbig.write(')---------------------------------------------------------------------' + '\n')
outputbig.write(' style (Cartesian, Asteroidal, Cometary) = Asteroidal' + '\n')
outputbig.write(' epoch (in days) = 0' + '\n')
outputbig.write(')---------------------------------------------------------------------' + '\n')

outputsmall = open(path + 'small.in1', 'w')
outputsmall.write(')O+_06 Small-body initial data  (WARNING: Do not delete this line!!)' + '\n')
outputsmall.write(") Lines beginning with `)' are ignored." + '\n')
outputsmall.write(')---------------------------------------------------------------------' + '\n')
outputsmall.write(' style (Cartesian, Asteroidal, Cometary) = Ast' + '\n')
outputsmall.write(')---------------------------------------------------------------------' + '\n')

pl = ['     Jup  m=0.000954786104  d=1.3', '     Sat  m=0.000285836787  d=0.7', '     Ura  m=4.37273165e-05  d=1.2', '     Net  m=5.17759138e-05  d=1.7']

k = 1
n = 0

hx = 0
hy = 0
hz = 0

for line in lines[0:4]:

	test = line.split()
	a = float(test[1])
	e = float(test[2])
	inc = float(test[3])*math.pi/180
	omega = float(test[4])*math.pi/180
	capom = float(test[5])*math.pi/180
	capm = float(test[6])*math.pi/180
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

h2 = np.float32(hx*hx + hy*hy + hz*hz)
h = np.float32(math.sqrt(h2))

if (hz > h):
		
	hz = h	
	hx = 0
	hy = 0

# math.acos(x): Return the arc cosine of x, in radians

incplano = np.float32(math.acos(hz/h))
fac = np.float32(math.sqrt(hx*hx + hy*hy)/h)
	
if (fac < 2**-127 or incplano == 0):

	alpha = 0

else:

	# math.atan2(y, x): Return atan(y/x), in radians. The result is between -pi and pi. 
	# The vector in the plane from the origin to point (x, y) makes this angle with the 
	# positive X axis. The point of atan2() is that the signs of both inputs are known 
	# to it, so it can compute the correct quadrant for the angle. For example, atan(1) 
	# and atan2(1, 1) are both pi/4, but atan2(-1, -1) is -3*pi/4. 
	
	alpha = np.float32(math.atan2(hx, -hy))

if (alpha < 0):

	alpha = np.float32(alpha + 2*math.pi)

# Aplicamos a transformacao ao plano iv4 a todos os corpos

n = 0

for line in lines:

	test = line.split()
	a = float(test[1])
	e = float(test[2])
	inc = float(test[3])*math.pi/180
	omega = float(test[4])*math.pi/180
	capom = float(test[5])*math.pi/180
	capm = float(test[6])*math.pi/180
	ialpha = detialpha(e)
		
	# el2xv. Nesses objetos não temos órbitas hiperbólicas
	
	x, y, z, vx, vy, vz = el2xv.orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm)
					
	# xv no iv4
	
	xi, yi, zi, vxi, vyi, vzi = rotation(np.float32(incplano), np.float32(alpha), np.float32(x), np.float32(y), np.float32(z), np.float32(vx), np.float32(vy), np.float32(vz))	

	if n < 4:

		gmsum = mpl[n] + gm

		# xv2el no iv4. Essa função retorna os ângulos em radianos
		
		_, ai, ei, inci, capomi, omegai, capmi = xv2el.orbel_xv2el(xi, yi, zi, vxi, vyi, vzi, gmsum)

		outputbig.write(pl[n] + '\n')
		outputbig.write(str("%6E" %ai) + "  " + str("%6E" %ei) + "  " + str("%6E" %(inci*180/math.pi)) + "  " + str("%6E" %(omegai*180/math.pi)) + "  " + str("%6E" %(capomi*180/math.pi)) + "  " + str("%6E" %(capmi*180/math.pi)) + '  0  0  0' + "\n")	

		n = n + 1

	else:

		gmsum = gm
	
		# xv2el no iv4. Essa função retorna os ângulos em radianos
		
		_, ai, ei, inci, capomi, omegai, capmi = xv2el.orbel_xv2el(xi, yi, zi, vxi, vyi, vzi, gmsum)

		outputsmall.write('p' + str(k) + '\n')
		outputsmall.write(str("%6E" %ai) + "  " + str("%6E" %ei) + "  " + str("%6E" %(inci*180/math.pi)) + "  " + str("%6E" %(omegai*180/math.pi)) + "  " + str("%6E" %(capomi*180/math.pi)) + "  " + str("%6E" %(capmi*180/math.pi)) + '  0  0  0' + "\n")	

		k = k + 1
	
		if k == 101:

			outputsmall.close()
			outputsmall = open(path + 'small.in' + str(j), 'w')
			outputsmall.write(')O+_06 Small-body initial data  (WARNING: Do not delete this line!!)' + '\n')
			outputsmall.write(") Lines beginning with `)' are ignored." + '\n')
			outputsmall.write(')---------------------------------------------------------------------' + '\n')
			outputsmall.write(' style (Cartesian, Asteroidal, Cometary) = Ast' + '\n')
			outputsmall.write(')---------------------------------------------------------------------' + '\n')
			
			k = 1

outputbig.close()
outputsmall.close()

# Para o pl9 consideramos os dados já no iv4

os.chdir(path)

os.system('cp big.in big.ina700')
os.system("echo '     pl9  r=1.00000E+00  d=1.00000E+00  m= 3e-5' >> big.ina700")
os.system("echo '700  0.914286  30  0  0  3  0  0  0' >> big.ina700")

os.system('cp big.in big.ina1500')
os.system("echo '     pl9  r=1.00000E+00  d=1.00000E+00  m= 5e-5' >> big.ina1500")
os.system("echo '1500  0.96  30  0  0  149  0  0  0' >> big.ina1500")
