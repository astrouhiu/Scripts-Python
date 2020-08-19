'''
- Salvamos dados de TNOs com a_inf au <= a <= a_sup au, q >= q_scat au e bem caracterizados (mais do que 1 opps.)
- Salvamos dados de TNOs com a_inf au <= a <= a_sup au, q >= q_det au e bem caracterizados (mais do que 1 opps.)
- Plotamos elementos angulares vs. a dos TNOs e Centauros bem caracterizados (mais do que 1 opps.), destacando os TNOs distantes
- Determinamos a magnitude visual do TNO real mais brilhante na região a_inf au <= a <= a_sup au e q >= q_scat
'''

from matplotlib.ticker import AutoMinorLocator
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex = True)
rc('legend', fontsize = 12)
rc('font', **{'weight': 'light','serif':['Times']}, size = 14)

# Dados de Entrada:
# =======================================================================================================================================
#
# file_TNOs:		rota a arquivo TNOs.txt
# file_Centaurs:	rota a arquivo Centaurs.txt
# a_inf:		limite inferior em 'a' (au) para aplicar OBIP
# a_sup:		limite superior em 'a' (au) para aplicar OBIP
# q_scat:		'q' (au) dos objetos espalhados
# q_det:		'q' (au) dos objetos destacados
# amin:		limite inferior em 'a' (au) para os TNOs distantes
# qmin:		limite inferior em 'q' (au) para os TNOs distantes
# aproximacao:		método usado para determinar a distância ao Sol dado como string, pode ser: 'iteracao' ou 'serie'

file_TNOs = '/home/jessi/Dropbox/MPC/Distantes/TNOs.txt'
file_Centaurs = '/home/jessi/Dropbox/MPC/Distantes/Centaurs.txt'
a_inf = 100
a_sup = 200
q_scat = 30
q_det = 40
amin = 250
qmin = 40
aproximacao = 'serie'
# =======================================================================================================================================

def zero360(angle):

	# angle é um elemento e está em graus	
	
	nper = int(angle/360)	
	angle = angle - nper*360

	if (angle < 0):

		angle = angle + 360

	return angle

def visual_mag(a, e, M, H, aproximacao):

	# a: semieixo maior em au
	# e: excentricidade
	# M: anomalia média em radianos
	# aproximacao: método usado para determinar a distância ao Sol
	# phi: anomalia verdadeira, em radianos
	# dS: distancia ao Sol em au
	# dE: distancia à Terra em au 
	# H: mag. absoluta
	# V: mag. visual

	if aproximacao == 'iteracao':

		epsilon = math.pow(10,-10)

		x_ant = (e*math.sin(M))/(1 - e*math.cos(M))
		x = e*math.sin(M + x_ant)	

		while sqrt((x - x_ant)**2) > epsilon:
	
			x_ant = x
			x = e*math.sin(M + x)

		E = M + x

		y = sqrt(1 - e**2)*math.sin(E)/(1 - e*math.cos(E))
		x = (math.cos(E) - e)/(1 - e*math.cos(E))
		phi = math.atan2(y, x)
		
		dS = a*(1 - e**2)/(1 + e*math.cos(phi))

	if aproximacao == 'serie':

		dS = a*( 1 + math.pow(e,2)/2 + (-e + 3*math.pow(e,3)/8)*math.cos(M) - math.pow(e,2)/2*math.cos(2*M) - 3*math.pow(e,3)/8*math.cos(3*M) )

	# A seguinte relação corresponde a uma observação em oposição, uma boa aproximação para a observação de objetos distantes

	dE = dS - 1
	
	V = H + 5*math.log10(dS*dE)

	return V

def dados(arquivo):

	# Extraimos dados de objetos bem caracterizados

	data = open(arquivo, 'r')
	lines = data.read().splitlines()
	data.close()

	name = []
	H = [] 
	a = []
	e = []
	q = []
	inc = []
	omega = []
	capom = []
	obar = []
	capm = []
	opps = []
	
	for line in lines[2:]:

		if (line[113:114] != "("):

			name.append(line[27:37])
			H.append(float(line[55:60]))
			a.append(float(line[103:111]))
			e.append(float(line[98:103]))
			q.append(float(line[38:45]))
			inc.append(float(line[92:97]))
			omega.append(float(line[80:85]))
			capom.append(float(line[86:91]))
			obar.append(omega[-1] + capom[-1])
			capm.append(float(line[73:78]))
			opps.append(float(line[114:118]))

	obar = list(map(lambda x: zero360(x), obar))
	
	return name, H, a, e, q, inc, omega, capom, obar, capm, opps

# ---------------------------------------------------------------------------------------------------------------------------------------
# Criamos arquivos espalhados-destacados.txt e destacados.txt

name_T, H_T, a_T, e_T, q_T, inc_T, omega_T, capom_T, obar_T, capm_T, opps_T = dados(file_TNOs)
name_C, H_C, a_C, e_C, q_C, inc_C, omega_C, capom_C, obar_C, capm_C, opps_C = dados(file_Centaurs)

# Fusionamos ambas as listas

name_all = name_C + name_T
H_all = H_C + H_T
a_all = a_C + a_T
e_all = e_C + e_T
q_all = q_C + q_T
inc_all = inc_C + inc_T
omega_all = omega_C + omega_T
capom_all = capom_C + capom_T
obar_all = obar_C + obar_T
capm_all = capm_C + capm_T
opps_all = opps_C + opps_T

# Salvamos dados de TNOs com a_inf au <= a <= a_sup au, q >= q_scat au e bem caracterizados

name_s = [x for x,y,z in zip(name_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
a_s = [x for x,y,z in zip(a_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
e_s = [x for x,y,z in zip(e_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
q_s = [x for x,y,z in zip(q_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
inc_s = [x for x,y,z in zip(inc_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
capom_s = [x for x,y,z in zip(capom_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
omega_s = [x for x,y,z in zip(omega_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
capm_s = [x for x,y,z in zip(capm_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]
H_s = [x for x,y,z in zip(H_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_scat]

with open("espalhados-destacados.txt", "w") as f:
	f.writelines(map("{:10s}  {:5.2f}  {:5.3f}  {:6.3f}  {:5.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:3.1f}\n".format, name_s, a_s, e_s, q_s, inc_s, capom_s, omega_s, capm_s, H_s)) 

# Salvamos dados de TNOs com a_inf au <= a <= a_sup au, q >= q_det au e bem caracterizados

name_d = [x for x,y,z in zip(name_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
a_d = [x for x,y,z in zip(a_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
e_d = [x for x,y,z in zip(e_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
q_d = [x for x,y,z in zip(q_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
inc_d = [x for x,y,z in zip(inc_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
capom_d = [x for x,y,z in zip(capom_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
omega_d = [x for x,y,z in zip(omega_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
capm_d = [x for x,y,z in zip(capm_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]
H_d = [x for x,y,z in zip(H_all, a_all, q_all) if y >= a_inf and y <= a_sup and z >= q_det]

with open("destacados.txt", "w") as f:
	f.writelines(map("{:10s}  {:5.2f}  {:5.3f}  {:6.3f}  {:5.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:3.1f}\n".format, name_d, a_d, e_d, q_d, inc_d, capom_d, omega_d, capm_d, H_d))

# ---------------------------------------------------------------------------------------------------------------------------------------
# Salvamos dados de TNOs distantes com a >= amin, q >= qmin e bem caracterizados 

name_dist = [x for x,y,z in zip(name_all, a_all, q_all) if y >= amin and z >= qmin]
a_dist = [x for x,y,z in zip(a_all, a_all, q_all) if y >= amin and z >= qmin]
e_dist = [x for x,y,z in zip(e_all, a_all, q_all) if y >= amin and z >= qmin]
q_dist = [x for x,y,z in zip(q_all, a_all, q_all) if y >= amin and z >= qmin]
inc_dist = [x for x,y,z in zip(inc_all, a_all, q_all) if y >= amin and z >= qmin]
capom_dist = [x for x,y,z in zip(capom_all, a_all, q_all) if y >= amin and z >= qmin]
omega_dist = [x for x,y,z in zip(omega_all, a_all, q_all) if y >= amin and z >= qmin]
obar_dist = [x for x,y,z in zip(obar_all, a_all, q_all) if y >= amin and z >= qmin]
opps_dist = [x for x,y,z in zip(opps_all, a_all, q_all) if y >= amin and z >= qmin]

with open("TNOs-distantes-helio-eclip.txt", "w") as f:
	f.writelines(map("{:10s}  {:6.1f}  {:5.3f}  {:6.3f}  {:5.1f}  {:5.1f}  {:5.1f}  {:5.1f}  {:5.0f}\n".format, name_dist, a_dist, e_dist, q_dist, inc_dist, capom_dist, omega_dist, obar_dist, opps_dist)) 

# ---------------------------------------------------------------------------------------------------------------------------------------
# Plot de todos os objetos com q >= qmin e bem caracterizados. Destacamos em vermelho os TNOs distantes (aqueles com a >= amin)

a_qgtqmin = [x for x,y in zip(a_all, q_all) if y >= qmin]
capom_qgtqmin = [x for x,y in zip(capom_all, q_all) if y >= qmin]
omega_qgtqmin = [x for x,y in zip(omega_all, q_all) if y >= qmin]
obar_qgtqmin = [x for x,y in zip(obar_all, q_all) if y >= qmin]

f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
f.subplots_adjust(hspace = 0.15)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible = False)
ax1.axvline(x = 250, color = 'gray', linestyle = '--', linewidth = 1)
ax2.axvline(x = 250, color = 'gray', linestyle = '--', linewidth = 1) 
ax3.axvline(x = 250, color = 'gray', linestyle = '--', linewidth = 1)
ax1.plot(a_qgtqmin, capom_qgtqmin, color = 'darkgrey', marker = 'o', markeredgewidth = 0.0, alpha = 0.5, linestyle = 'None',  markersize = 4)
ax2.plot(a_qgtqmin, omega_qgtqmin, color = 'darkgrey', marker = 'o', markeredgewidth = 0.0, alpha = 0.5, linestyle = 'None',  markersize = 4)
ax3.plot(a_qgtqmin, obar_qgtqmin, color = 'darkgrey', marker = 'o', markeredgewidth = 0.0, alpha = 0.5, linestyle = 'None',  markersize = 4)
ax1.plot(a_dist, capom_dist, color = 'lightcoral', marker = 'o', markeredgewidth = 0.0, linestyle = 'None', markersize = 5)
ax2.plot(a_dist, omega_dist, color = 'lightcoral', marker = 'o', markeredgewidth = 0.0, linestyle = 'None', markersize = 5)
ax3.plot(a_dist, obar_dist, color = 'lightcoral', marker = 'o', markeredgewidth = 0.0, linestyle = 'None', markersize = 5)
for i in range(len(a_dist)):
	ax1.text(a_dist[i] + 14, capom_dist[i], str("%.0f"%q_dist[i]), fontsize = 5, va = 'center')
	ax2.text(a_dist[i] + 14, omega_dist[i], str("%.2f"%e_dist[i]), fontsize = 5, va = 'center')
	ax3.text(a_dist[i] + 14, obar_dist[i], str("%.0f"%inc_dist[i]), fontsize = 5, va = 'center')
ax1.text((1300 - 30)/2, 290, r"$q\geq40$ au", ha = "center")
ax1.set_ylabel('$\Omega$ (deg)')
ax2.set_ylabel('$\omega$ (deg)')
ax3.set_ylabel('$\\varpi$ (deg)')
ax3.set_xlabel('$a$ (au)')
ax3.set_ylim(0, 360)
ax3.set_xlim(30, 1300)
ax1.minorticks_on()
ax2.minorticks_on()
ax3.minorticks_on()
ax1.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
ax1.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
ax2.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
ax2.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
ax3.tick_params(direction = 'in', bottom = True, top = True, left = True, right = True, which = 'major')
ax3.tick_params(direction = 'in', bottom = False, top = False, left = True, right = True, which = 'minor')
ax1.set_xticks(np.linspace(200, 1200, 6))
ax1.set_yticks(np.linspace(0, 300, 4))
plt.savefig('TNOs-distantes-MPC-helio-eclip.pdf', format = 'pdf', bbox_inches = 'tight')

# ---------------------------------------------------------------------------------------------------------------------------------------
# Determinamos a mag. visual do TNO mais brilhante (min(V)) entre aqueles com a_inf au <= a <= a_sup au e q >= q_scat au, usando dados de
# H (ao invés de raio e albedo) do MPC

# Ângulos devem estar em radianos

capm_s2 = [x*math.pi/180 for x in capm_s]
aproximacao = [aproximacao]*len(capm_s2)
V = list(map(visual_mag, a_s, e_s, capm_s2, H_s, aproximacao))

print('******* Amostra Observacional *******', 'Dados segundo o Minor Planet Center', sep = "\n")
print('- Núm. de objetos com', a_inf, 'au <= a <=', a_sup, 'au e q >=', q_scat, 'au =', len(a_s))
print('- Núm. de objetos com', a_inf, 'au <= a <=', a_sup, 'au e q >=', q_det, 'au =', len(a_d))
print('- A magnitude visual do objeto mais brilhante dentre', a_inf, 'au <= a <=', a_sup, 'au e q >=', q_scat, 'au é', min(V), 'correspondente ao objeto', name_s[V.index(min(V))])
print('- A magnitude visual do objeto mais fraco dentre', a_inf, 'au <= a <=', a_sup, 'au e q >=', q_scat, 'au é', max(V), 'correspondente ao objeto', name_s[V.index(max(V))])
print('- Núm. de objetos com', a_inf, 'au <= a <=', a_sup, 'au e q >=', q_scat, 'au e V < 21 =', len([x for x in V if x < 21]))
