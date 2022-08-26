import math
import matplotlib.pyplot as plt
import numpy as np
from sympy import*
import pandas as pd
#from IPython.display import display

#Discretización  de la función Fuerza vs Tiempo
def Fuerza(tiempo):
    tuplas = [(0.0,0.0),(0.005,14),(0.010,29),(0.015,40),(0.020,52),(0.025,63),(0.03,73),(0.035,82),(0.04,89),(0.045,94),(0.05,99),(0.055,103),(0.06,106),(0.065,109),(0.07,110),(0.075,111),(0.08,111),(0.085,110),(0.09,109),(0.095,107),(0.100,105),(0.105,102),(0.110,100),(0.115,98),(0.120,96),(0.125,92),(0.130,89),(0.135,87),(0.140,85),(0.145,83),(0.150,82),(0.155,81),(0.160,80),(0.165,80),(0.170,80),(0.175,80),(0.180,81),(0.185,82),(0.190,83),(0.195,88),(0.200,91),(0.205,94),(0.210,96),(0.215,98),(0.220,100),(0.225,101),(0.230,103),(0.235,104),(0.240,105),(0.245,106),(0.250,107),(0.255,107),(0.260,107),(0.265,106),(0.270,105),(0.275,103),(0.280,100),(0.285,97),(0.290,92),(0.295,84),(0.300,76),(0.305,65),(0.310,56),(0.315,46),(0.320,35),(0.325,22),(0.330,0.0),(0.5,0.0)]
    tiempos = [t[0] for t in tuplas]
    fuerzas = [f[1] for f in tuplas]
    while(tiempo > 0.5): tiempo -= 0.5
    t0 = tiempo

    if t0 in tiempos:
        return fuerzas[tiempos.index(t0)]
    else:
        tup2 = [tup for tup in tuplas if (tup[0]>=t0)][0]
        tup1 = tuplas[tuplas.index(tup2)-1 if tuplas.index(tup2)!= 0 else 0]
        
        #Función de Interpolación
        # y = [(y2 - y1)/(x2 - x1)]*(x0 - x1) + y1
        t2,f2 = tup2
        t1,f1 = tup1

        return round((f2-f1)*(t0-t1)/(t2-t1)+f1,3)

def angle_time(angle):
    while(angle>360): angle -= 360
    return round(angle*(0.5/360.0),3)

'''def angle_rad(angle):
    return angle*(math.pi/180.0)'''


# Gráfica de la Fuerza vs Tiempo
init_angle = 0
finish_angle = 360
step_angle = 2
ángulos = np.arange(init_angle, finish_angle, step_angle)
tiempos = np.array(list(map(angle_time,ángulos)))
fuerzas = np.array(list(map(Fuerza,tiempos)))

fig1, ax1_1 = plt.subplots()
fig1_2, ax1_2 = plt.subplots()

ax1_1.plot(ángulos,fuerzas, color='#3F0B7D', linestyle='solid', label='Fuerza')
ax1_1.set_title("Fuerza humana al caminar")
ax1_1.set_ylabel("Fuerza [lb]")
ax1_1.set_xlabel("Tiempo [s]")
ax1_1.legend()
ax1_1.grid(True)


# Gráfica Presión-volumen del sistema
# Relación Presión/Volumen: P/V = 30
def Presión_volumen(volumen, relación=30, min_val=0, max_val=100):
    return min(max_val, max(min_val, relación*volumen))

volumenes = np.array(range(0,4,1))
presiones_vol = np.array(list(map(Presión_volumen,volumenes)))

ax1_2.plot(volumenes, presiones_vol,color='#367D0B', linestyle='dotted', label='Presión')
ax1_2.set_title("¨Presión-Volumen del Sistema")
ax1_2.set_ylabel("Presión de la bomba [psi]")
ax1_2.set_xlabel("Volumen desplazado [plug^3]")
ax1_2.legend()
ax1_2.grid(True)

# Se busca determinar el Desplazamiento
# Como volumen desplazado se considera un cilindro
# r: Radio del cilindro
# x: distancia de desplazamiento
# Volumen_des = Area_cir*desplazamiento
# Volumen_des = pi*r^2*x
# despejando el Desplazamiento
# x = V/pi*r^2

# Se conoce la relación Presión/Volumen=30
# V = P/30
# La presión es equivalente a un fuerza aplicada sobre un área.
# Presión = Fuerza/Área
# se conoce la Fuerza humana al caminar
# Área transversal del cilindro desplazado: A=pi*r^2

def Área(radio):
    return math.pi*radio**2

def Presión(fuerza, área):
    return fuerza/área

def Volumen_des(presión, relación=30, min_val=0, max_val=100):
    return min(max_val, max(min_val, presión/relación))

def Desplazamiento(volumen_des,área):
    return volumen_des/área

radio = 1.181
#Radio del ayudante: 1 in ; 2.54 cm

áreas = np.repeat(Área(radio),len(fuerzas))
#Se asume el talón promedio en un radio de 3 cm = 1.18 in
presiones = np.array(list(map(Presión, fuerzas, áreas)))
volumenes_des = np.array(list(map(Volumen_des, presiones)))
desplazamientos = np.array(list(map(Desplazamiento, volumenes_des, áreas)))
#desplazamientos = np.array(list(map(desplazamiento, fuerzas)))

fig2, (ax2_1, ax2_2, ax2_3) = plt.subplots(nrows=3, ncols=1, sharex=True)

fig3, ax3 = plt.subplots()
ax3.plot(ángulos, desplazamientos, color='#901B12')
ax3.set_title("Desplazamiento vs Ángulo")
ax3.set_xlabel("Ángulo [°]")
ax3.set_ylabel("Desplazamiento [in]")
ax3.grid(True)

ax2_1.plot(ángulos, presiones)
ax2_1.set_title("Presión vs Ángulo")
#ax2_1.set_xlabel("Ángulo [°]")
ax2_1.set_ylabel("Presión [psi]")
ax2_1.grid(True)

ax2_2.plot(ángulos, volumenes_des)
ax2_2.set_title("Volumen vs Ángulo")
ax2_2.set_xlabel("Ángulo [°]")
ax2_2.set_ylabel("Volumen [pulg^3]")
ax2_2.grid(True)

ax2_3.plot(ángulos, desplazamientos, color='#921B18')
ax2_3.set_title("Desplazamiento vs Ángulo")
ax2_3.set_xlabel("Ángulo [°]")
ax2_3.set_ylabel("Desplazamiento [in]")
ax2_3.grid(True)

fig1.savefig("Figura1")
fig2.savefig("Figura2")

datos = {'col1': ángulos, 'col2': desplazamientos}
dtf = pd.DataFrame(data = datos)
dtf.columns = ['Ángulos', 'Desplazamientos']

for i in range(len(fuerzas)):
    print("{}  - {}  -  {}".format(ángulos[i], desplazamientos[i], fuerzas[i]))

plt.tight_layout()
plt.show()