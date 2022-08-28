import numpy as np
import cmath
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

#Tabla de Magnitudes y Ángulos

#Elemento -- |r_i| -- θ_i
# 1          9.625    317°
# 2          2        θ_2
# 3          8.375    θ_3
# 4          4.187    θ_4

r1 = 9.625 #(pulg)
r2 = 2 #(pulg)
r3 = 8.375 #(pulg)
r4 = 4.187 #(pulg)

θ1 = 0 #(rad)
θ2 = np.array([(math.pi/180)*x for x in range(0, 361, 2)]) #(rad)
ω2 = 52.360 #(rad/s)
α2 = 0 #(rad/s^2)
tiempos = θ2/ω2
length = len(θ2)

#r1_x = r1*math.cos(θ1)
#r1_y = r1*math.sin(θ1)
r1_x = r1

r2_x = np.array([r2*math.cos(i2) for i2 in θ2]) 
r2_y = np.array([r2*math.sin(i2) for i2 in θ2]) 

"""x = np.array([r3, r4])
y = np.array([r3, r4])

p_coe = [x, y]
p_solu = [r1_x-r2_x, r1_y-r2_y]
"""
#POSICIÓN

#(x): r3 * Cosθ3 - r4 * Cosθ4 = r1 * Cosθ1 - r2 * Cosθ2
#(y): r3 * Senθ3 - r4 * Senθ4 = r1 * Senθ1 - r2 * Senθ2

#(x): r3 * Cosθ3 = r4 * Cosθ4 - r2 * Cosθ2 + r1
#(y): r3 * Senθ3 = r4 * Senθ4 - r2 * Senθ2

#Se eleva al cuadrado ambas ecuaciones y se suman
# r3^2 = (-r2*Senθ2 + r4*Senθ4)^2 + (-r2*Cosθ2 + r4*Cosθ4 + r1)^2
# r3^2 = r2^2 + r4^2 + r1^2 -2*r2*r1*Cosθ2 + 2*r4*r1*Cosθ4 - 2*r2*r4*(Senθ2*Senθ4 + Cosθ2*Cosθ4)

#Dividiendo para 2*r2*r4
# (r1/r2)*Cosθ4 - (r1/r4)*Cosθ2 + (r2^2 - r3^2 + r4^2 + r1^2)/(2*r2*r4) = Senθ2*Senθ4 + Cosθ2*Cosθ4

#Se definen las constantes:
k1 = r1/r2
k2 = r1/r4
k3 = (r2**2 -r3**2 + r4**2 + r1**2)/(2*r2*r4)

#Se aplica Identidades Semiangulares para reemplazar Senθ4 y Cosθ4 en terminos de Tanθ4
# Senθ4 = (2*Tan(θ4/2))/(1 + Tan^2(θ4/2))
# Cosθ4 = (1 - Tan^2(θ4/2))/(1 + Tan^2(θ4/2))

#Las longitudes de los eslabones y el valor de entrada θ2 se reunen como constantes A, B y C.
# A*Tan^2(θ4/2) + B*Tan(θ4/2) + C = 0

def θ4(β2):
    A = math.cos(β2) - k1 -k2*math.cos(β2) + k3
    B = -2*math.sin(β2)
    C = k1 - (k2 + 1)*math.cos(β2) + k3

    # calculating  the discriminant
    dis = (B**2) - (4 * A*C)
    
    if dis >= 0:
        # find two results
        ans1 = (-B - math.sqrt(abs(dis)))/(2*A)
        ans2 = (-B + math.sqrt(abs(dis)))/(2*A)

        #Ángulo θ4 - Configuración Abierta:
        θ4_a = 2*math.atan(ans1)
        
        #Ángulo θ4 - Configuración Cerrada:
        θ4_c = 2*math.atan(ans2)

        #return(math.pi - θ4_a, θ4_c + math.pi)
        return(math.pi - θ4_a)
    else:
        print("(raices conjugadas complejas) -> La Cadena Cinemática está Abierta")
        return 0

#Determinando θ4:
#(x): r4 * Cosθ4 = r2 * Cosθ2  + r3 * Cosθ3 - r1
#(y): r4 * Senθ4 = r2 * Senθ2  + r3 * Senθ3

#Siguiendo el anterior proceso.
# k1*Cosθ3 + k4*Cosθ2 + k5 = Cosθ2*Cosθ3 + Senθ2*Senθ3
k4 = r1/r3
k5 = (r4**2 - r1**2 - r2**2 - r3**2)/(2*r2*r3)

#D*Tan^2(θ3/2)  E*Tan(θ3/2) + F = 0

def θ3(β2):
    D = math.cos(β2) - k1 -k4*math.cos(β2) + k5
    E = -2*math.sin(β2)
    F = k1 + (k4 - 1)*math.cos(β2) + k5

    # calculating  the discriminant
    dis = (E**2) - (4 * D*F)
    
    if dis >= 0:
        # find two results
        ans1 = (-E + math.sqrt(abs(dis)))/(2*D)
        ans2 = (-E - math.sqrt(abs(dis)))/(2*D)

        #Ángulo θ3 - Configuración Abierta:
        θ3_a = 2*math.atan(ans1)
        
        #Ángulo θ3 - Configuración Cerrada:
        θ3_c = 2*math.atan(ans2)
        #return(θ3_a, θ3_c)
        return(θ3_a)
    else:
        print("(raices conjugadas complejas) -> La Cadena Cinemática está Abierta")
        return 0

θ4 = np.array(list(map(θ4, θ2)))
θ3 = np.array(list(map(θ3, θ2)))

#Extremos
θ3_min = min(θ3)
θ3_max = max(θ3)
θ4_min = min(θ4)
θ4_max = max(θ4)

#Gráfica de Posiciones
fig1, ax1 = plt.subplots(figsize=(10, 8))
# Set axis ranges; by default this will put major ticks every 25.
ax1.set_xlim(0, 0.12)
ax1.set_ylim(-1, 6.5)

# Change major ticks to show every 20.
ax1.xaxis.set_major_locator(MultipleLocator(0.012))
ax1.yaxis.set_major_locator(MultipleLocator(0.5))

# Change minor ticks to show every 5. (20/4 = 5)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))

# Turn grid on for both major and minor ticks and style minor slightly
# differently.
ax1.grid(which='major', color='#CCCCCC', linestyle='--')
ax1.grid(which='minor', color='#CCCCCC', linestyle=':')

ax1.plot(tiempos, θ2, color='#C9002D', linestyle='solid', label='θ2')
ax1.plot(tiempos, θ3, color='#e22e13', linestyle='dotted', label='θ3', lw=2)
#ax1.axhline(θ3_min, color='#C0C0C0', lw=1)
#ax1.axhline(θ3_max, color='#C0C0C0', lw=1)
ax1.plot(tiempos, θ4, color='#ffb52b', linestyle='dotted', label='θ4', lw=2)
#ax1.axhline(θ4_min, color='#C0C0C0', lw=1)
#ax1.axhline(θ4_max, color='#C0C0C0', lw=1)
ax1.axhline(0, color='#303030', lw=1)
ax1.set_title("Posiciones")
ax1.set_ylabel("Ángulo [rad]")
ax1.set_xlabel("Tiempo [s]")
ax1.legend()
ax1.grid(True)

# Funciones de Velocidades Angulares
def ω3(β2, β3, β4):
    return (r2*ω2/r3)*(math.sin(β4-β2)/math.sin(β3-β4))

def ω4(β2, β3, β4):
    return (r2*ω2/r4)*(math.sin(β2-β3)/math.sin(β4-β3))

# Funciones de Velocidades Lineales
def V(r, θ, ω):
    return (r*ω)*np.array([-math.sin(θ), math.cos(θ)])

#Evaluación
ω3 = np.array(list(map(ω3, θ2, θ3, θ4)))
ω4 = np.array(list(map(ω4, θ2, θ3, θ4)))

#Gráfica de Velocidades
fig2, ax2 = plt.subplots(figsize=(10,8))

# Set axis ranges; by default this will put major ticks every 25.
ax2.set_xlim(0, 0.12)
ax2.set_ylim(-40, 60)

# Change major ticks to show every 20.
ax2.xaxis.set_major_locator(MultipleLocator(0.01))
ax2.yaxis.set_major_locator(MultipleLocator(10))

# Change minor ticks to show every 5. (20/4 = 5)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))

# Turn grid on for both major and minor ticks and style minor slightly
# differently.
ax2.grid(which='major', color='#CCCCCC', linestyle='--')
ax2.grid(which='minor', color='#CCCCCC', linestyle=':')

ax2.plot(tiempos, np.repeat(ω2,length), color='#0b1180', linestyle='solid', label='ω2')
ax2.plot(tiempos, ω3, color='#a500a9', linestyle='dotted', label='ω3', lw=2)
#ax2.axhline(ω3_min, color='gray', ls='--', lw=1)
#ax2.axhline(ω3_max, color='gray', ls='--', lw=1)
ax2.plot(tiempos, ω4, color='#00a5ba', linestyle='dotted', label='ω4', lw=2)
#ax2.axhline(ω4_min, color='gray', ls='--', lw=1)
#ax2.axhline(ω4_max, color='gray', ls='--', lw=1)
ax2.axhline(0, color='#303030', lw=1)
ax2.set_title("Velocidades Angulares")
ax2.set_ylabel("Velocidad [rad/s]")
ax2.set_xlabel("Tiempo [s]")
ax2.legend()
ax2.grid(True)



#Aceleraciones
# A_A + A_BA - A_B = 0

# A_A = (At_A + An_A)    = [r2*α2*(-Sen(θ2)+jCos(θ2) - r2*ω2^2*(Cos(θ2)+jSen(θ2))]
# A_BA = (At_AB + An_AB) = [r3*α3*(-Sen(θ3)+jCos(θ3) - r3*ω3^2*(Cos(θ3)+jSen(θ3))]
# A_B = (At_B + An_B)    = [r4*α4*(-Sen(θ4)+jCos(θ4) - r4*ω4^2*(Cos(θ4)+jSen(θ4))]

#Separando por componentes
#(x): - r2*α2*Sen(θ2) - r2*ω2^2*Cos(θ2) - r3*α3*Sen(θ3) - r3*ω3^2*Cos(θ3) + r4*α4*Sen(θ4) + r4*ω4^2*Cos(θ4)
#(y): - r2*α2*Cos(θ2) - r2*ω2^2*Sen(θ2) - r3*α3*Cos(θ3) - r3*ω3^2*Sen(θ3) - r4*α4*Cos(θ4) + r4*ω4^2*Sen(θ4)

#Funciones de Aceleraciones Angulares
def α3(β2, β3, β4, ώ3, ώ4):
    #Se definen las constantes:
    Aa = r4*math.sin(β4)
    Ba = r3*math.sin(β3)
    Ca = r2*α2*math.sin(β2) + r2*ω2**2*math.cos(β2) + r3*ώ3**2*math.cos(β3) - r4*ώ4**2*math.cos(β4)
    Da = r4*math.cos(β4)
    Ea = r3*math.cos(β3)
    Fa = r2*α2*math.cos(β2) - r2*ω2**2*math.sin(β2) - r3*ώ3**2*math.sin(β3) + r4*ώ4**2*math.sin(β4)
    return (Ca*Da - Aa*Fa)/(Aa*Ea - Ba*Da)

def α4(β2, β3, β4, ώ3, ώ4):
    #Se definen las constantes:
    Aa = r4*math.sin(β4)
    Ba = r3*math.sin(β3)
    Ca = r2*α2*math.sin(β2) + r2*ω2**2*math.cos(β2) + r3*ώ3**2*math.cos(β3) - r4*ώ4**2*math.cos(β4)
    Da = r4*math.cos(β4)
    Ea = r3*math.cos(β3)
    Fa = r2*α2*math.cos(β2) - r2*ω2**2*math.sin(β2) - r3*ώ3**2*math.sin(β3) + r4*ώ4**2*math.sin(β4)
    return (Ca*Ea - Ba*Fa)/(Aa*Ea - Ba*Da)

#Funciones de Aceleraciones Lineales
def A(r, θ, ω, α):
    At = r*α*np.array([-math.sin(θ), math.cos(θ)])
    An = -r*ω**2*np.array([math.cos(θ), math.sin(θ)])
    return np.array([At[0] + An[0], At[1] + An[1]])

#Evaluación
α3 = np.array(list(map(α3, θ2, θ3, θ4, ω3, ω4)))
α4 = np.array(list(map(α4, θ2, θ3, θ4, ω3, ω4)))

A2 = np.array(list(map(A,np.repeat(r2, length), θ2, np.repeat(ω2, length), np.repeat(α2, length))))
A3 = np.array(list(map(A,np.repeat(r3, length), θ3, ω3, α3)))
A4 = np.array(list(map(A,np.repeat(r4, length), θ4, ω4, α4)))


#Gráfica de Aceleraciones
fig3, ax3 = plt.subplots(figsize=(10, 8))

# Set axis ranges; by default this will put major ticks every 25.
ax3.set_xlim(0, 0.12)
ax3.set_ylim(-2500, 2000)

# Change major ticks to show every 20.
ax3.xaxis.set_major_locator(MultipleLocator(0.01))
ax3.yaxis.set_major_locator(MultipleLocator(250))

# Change minor ticks to show every 5. (20/4 = 5)
ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
ax3.yaxis.set_minor_locator(AutoMinorLocator(5))

# Turn grid on for both major and minor ticks and style minor slightly
# differently.
ax3.grid(which='major', color='#CCCCCC', linestyle='--')
ax3.grid(which='minor', color='#CCCCCC', linestyle=':')

ax3.plot(tiempos, α3, color='#00a980', linestyle='dotted', label='α3', lw=2)
#ax3.axhline(α3_min, color='gray', ls='--', lw=1)
#ax3.axhline(α3_max, color='gray', ls='--', lw=1)
ax3.plot(tiempos, α4, color='#8ca916', linestyle='dotted', label='α4', lw=2)
#ax3.axhline(α4_min, color='gray', ls='--', lw=1)
#ax3.axhline(α4_max, color='gray', ls='--', lw=1)
ax3.axhline(0, color='#303030', lw=1)
ax3.set_title("Aceleraciones Angulares")
ax3.set_ylabel("Aceleración [rad/s^2]")
ax3.set_xlabel("Tiempo [s]")
ax3.legend()
ax3.grid(True)

fig1.savefig("Figura1")
fig2.savefig("Figura2")
fig3.savefig("Figura3")


#ANÁLISIS DINÁMICO

g = 32.2 #[lb/s^2] Aceleración gravitacional
peso_combinado = 29 #[lb] Magnitud del peso combinado del peine y la barra de enlace
fuerza_externa = 540 #[lb] Magnitud de la FUerza de batido aplicada al peine
r5 = 3.75 #[pulg] Distancia desde la junta B hasta el Peine, donde se aplica la fuerza externa.

#La cadena cinemática de 4 barras está formada por eslabones con dimensiones:
ancho = 2 #[pulg]
grosor = 1 #[pulg]
sección_transversal = ancho*grosor #(pulg^2)
#Se asume que el material es Acero Inoxidable AISI 304, cuya densidad es:
#Densidad del Acero AISI 304: 8.027 g/cm^3
ρ = 8.027 * (2.2/1000)*(2.54**3/1) #(lb/pulg^3)

#Masa de los eslabones - se asume que una distribución uniforme de la masa en los eslabones.
m2 = r2*sección_transversal*ρ #[blobs]
m3 = r3*sección_transversal*ρ #[blobs]
m4 = r4*sección_transversal*ρ #[blobs]

#Calculo de los Momentos de Inercia
Ig2 = (m2/12)*(r2**2 + ancho**2)
Ig3 = (m3/12)*(r3**2 + ancho**2)
Ig4 = (m4/12)*(r4**2 + ancho**2)

#Se establece un Sistema Coordenado xy LNCS (Sistema de Coordenasdas rotatorio Local No insertado)
#en el centro de cada eslabón
for i in range(length):
    r12 = np.array([-(r2/2)*math.cos(θ2[i])        , -(r2/2)*math.sin(θ2[i])])
    # r12x = -(r2/2)*math.cos(θ2[i]) 
    # r12y = -(r2/2)*math.sin(θ2[i])
    # r12 = np.array([r12x, r12y])

    r32 = np.array([-(r2/2)*math.cos(θ2[i]+math.pi), -(r2/2)*math.sin(θ2[i]+math.pi)])
    # r32x = -(r2/2)*math.cos(θ2[i]+math.pi)
    # r32y = -(r2/2)*math.sin(θ2[i]+math.pi)
    # r32 = np.array([r32x, r32y])

    r23 = np.array([-(r3/2)*math.cos(θ3[i])        , -(r3/2)*math.sin(θ3[i])])
    # r23x = -(r3/2)*math.cos(θ3[i])
    # r23y = -(r3/2)*math.sin(θ3[i])
    # r23 = np.array([r23x, r23y])

    r43 = np.array([-(r3/2)*math.cos(θ3[i]-math.pi), -(r3/2)*math.sin(θ3[i]-math.pi)])
    # r43x = -(r3/2)*math.cos(θ3[i]-math.pi)
    # r43y = -(r3/2)*math.sin(θ3[i]-math.pi)
    # r43 = np.array([r43x, r43y])

    r34 = np.array([-(r4/2)*math.cos(θ4[i]+math.pi), -(r4/2)*math.sin(θ4[i]+math.pi)])
    # r34x = -(r4/2)*math.cos(θ4[i]+math.pi)
    # r34y = -(r4/2)*math.sin(θ4[i]+math.pi)
    # r34 = np.array([r34x, r34y])

    r14 = np.array([-(r4/2)*math.cos(θ4[i])        , -(r4/2)*math.sin(θ4[i])])
    # r14x = -(r4/2)*math.cos(θ4[i])
    # r14y = -(r4/2)*math.sin(θ4[i])
    # r14 = np.array([r14x, r14y])

    #r54 = np.array([-(r5+r4/2)*math.cos(θ4[i])     , -(r5+r4/2)*math.sin(θ4[i])])
    r54x = -(r5+r4/2)*math.cos(θ4[i])
    r54y = -(r5+r4/2)*math.sin(θ4[i])
    r54 = np.array([r54x, r54y])


#Cálculo de las componentes x y de la aceleración de los CG de todos los eslabones móviles
#en el Sistema Coordenado Global (CGS)
Ag2 = np.array(list(map(A,np.repeat(r2/2, length), θ2, np.repeat(ω2, length), np.repeat(α2,length))))
Ag3 = np.array(list(map(A,np.repeat(r3/2, length), θ3, ω3, α3)))
Ag4 = np.array(list(map(A,np.repeat(r4/2, length), θ4, ω4, α4)))


#Se separan los eslabones y se realiza un DCL para cada uno.

#ESLABÓN 2
#Fuerza de la Bancada sobre el Eslabón 2
# F12 = ( F12x, F12y)
# F12x = 0 #[N] (i)
# F12y = 0 #[N] (j)
F12 = np.array(2)
#Fuerza del Eslabón 3 sobre el Eslabón 2
# F32 = ( F32x, F32y)
# F32x = 0 #[N] (i)
# F32y = 0 #[N] (j)
F32 = np.array(2)
#Par de Torsión de entrada
T12 = 0 #[N*m] (k)
#Peso del Eslabón 2
W2 = m2*g #[N] (j)

#Eslabón 3
#Fuerza del Eslabón 2 sobre el Eslabón 3
# F23 = ( F32x, F32y)
# F23x = -F32x #[N] (i)
# F23y = -F32y #[N] (j)
F23 = np.array([-F32[0], -F32[1]])
#Fuerza del Eslabón 4 sobre el Eslabón 3
# F43 = ( F43x, F43y)
# F43x = 0 #[N] (i)
# F43y = 0 #[N] (j)
F43 = np.array(2)
#Peso del Eslabón 3
W3 = m3*g #[N] (j)

#ESLABÓN 4
#Fuerza del Eslabón 3 sobre el Eslabón 4
# # F34 = ( F34x, F34y)
# F34x = -F43x #[N] (i)
# F34y = -F43y #[N] (j)
F34 = np.array([-F43[0], -F43[1]])
#Fuerza de la Bancada sobre el Eslabón 4
# # F14 = ( F14x, F14y)
# F14x = 0 #[N] (i)
# F14y = 0 #[N] (j)
F14 = np.array(2)
#Peso del Eslabón 4
W4 = m4*g #[N] (j)
#Peso de la Barra de Enlace
Wc = peso_combinado #[N] (j)
#Fuerza Externa de Batido - se considera que siempre es horizontal
# Fb = ( Fbx, Fby)
Fb = np.array([-fuerza_externa, 0])
Fbx = Fb[0] #[N] (i)
Fby = Fb[1] #[N] (j)

#Se aplican las Ecuaciones de Newton-Euler

# Eslabón 2
#   F12x + F32x = m2*Ag2x
#   F12y + F32y - W2 = m2*Ag2y
#   T12 + (r12x*F12y - r12y*F12x) + (r32x*F32y - r32y*F32x) = Ig2*α2

# Eslabón 3
#   F43x - F32x = m3*Ag3x
#   F43y - F32y -W3 = m3*Ag3y
#   (r43x*F43y - r43y*F43x) - (r23x*F32y - r23y*F32x) = Ig3*α3

# Eslabón 4
#   F14x - F43x + Fbx = m4*Ag4x
#   F14y - F43y + Fby -W4 = m4*Ag4y
#   (r14x*F14y - r14y*F14x) - (r34x*F43y - r34y*F43x) + (r54x*Fby - r54y*Fbx) = Ig4*α4

#Desarrollando la Matriz

#      F12x    F12y    T12     F32x    F32y    F43x    F43y    F14x   F14y

#   |   1       0       0       1       0       0       0       0      0   |         |   F12x   |
#   |   0       1       0       0       1       0       0       0      0   |         |   F12y   |
#   | -r12y    r12x     1     -r32y    r32x     0       0       0      0   |         |   T12    |
#   |   0       0       0      -1       0       1       0       0      0   |         |   F32x   |
#   |   0       0       0       0      -1       0       1       0      0   |    X    |   F32y   |   =
#   |   0       0       0      r23y   -r23x   -r43y    r43x     0      0   |         |   F43x   |
#   |   0       0       0       0       0      -1       0       1      0   |         |   F43y   |
#   |   0       0       0       0       0       0      -1       0      1   |         |   F14x   |
#   |   0       0       0       0       0      r34y   -r34x   -r14y   r14x |         |   F14y   |

#       |   m2*Ag2x   |
#       |  m2*Ag2y+W2 |
#       |   Ig2*α2    |
#       |   m3*Ag3x   |
#   =   |  m3*Ag3y+W3 |
#       |   Ig3*α3    |
#       | m4*Ag4x-Fbx |
#       |m4*Ag4y-Fby+W4|
#  |Ig4*α4-r54x*Fby+r54y*Fbx|

def matriz(r12, r32, r23, r43, r34, r14, F12, T12, F32, F43, F14, Ag2, Ag3, Ag4 ):
    r12x, r12y = r12
    r32x, r32y = r32
    r23x, r23y = r23
    r43x, r43y = r43
    r34x, r34y = r34
    r14x, r14y = r14
    F12x, F12y = F12
    F32x, F32y = F32
    F43x, F43y = F43
    F14x, F14y = F14
    Ag2x, Ag2y = Ag2
    Ag3x, Ag3y = Ag3
    Ag4x, Ag4y = Ag4

    M = np.array([
    [   1,      0,      0,      1,      0,      0,      0,      0,     0,  ],
    [   0,      1,      0,      0,      1,      0,      0,      0,     0,  ],
    [ -r12y,   r12x,    1,    -r32y,   r32x,    0,      0,      0,     0,  ],
    [   0,      0,      0,     -1,      0,      1,      0,      0,     0,  ],
    [   0,      0,      0,      0,     -1,      0,      1,      0,     0,  ],
    [   0,      0,      0,     r23y,  -r23x,  -r43y,   r43x,    0,     0,  ],
    [   0,      0,      0,      0,      0,     -1,      0,      1,     0,  ],
    [   0,      0,      0,      0,      0,      0,     -1,      0,     1,  ],
    [   0,      0,      0,      0,      0,     r34y,  -r34x,  -r14y,  r14x ]
    ])

    Incognitas = np.array([
        F12x,
        F12y,
        T12,
        F32x,
        F32y,
        F43x,
        F43y,
        F14x,
        F14y
    ])

    Soluciones = np.array([
        m2*Ag2x,
        m2*Ag2y+W2,
        Ig2*α2,   
        m3*Ag3x,  
        m3*Ag3y+W3,  
        Ig3*α3,
        m4*Ag4x-Fbx,
        m4*Ag4y-Fby+W4,
        Ig4*α4-r54x*Fby+r54y*Fbx
    ])


plt.tight_layout()
plt.show()