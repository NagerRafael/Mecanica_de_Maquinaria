import numpy as np
import math as mt

g = 32.2 #(ft/s^2)

#Valores del Sistema
r_O2A = 3 #(in) distancia desde el Eje de la Biela (O2) hasta la Junta A

ang_AB = 79.82
r_AB  = 12#(in) distancia desde la junta A hasta la junta B en la corredera
R_AB = np.array([r_AB*mt.sin(ang_AB), -r_AB*mt.cos(ang_AB)])

r_Ag3 =4.5#(in) distancia desde la junta A hasta el centro de masa del Eslabón 3
R_Ag3 = np.array([r_Ag3*mt.sin(ang_AB), -r_Ag3*mt.cos(ang_AB)])

W3 = 3.4 #(lb) peso del Eslabón 3
W4 = 2.86#(lb) peso del Eslabón 4

Ig2 = 0.352 #(in*Lb*s^2) Momento de Inercia del Eslabón 2
Ig3 = 0.108 #(in*Lb*s^3) Momento de Inercia del Eslabón 3

w2 = np.array([0, 0, +210]) #(rad/s)[k] Vector de velocidad Ángular del Eslabón 2

α2 = np.array([0, 0, 0]) #(rad/s^2)[k] Vector de Aceleracion Ángular del Eslabón 2
α3 = np.array([0, 0, -7670]) #(rad/s^3)[k] Vector de Aceleracion Ángular del Eslabón 3

Ace3 = np.array([-7820, +4876, 0]) #(ft/s^2) [i j] Aceleración Lineal del Eslabón 3
Ace4 = np.array([-7850, 0, 0]) #(ft/s^2) [i] Aceleracion Lineal de la Corredera

#DCL Corredera
# SUM Fx = m4 * Ace4x
# Bx = m4 * Ace4x
Bx = (W4/g)*Ace4[0]

# SUM Fy = m4 * Ace4y
# +By - F14(j) -W4 = m4 * Ace4y
#By = F14 + W4

#DCL Eslabón 3
# SUM Fx = m3 * Ace3x
# -Bx +Ax = m3 * Ace3x
Ax = (W3/g)*Ace3[0] + Bx

# SUM Fy = m3 * Ace3y
# +Ay -By -W3 = m3 * Ace3y

#Aplicando el Teorema de Ejes Paralelos
# SUM M_A = Icm3*α3 + r_cm3/A*m3*Ace3
# R_AB[1]*Bx + R_AB[0]*By + R_Ag3[0]*W3 = Icm3*α3 + R_Ag3*(m3*Ace3)

A = np.array([2,4,6])
B = np.array([1,3,5])
print(np.cross(A, B))
