import math
import numpy as np
# Distancia del cojinete 1 hacia el Peso 2
l12 = 1.75 # [pulg]
# Distancia del cojinete 4 hacia el Peso 3
l43 = 1.75 # [pulg]
# Distancia del Cojinete 1 hacia el Cojinete 4
l14 = 10.0 # [pulg]

# Distancia del Peso 3 hacia el Peso 2
l32 = l14 - l12 - l43 # [pulg]
# Distancia del cojinete 4 hacia el Peso 2
l42 = l14 - l12 # [pulg]
g = 39.37 # [pulg/s^2] 

print("Ingrese las Variables de Entrada")
rpm = float(input('Velocidad de operación (rpm): '))
print("\nFuerza Desbalanceada 1: ")
f1 = float(input('   magnitud (lb): '))
θ1 = float(input('   ángulo (°): '))
print("\nFuerza Desbalanceada 2: ")
f2 = float(input('   magnitud (lb): '))
θ2 = float(input('   ángulo (°): '))
Y = float(input('\nPeso de Balanceo (lb): \n'))

m2 = m3 = Y/g
e2 = 0
e3 = 0
w = rpm*(2*math.pi/60)

m3_xe3_x = (-(f1 * math.cos(θ1)*l12) - (f2*math.cos(θ2)*l42)) / (l32*w**2)
m3_ye3_y = (-(f1 * math.sin(θ1)*l12) - (f2*math.sin(θ2)*l42)) / (l32*w**2)
θ3 = math.atan(m3_ye3_y/m3_xe3_x)*(180/math.pi)
m3e3 = math.sqrt((m3_xe3_x**2)+(m3_ye3_y**2))
e3 = m3e3/m3
print('θ3: ',θ3,'°')
print('e3: ',e3,'pulg')

m2_xe2_x = (-f1 * math.cos(θ1) - f2*math.cos(θ2) - m3_xe3_x*w**2) / w**2
m2_ye2_y = (-f1 * math.sin(θ1) - f2*math.sin(θ2) - m3_ye3_y*w**2) / w**2
θ2 = math.atan(m2_ye2_y/m2_xe2_x)*(180/math.pi)
m2e2 = math.sqrt((m2_xe2_x**2)+(m2_ye2_y**2))
e2 = m2e2/m2
print('θ2: ',θ2,'°')
print('e2: ',e2,'pulg')