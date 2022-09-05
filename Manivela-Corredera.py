import numpy as np
import math

# Libro de Robert Norton - pág 132
# SOLUCIÓN DE POSICIÓN DE UN MECANISMO
# DE CUATRO BARRAS DE MANIVELA-CORREDERA
# El caso general es un
# mecanismo manivela-corredera de cuatro barras descentrado
# El término descentrado significa que el eje de la corredera extendido no pasa por el pivote de la manivela.

def ManivelaCorredera_Pos(a,b,c,θ2):

    # Se definen cuatro vectores, R1, R2, R3 y R4
    # con R1 dispuesto paralelamente al eje del deslizamiento y R4 perpendicular.

    # R1 y R4 son componentes ortogonales del vector Rs que va desde O2 hasta B
    # R1 representa la posición de la corredera con magnitud d
    # R4 es ortogonal a R1 y define el descentrado de magnitud constante del mecanismo.
    # R3 es el vector posición del acoplador con raíz en la corredera, punto B

    # Se definen los vectores que permitan medir los ángulos de los eslabones
        # R1: ( d ; θ1=0° ) con pivote en O2
        # R2: ( a ; θ2) con pivote en O2
        # R3: ( b ; θ3) con pivote en B
        # R4: ( c ; θ4) con pivote en O'4 (es un punto de referencia, mas no es un eje de rotación)

    # De este modo se puede definir la siguiente ecuación de Lazo Vectorial:
    # R2 = R1 + R4 + R3
    # R2 - R3 - R4 - R1 = 0

    # Reescribiendo en forma polar :
    # a*e^(jθ2) - b*e^(jθ3) - c*e^(jθ4) - d*e^(jθ1) = 0

    # Al sustituir por las equivalentes de Euler:
    # a*(Cos(θ2) + jSen(θ2)) - b*(Cos(θ3) + jSen(θ3)) - c*(Cos(θ4) + jSen(θ4)) - d*(Cos(θ1) + jSen(θ1)) = 0

    # Separando la parte Real de la Imaginaria
    # Real: a*Cos(θ2) - b*Cos(θ3) - c*Cos(θ4) - d*Cos(θ1) = 0
    # Imaginario: a*Sin(θ2) - b*Sin(θ3) - c*Sin(θ4) - d*Sin(θ1) = 0

    # Sabiendo que θ1 = 0°
    # a*Cos(θ2) - b*Cos(θ3) - c*Cos(θ4) - d = 0
    # a*Sin(θ2) - b*Sin(θ3) - c*Sin(θ4) = 0

    # Las Incognitas son: [ d , θ3]
    # La longitud del eslabón d y el ángulo del eslabón θ3

    # La Variable Independiente es: [ θ2 ]
    # El ángulo de la Manivela θ2

    # Los Valores conocidos son: [ a, b, c, θ4 ]
    # Las longitudes de los dos eslabones, la distancia de Descentrado y el ángulo θ4

    # Debido al Sistema de Coordenadas paralelo y perpendicular al Bloque deslizante:
    # θ1 = 0° y θ4 = 90°

    # Reemplazando en las ecuaciones anteriores
    # d = a*Cos(θ2) - b*Cos(θ3)
    # θ3 = arcSin( (a*Sin(θ2) - c) / b )

    # Observese que si:
    # [ a_y < c ] -> θ3 se encuentra en el tercer o cuarto cuadrante
    # [ a_y > c ] -> θ3 se encuentra en el primer o segundo cuadrante

    # Obsérvese que de nuevo existen dos soluciones válidas correspondientes a los dos circuitos del
    # mecanismo. La función arcoseno es de valores múltiples.
    # Su evaluación dará un valor entre ±90°, que representa sólo un circuito del mecanismo.

    # Suponiendo que se evalua dentro del arcSin(x) un valor negativo [-1<x<0]
    # Se obtienen 2 ángulos cuyo seno, Sen(θ3_1)  Sen(θ3_2)es el mismo valor negativo (x)
    θ3_1 = math.asin( (a*math.sin(θ2) - c) / b ) # Python entrega el ángulo entre 0° y -90°
    # Para calcular el segundo ángulo:
    θ3_2 = math.asin( -(a*math.sin(θ2) - c) / b ) + math.pi
    # Al multiplicar para -1 el argumento de arcSen() se obtiene los ángulos sumplementarios
    # En este caso python entrega el ánuglo entre 0° y 90°, al cual se le suma π (180°) y se obtiene el segundo ángulo buscado

    d_1 = a*math.cos(θ2) - b*math.cos(θ3_1)
    d_2 = a*math.cos(θ2) - b*math.cos(θ3_2)

    return ((d_1, d_2), (θ3_1, θ3_2))

def radToDeg(rad): #Función para convertir de rad a grados hexadecimales
    return rad*(180.0/math.pi)


# a = 40 # [mm]
# b = 120 #[mm]
# c = -20 #[mm]
# θ2 = 60*(math.pi/180) # [rad]
# print('Ejemplo de prueba 4-2 (Libro de Robert Norton - pág 134)')
# print('a: {:>8.3f}\nb: {:>8.3f}\nc: {:>8.3f}\nθ2: {:>7.3f}'.format(a,b,c,θ2))

# lvc = ManivelaCorredera_Pos(a,b,c,θ2)

# print("\nLVC",lvc)
# print("\nConfiguración: [ Abierta ; Cruzada ]")
# print("d: [{:6.3f} mm ; {:6.3f} mm]\nθ3: [{:6.3f}° ; {:6.3f}°]".format(lvc[0][0], lvc[0][1], radToDeg(lvc[1][0]), radToDeg(lvc[1][1])))


# a = 3 # [in]
# b = 12 #[in]
# c = 0 #[in]
# θ2 = 45*(math.pi/180) # [rad]
# print('Ejercicio en CLase 14 (Livingston)')
# print('a: {:>8.3f}\nb: {:>8.3f}\nc: {:>8.3f}\nθ2: {:>7.3f}'.format(a,b,c,θ2))

# lvc = ManivelaCorredera_Pos(a,b,c,θ2)

# print("\nLVC",lvc)
# print("\nConfiguración: [ Abierta ; Cruzada ]")
# print("d: [{:6.3f} mm ; {:6.3f} mm]\nθ3: [{:6.3f}° ; {:6.3f}°]".format(lvc[0][0], lvc[0][1], radToDeg(lvc[1][0]), radToDeg(lvc[1][1])))