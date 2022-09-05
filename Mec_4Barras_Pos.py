import numpy as np
import math

n = 0

#VALORES DEFINIDOS - Longitud de los Eslabones

# a = float(input('Longitud de la Manivela (eslabón a): '))
# b = float(input('Longitud de la Biela (eslabón b): '))
# c = float(input('Longitud del Balancí (eslabón c): '))
# d = float(input('Distancia de las Bancadas (eslabón d): '))

# print('\n{:_^25s}\n'.format('VALORES INGRESADOS'))
# for X in [a,b,c,d]:
#     n+=1
#     print("Eslabón {:d} = {:6.3f} [pulg]".format(n,X))

# n = 0
# print('\n{:_^25s}\n'.format('MODOS DE EJECUCIÓN'))
# for op in ['Evaluar un solo valor de θ',
#            'Evaluar un conjunto de valores de θ',
#            'Evaluar un rango de valores de θ',
#            'Evaluar una rotación de 360°']:
#     n+=1
#     print('{:d}) {:s}'.format(n,op))
# md = int(input('<Ingresar Opción>: '))

# MÉTODO ALGEBRAICO
def Pos_Algeb(a, b, c, d, θ2):
    # Coordenadas de la Junta A:
    Ax = a*math.cos(θ2)
    Ay = a*math.sin(θ2)

    # Las coordenadas del punto B se encuentran con las ecuaciones de los círculos en torno a A y O4.
    # b^2 = (Bx - Ax)^2 + (By - Ay)^2
    # c^2 = (Bx - d)^2 + By^2

    # Si se resta b^2 - c^2 se obtiene una ecuación para Bx
    # Bx = (a2 - b2 + c^2 - d^2) / 2*(Ax - d) - 2*Ay*By / 2*(Ax - d)
    # sea: S = (a2 - b2 + c^2 - d^2) / 2*(Ax - d)
    # Bx = S - 2*Ay*By / 2*(Ax - d)

    # Si se sustituye Bx en c^2:
    # By^2 + [S - 2*Ay*By / 2*(Ax - d) - d]^2 - c^2 = 0

    # Se expande la ecuación para agruparla con By
    # By^2 + [ S^2 - 2*S*Ay*By/(Ax-d) - 2*S*d + Ay^2*By^2/(Ax-d) + 2*d*Ay*By/(Ax-d) + d^2 ] - c^2 = 0
    # (S^2 - 2*S*d + d^2) = (S - d)^2
    # By^2*(1 + Ay^2/(Ax-d)) + By*(2*Ay*(d-S)/(Ax-d)) + ((S - d)^2 - c^2)

    # Se definen las siguentes constantes:
    S = (a**2 - b**2 + c**2 - d**2) / 2*(Ax - d)
    P = ((Ay**2) / ((Ax - d)**2)) + 1
    Q = (2*Ay*(d - S)) / (Ax - d)
    R = (S - d)**2 - c**2

    # Esta es una cuadrática en función de By^2, se resulven sus raices
    # By = (-Q +- sqrt(Q^2 - 4PR)) / 2*P

    # Se obtendrá una respuesta para el Circuito Abierto y otra para el Circuito Cruzado
    # calculating  the discriminant
    dis = (Q**2) - (4*P*R)
    print(dis)
    if dis >= 0:
        # find two results
        By_c = (-Q - math.sqrt(abs(dis)))/(2*P) # Circuito Cruzado
        By_a = (-Q + math.sqrt(abs(dis)))/(2*P) # Circuito Abierto

        # Se puede determinar Bx
        Bx_c = S - 2*Ay*By_c / 2*(Ax - d)
        Bx_a = S - 2*Ay*By_a / 2*(Ax - d)
        
        # Recordando el triángulo pitagórico de b y c:
        # b^2 = (Bx - Ax)^2 + (By - Ay)^2
        # c^2 = (Bx - d)^2 + By^2
        # Los ángulos de los eslabones 3 y 4 para esta posición son:
        θ3_a = math.atan((By_a - Ay) / (Bx_a - Ax))
        θ3_c = math.atan((By_c - Ay) / (Bx_c - Ax))
        
        θ4_a = math.atan(By_a / (Bx_a - d))
        θ4_c = math.atan(By_c / (Bx_c - d))
        
        return((θ3_a, θ3_c), (θ4_a, θ4_c))

    else:
        return print("<raices complejas conjugadas>\nLa Cadena Cinemática está Abierta")


# MÉTODO DE lAZO VECTORIAL CERRADO (Raven) - (Libro de Robert Norton - pág 127)
def Pos_LVC(a, b, c, d, θ2):
    # Se definen los vectores que permitan medir los ángulos de los eslabones
    # R1: ( d ; θ1=0° ) con pivote en O2
    # R2: ( a ; θ2) con pivote en O2
    # R3: ( b ; θ3) con pivote en A
    # R4: ( c ; θ4) con pivote en O4

    # R2 + R3 = R1 + R4 -> R2 + R3 - R4 - R1 = 0 

    # Reescribiendo en forma polar :
    # a*e^(jθ2) + b*e^(jθ3) - c*e^(jθ4) - d*e^(jθ1) = 0

    # Al sustituir por las equivalentes de Euler:
    # a*(Cos(θ2) + jSen(θ2)) + b*(Cos(θ3) + jSen(θ3)) - c*(Cos(θ4) + jSen(θ4)) - d*(Cos(θ1) + jSen(θ1)) = 0

    # Separando la parte Real de la Imaginaria
    # Real: a*Cos(θ2) + b*Cos(θ3) - c*Cos(θ4) - d*Cos(θ1) = 0
    # Imaginario: a*Sin(θ2) + b*Sin(θ3) - c*Sin(θ4) - d*Sin(θ1) = 0

    # Sabiendo que θ1 = 0°
    # a*Cos(θ2) + b*Cos(θ3) - c*Cos(θ4) - d = 0
    # a*Sin(θ2) + b*Sin(θ3) - c*Sin(θ4) = 0

    # Se tiene dos incognitas, se despeja una e funcion de la otra
    # b*Cos(θ3) = - a*Cos(θ2) + c*Cos(θ4) + d
    # b*Sin(θ3) = - a*Sin(θ2) + c*Sin(θ4) 

    # Se elevan al cuadrado las dos ecuaciones y se suman
    # [ b*Cos(θ3) ]^2 = [ - a*Cos(θ2) + c*Cos(θ4) + d ]^2
    # [ b*Sin(θ3) ]^2 = [ - a*Sin(θ2) + c*Sin(θ4) ]^2
    # b^2*(Cos^2(θ3)+Sin^2(θ3)) = [-a*Cos(θ2)+c*Cos(θ4)+d]^2 + [-a*Sin(θ2)+c*Sin(θ4)]^2

    # Se expande la ecuación y se asocian terminos similares
    # b^2 = a^2 + c^2 + d^2 - 2*a*d*Cos(θ2) + 2*c*d*Cos(θ4) -2*a*c*[Cos(θ2)*Cos(θ4) + Sin(θ2)*Sin(θ4)]

    # Se divide para 2*a*c y se reordena
    # (d/a)*Cos(θ4) - (d/c)*Cos(θ2) + (a^2-b^2+c^2+d^2)/(2*a*c) = Cos(θ2)*Cos(θ4) + Sin(θ2)*Sin(θ4)

    # Se definen las siguientes constantes:
    k1 = d/a
    k2 = d/c
    k3 = (a**2-b**2+c**2+d**2)/(2*a*c)

    # Se obtiene:
    # k1*Cos(θ4) - k2*Cos(θ2) + k3 = Cos(θ2)*Cos(θ4) + Sin(θ2)*Sin(θ4)

    # /Si se sustituye Cos(θ2-θ4) = Cos(θ2)*Cos(θ4) + Sin(θ2)*Sin(θ4)
    # ECUACIÓN DE FREUDENSTEIN
    # k1*Cos(θ4) - k2*Cos(θ2) + k3 = Cos(θ2-θ4)
    # /

    # Se aplican las Identidades Semiangulares
    # para convertir los terminos Sen(θ4) y Cos(θ4) en Tan(θ4)
    # Sin(θ4) = 2*Tan(θ4/2) / (1 + Tan^2(θ4/2))
    # Cos(θ4) = (1 - Tan^2(θ4/2)) / (1 + Tan^2(θ4/2))

    # Definiendo las siguientes constantes:
    A = math.cos(θ2) - k1 - k2*math.cos(θ2) + k3
    B = -2*math.sin(θ2)
    C = k1 - (k2+1)*math.cos(θ2) + k3

    # De este modo de sedifine la siguiente ecuación cuadrática:
    # A*Tan^2(θ4/2) + B*Tan^2(θ4/2) + C = 0

    # Las rices de esta ecuación son:
    # Tan(θ4/2) = (-B +- sqrt(B^2 - 4*A*C)) / 2*A

    # Donde se puede despejar θ4:
    # θ4_1,2 = 2*arcTan[(-B +- sqrt(B^2 - 4*A*C)) / 2*A] 

    # calculating  the discriminant
    dis = (B**2) - (4 * A*C)
    
    if dis >= 0:
        # find two results
        ans1 = (-B - math.sqrt(abs(dis)))/(2*A)
        ans2 = (-B + math.sqrt(abs(dis)))/(2*A)

        #En el mecanismo de cuatro barras, la solución negativa da θ4 para la configuración
        # abierta, y la positiva da θ4 para la confi guración cruzada.

        #Ángulo θ4 en Configuración Abierta:
        θ4_a = 2*math.atan(ans1)
        
        #Ángulo θ4 en Configuración Cruzada:
        θ4_c = 2*math.atan(ans2)

        # Los términos cruzado y abierto están basados en la suposición
        # de que el eslabón de entrada 2, para el cual θ2 está definido, se encuentra en el primer cuadrante
        # (es decir, 0 < θ2 < π/2). Un mecanismo de Grashof se define entonces como cruzado si los
        # dos eslabones adyacentes al eslabón más corto se cruzan entre sí, y como abierto si no lo hacen en
        # esta posición.

        # Para encontrar el ángulo θ3 se recuerda las ecuaciones:
        # b*Cos(θ3) = - a*Cos(θ2) + c*Cos(θ4) + d
        # b*Sin(θ2) = - a*Sin(θ2) + c*Sin(θ4) 
    
        # Se despeja para θ4
        # c*Cos(θ4) = - a*Cos(θ2) + b*Cos(θ3) + d
        # c*Sin(θ4) = - a*Sin(θ2) + b*Sin(θ3) 
        
        # Se elevan al cuadrado las dos ecuaciones y se suman
        # [ c*Cos(θ4) ]^2 = [ - a*Cos(θ2) + b*Cos(θ3) + d ]^2
        # [ c*Sin(θ4) ]^2 = [ - a*Sin(θ2) + b*Sin(θ3) ]^2
        # c^2*[Cos^2(θ4)+Sin^2(θ4)] = [-a*Cos(θ2)+b*Cos(θ3)+d]^2 + [-a*Sin(θ2)+b*Sin(θ3)]^2

        # Se expande la ecuación y se asocian terminos similares
        # c^2 = a^2 + b^2 + d^2 - 2*a*d*Cos(θ2) + 2*b*d*Cos(θ3) -2*a*b*[Cos(θ2)*Cos(θ3) + Sin(θ2)*Sin(θ3)]

        # Se divide para 2*a*b y se reordena
        # (d/a)*Cos(θ3) - (d/b)*Cos(θ2) + (a^2+b^2-c^2+d^2)/(2*a*b) = Cos(θ2)*Cos(θ3) + Sin(θ2)*Sin(θ3)

        # Se definen las siguientes constantes:
        k4 = d/b
        k5 = (-a**2-b**2+c**2-d**2)/(2*a*b)

        # Se obtiene:
        # k1*Cos(θ3) + k4*Cos(θ2) + k5 = Cos(θ2)*Cos(θ3) + Sin(θ2)*Sin(θ3)

        # Se aplican las Identidades Semiangulares

        # Definiendo las siguientes constantes:
        D = math.cos(θ2) - k1 + k4*math.cos(θ2) + k5
        E = -2*math.sin(θ2)
        F = k1 + (k4-1)*math.cos(θ2) + k5

        # De este modo de se define la siguiente ecuación cuadrática:
        # D*Tan^2(θ4/2) + E*Tan^2(θ4/2) + F = 0

        # Las rices de esta ecuación son:
        # Tan(θ3/2) = (-E +- sqrt(E^2 - 4*D*F)) / 2*D

        # Donde se puede despejar θ3:
        # θ3_1,2 = 2*arcTan[(-E +- sqrt(E^2 - 4*D*F)) / 2*D] 

        # calculating  the discriminant
        dis = (E**2) - (4 * D*F)
    
        if dis >= 0:
            # find two results
            ans1 = (-E - math.sqrt(abs(dis)))/(2*D)
            ans2 = (-E + math.sqrt(abs(dis)))/(2*D)

            #Ángulo θ3 en Configuración Abierta:
            θ3_a = 2*math.atan(ans1)
            
            #Ángulo θ3 en Configuración Cruzada:
            θ3_c = 2*math.atan(ans2)

            #return(math.pi - θ4_a, θ4_c + math.pi)
            # print("{:6.3f} {:6.3f} {:6.3f}".format(A,B,C))
            return((θ3_a, θ3_c), (θ4_a, θ4_c))

        else:
            print("<--raices conjugadas compleja-->\nLos Eslabones no se Conectan\no\nMecanismo de No Grashof fuera de la posición límite estacionaria")

    else:
        return print("<--raices conjugadas compleja-->\nLos Eslabones no se Conectan\no\nMecanismo de No Grashof fuera de la posición límite estacionaria")

def radToDeg(rad): #Función para convertir de rad a grados hexadecimales
    return rad*(180.0/math.pi)


# Función de Velocidad Angular
def ω(a, b, c, β2, β3, β4, ω2):
    ω3 = (a*ω2/b)*(math.sin(β4-β2)/math.sin(β3-β4))
    ω4 = (a*ω2/c)*(math.sin(β2-β3)/math.sin(β4-β3))
    return (ω3, ω4)

#Función de Aceleración Angular
def α(a, b, c, β2, β3, β4, ώ2, ώ3, ώ4, α2):
    #Se definen las constantes:
    A = c*math.sin(β4)
    B = b*math.sin(β3)
    C = a*α2*math.sin(β2) + a*ώ2**2*math.cos(β2) + b*ώ3**2*math.cos(β3) - c*ώ4**2*math.cos(β4)
    D = c*math.cos(β4)
    E = b*math.cos(β3)
    F = a*α2*math.cos(β2) - a*ώ2**2*math.sin(β2) - b*ώ3**2*math.sin(β3) + c*ώ4**2*math.sin(β4)
    α3 = (C*D - A*F) / (A*E - B*D)
    α4 = (C*E - B*F) / (A*E - B*D)
    return (α3, α4)

#Función de la Aceleración Lineal
def A(r, θ, ω, α):
    At = r*α*np.array([-math.sin(θ), math.cos(θ)])
    An = -r*ω**2*np.array([math.cos(θ), math.sin(θ)])
    return np.array([At[0] + An[0], At[1] + An[1]])


# Pruebas:

# a = 40 #[mm]
# b = 120 #[mm]
# c = 80 #[mm]
# d = 100 #[mm]
# θ2 = 40*(math.pi/180) #[rad]
# print('Ejemplo de prueba 4-1 (Libro de Robert Norton - pág 131)')
# print('a: {:>8.3f}\nb: {:>8.3f}\nc: {:>8.3f}\nd: {:>8.3f}\nθ2: {:>7.3f}'.format(a,b,c,d,θ2))
# # alg = Pos_Algeb(a,b,c,d,θ2)
# lvc = Pos_LVC(a,b,c,d,θ2)

# # print("\nALG",alg)
# # print("\nConfiguración: [ Abierta ; Cruzada ]")
# # print("θ3: [{:6.3f} ; {:6.3f}]\nθ4: [{:6.3f} ; {:6.3f}]".format(radToDeg(alg[0][0]), radToDeg(alg[0][1]), radToDeg(alg[1][0]), radToDeg(alg[1][1])))

# print("\nLVC",lvc)
# print("\nConfiguración: [ Abierta ; Cruzada ]")
# print("θ3: [{:6.3f} ; {:6.3f}]\nθ4: [{:6.3f} ; {:6.3f}]".format(radToDeg(lvc[0][0]), radToDeg(lvc[0][1]), radToDeg(lvc[1][0]), radToDeg(lvc[1][1])))


# a = 0.2 #[m]
# b = 0.3 #[m]
# c = 0.2 #[m]
# d = 0.3 #[mm]
# θ2 = 135*(math.pi/180) #[rad]
# print('Ejemplo de Lección (Livingston)')
# print('a: {:>8.3f}\nb: {:>8.3f}\nc: {:>8.3f}\nd: {:>8.3f}\nθ2: {:>7.3f}'.format(a,b,c,d,θ2))

# lvc = Pos_LVC(a,b,c,d,θ2)

# print("\nLVC",lvc)
# print("\nConfiguración: [ Abierta ; Cruzada ]")
# print("θ3: [{:6.3f} ; {:6.3f}]\nθ4: [{:6.3f} ; {:6.3f}]".format(radToDeg(lvc[0][0]), radToDeg(lvc[0][1]), radToDeg(lvc[1][0]), radToDeg(lvc[1][1])))



# a = 0.2 # [m]
# b = 0.3 #[m]
# c = 0.2 #[m]
# d = 0.3 #[m]
# θ2 = 135*(math.pi/180) # [rad]
# α2 = -1 # [rad/s^2]
# ω2 = -2 # [rad/s]
# print('Lección 3 (Livingston)')
# print('a: {:>8.3f}\nb: {:>8.3f}\nc: {:>8.3f}\nθ2: {:>7.3f}'.format(a,b,c,θ2))

# lvc = Pos_LVC(a,b,c,d,θ2)
# θ3 = lvc[0][0]
# θ4 = lvc[1][0]
# print("\nLVC",lvc)
# print("\nConfiguración: [ Abierta ; Cruzada ]")
# print("d: [{:6.3f} mm ; {:6.3f} mm]\nθ3: [{:6.3f}° ; {:6.3f}°]".format(lvc[0][0], lvc[0][1], radToDeg(lvc[1][0]), radToDeg(lvc[1][1])))

# ω3, ω4 = ω(a, b, c, 135, 0, 135, -2)
# α3, α4 = α(a, b, c, 135, 0, 135, -2, ω3, ω4, -1)
# print("\nω",ω3,ω4)
# print("ω3: {:6.3f} [rad/s]\nω4: {:6.3f} [rad/s]".format(ω3, ω4))

# print("\nα",α3,α4)
# print("α3: {:6.3f} [rad/s^2]\nα4: {:6.3f} [rad/s^2]".format(α3, α4))

# A2 = A(0, θ2, ω2, α2)
# print("\nA",A2)
# print("A2_x: {:6.3f} [m/s^2]\nA2_y: {:6.3f} [m/s^2]\n".format(A2[0], A2[1]))
# A3 = A(b, θ3, ω3, α3)
# print("\nA3",A3)
# print("A3_x: {:6.3f} [m/s^2]\nA3_y: {:6.3f} [m/s^2]\n".format(A3[0], A3[1]))
# A4 = A(0, θ4, ω4, α4)
# print("\nA4",A4)
# print("A4_x: {:6.3f} [m/s^2]\nA4_y: {:6.3f} [m/s^2]\n".format(A4[0], A4[1]))

print(A(3, 30*np.pi/180, 25, -40))