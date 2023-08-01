'''

MONTECARLO del crecimiento de una interfase

Una de las fases está formada por células, la otra corresponde a la superficie donde se adhieren (ej, poliestireno de la cápsula de
petri) y medio de cultivo.

El modelo es SIMPLE.

Cada espacio puede estar vacío u ocupado por una o más células.

Las células se comportan como adherentes epiteliales, forman monocapas.

Es decir que cada espacio corresponde a un tamaño definido (ej 30um x 30um), y el tamaño de las células en ese espacio depende
del número de células (si hay una es de 30umx30um, si hay 4 son de 15umx15um, y así).

Las células pueden moverse y duplicarse.

EN PRINCIPIO, las células se moveran a espacios adyacentes (es decir, una célula puede moverse a algunos de los cuadrados que la rodean, y
no a alguno que está a dos de distancia). Esto podría llegar a ajustarse por si se quieren estudiar otros parámetros como velocidad
de las células y direccionalidad.

SIEMPRE se moverán con mayor probabilidad desde lugares de mayor DENSIDAD a lugares de menor densidad.

El grado de avance del frente va a estar relacionado con el número de movimientos y de duplicaciones por cada paso de montecarlo, y
eso debería guardar una relación con los parámetros experimentales.

Ej: un paso de montecarlo son 24 hs. (((TENER EN CUENTA QUE DEPENDE DEL factor_intentos_por_paso)))
Entonces, si la línea celular tiene un tiempo de duplicación de un día, debería ser muy probable que todas las células se hayan
duplicado en un paso de montecarlo.
Entonces, si el avance de frente corresponde con 120 um, la simulación debería haber avanzado 4 casilleros en promedio.

El número de eventos dentro de un PASO DE MONTECARLO se ajustará con el número de células en cultivo.

Es de esperar que a densidades celulares altas, las probabilidades de duplicación sean muy bajas.-> podría pensarse para calcular una probabilidad de duplicación


'''

# import math #importa el módulo math para algunas funciones

from math import sqrt, exp #importamos la función raiz cuadrada del módulo math

import random #importa el módulo random para los números aleatorios

import numpy #importa el módulo numpy para manejar los arrays,calcular el gradiente y módulos

from time import time

import matplotlib.pyplot as plt

# import os

### DATOS DE LA MATRIZ INICIAL (FRENTE)

numero_pasos = 1000 # establece el número de pasos de Montecarlo que hará la simulación

densidad_maxima = 8 # establece un máximo de ocupación de un lugar. Debería relacionarse con los datos experimentales

numero_columnas = 64 # es el numero de columnas (largo, l) de la interfase (o del frente)

numero_filas = 30 #es el numero de filas (alto, h) de la interfase (perpendicular al frente)

pendiente = 0 #es la pendiente de la interface, la tg del ángulo entre el flujo y el frente. 0 da un frente normal. 1 da un frente a 45Âº. Si se quiere aumentar mucho la pendiente, recordar aumentar el nÃºmero de filas.

filas_iniciales_del_frente = 20 # indica el número incial de filas ocupadas por células en el frente, debería corresponder al menos a "frente efectivo"

#ya que cuando comience a avanzar el frente

densidad_inicial = 5 # indica la densidad inicial de células en el cultivo /// PODRíA LLEGAR A ALEATORIZARSE UN POCO///

movilidad_celulas = 1 #indica cuantos casilleros puede moverse una célula

beta = 30 # factor dentro de exp(-beta*deltaE) para determinar si una célula se mueve o no, el valor debería ir ajustándose según el experimento

factor_duplicación = 0.05 #indica que porcentaje de células se duplica por cada paso de montecarlo

numero_celulas = 0 #creo la variable donde guardaré el número de células

adhesion_celular = 1 #tiene en cuenta que son células epiteliales que forman una monocapa, con poca tendencia a soltarse unas de otras

#una adhesión de 1 corresponde con células que siempre se van a mover manteniendo al menos 1 primer vecino, si es 0.8, un 20% de los movimientos que lleven a perder todos los primeros vecinos son aceptados

#luego podría escribirse una función más compleja para tener en cuanta la adhesión entre células, por ejemplo tendiendo a maximizar el número de vecinos

factor_intentos_por_paso = 0.1 #es un índica de cuantos intentos de movimientos/duplicación habrá por paso de montecarlo

#un factor de 1 indica que habrá un intento por célula, un factor de 0.5 significa que habrán la mitad de intentos que de células

cada_cuanto_guarda = 50 #establece cada cuantos pasas va a guardar el frente y la imagen

celulas_restar = 0 #sirve para ir restando las celulas que vamos dejando atrás

empezar_por_fila = 0 #sirve para ir corriendo la fila donde se empiezan a dar los eventos, a medida que avanza el frente no tiene sentido seguir estudiando las células de atrás

### CREACIÓN DE LA MATRIZ INICIAL (FRENTE)

frente = [[0] * numero_columnas for i in range(numero_filas)] #crea una matriz llena de ceros con ese número de columnas y ese número de filas

#frente es la matriz o lista que tendrá siempre la información del frente de cultivo "real"

frente2 = [[0] * numero_columnas for i in range(numero_filas)]

#frente2 es la matriz o lista que tendrá siempre la información de la configuración "posible"

for j in range(numero_columnas): #barre en las columnas, es decir, a lo largo del frente
    for i in range(int(j*pendiente)+filas_iniciales_del_frente): #barre en las filas iniciales ocupadas más la corrección por inclinación
        numero_celulas = numero_celulas + densidad_inicial
        frente[i][j] = densidad_inicial ### asigna a cada elemento i,j el valor de la densidad inicial
        frente2[i][j] = densidad_inicial ### asigna a cada elemento i,j el valor de la densidad inicial



### DETERMINACIÓN de h(i) sin OH

"""
def ocho_vecinos (numero_columnas, numero_filas, frente):
    for i in (-1,0,1):
        for j in (-1,0,1):
            if not frente[i][j]>0:
                return False
            else:
                return True

#Podría servir para saber si se está analizando una posición del frente o una isla que se desprendió, pero aún tengo que pensarlo
"""

"""
def calculo_h(numero_columnas ,numero_filas , frente):
    hco=[0]*numero_columnas    #creo una matriz donde guardará las alturas con overhang
    hso=[0]*numero_columnas    #creo una matriz donde guardará las alturas sin overhang
    for j in range(numero_columnas): #van a barrerse una por una las columnas de la matriz frente
        count = 0
        for i in range(numero_filas): #van a barrerse todas las filas para ir viendo cuales están ocupadas y cuales no
            if frente[i][j]>0: #si el casillero tiene un nÃºmero mayor a 0, estÃ¡ ocupado, más allá de si está en contacto o no con el frente !!!!!!!!!!
                count = count + 1
                hco[j]=i+1 #sin importar el número de casilleros vacíos, el valor de la altura es el de la mayor fila ocupada
        hso[j]=count #es el número de casilleros ocupados en la columna, sin importar si esta o no conectados al frente !!!!
    return (hco, hso) #la función devuelve una lista con los valores de h sin OH y con OH para cada columna

Se guardarán las alturas sin OH para cada i. Como es sin OH, corresponde al número de filas (casilleros) ocupados y que están

en contacto con el frente a través de vecinos (es decir, que no son ocupados por células que se escaparon del frente)

/// COMO HAGO ESTO???????? ///

, independientemente de si hay o no hay huecos. La matriz tiene un número de filas igual al número de columnas del frente
"""

### DUPLICACIóN Y MIGRACIóN

""" En cada evento se sorterá una posición.

HABRía QUE PENSAR COMO HACER PARA QUE SE SORTEEN CON MAYOR PROBABILIDAD LAS POSICIONES MáS OCUPADAS (QUIZá ASIGNáNDOLE A CADA CéLULA/OBJETO UN
NúMERO DE IDENTIFICACIóN Y LUEGO SORTEAR EN BASE A ESO????)

En cada posición ocurrirá una división o se moverá una célula. La probabilidad de que ocurrá uno u otro está dado por Probabilidad_D_M

HABRá QUE TENER EN CUENTA EL NúMERO DE CéLULAS EN UN CASILLERO??? MáS CéLULAS INDICAN MAYOR PROBABILIDAD DE DUPLICACIóN O QUE LA DENSIDAD CELULAR
SEA MAYOR REDUCE EL TIEMPO DE DUPLICACIóN???????

Si sale duplicación, habrá que sortear nuevamente para ver si la duplicación ocurre o no, en función del número de células (pero mÃ¡s permisivo que la migración) Pensar


Si sale sorteado migración, se sorteará una posición al azar respecto a la celula:

1       2      3

4     celula   5

6       7      8

Luego, se evaluará si la célula se mueve.

Movimiento que disminuye energía: siempre aceptado

Movimiento que aumenta energía: aceptado con cierta probabilidad exp(-B*DE)
                                rechazar con probabilidad 1-exp(-B*DE)

La E del sistema está dada por el módulo del gradiente de densidad celular, el sistema deberÃ­a tender a minimizar el gradiente.

Para evaluar cada movimiento habrÃ­a que calcular un gradiente local (con los vecinos de la cÃ©lula) o general (de todo el sistema o matriz)

B habla de que tan permisivo soy con los movimientos "raros" o en contra del gradiente. A Mayor B mÃ¡s permisivo soy con esos movimientos
"""

#MANEJO DE ARRAYS

def calcular_módulo(lista): #defino la función que calcula el gradiente del módulo del frente y lo devuelve
    array = numpy.asarray(lista) #convierte una lista en un array
    gradiente = numpy.gradient(array) #calcula el gradiente del array
    return numpy.linalg.norm(gradiente) #la función devuelve el módulo del gradiente

def print_array(lista): #función que toma una lista, la transforma en array y la imprime
    print(numpy.asarray(lista))

#CONTEO DE NUMERO DE VECINOS

def numero_vecinos(numero_filas, numero_columnas, nueva_fila, nueva_columna, frente): #FUNCIÓN QUE CALCULA EL NÚMERO DE PRIMEROS Y SEGUNDOS VECINOS

                    numero_primeros_vecinos = 0 #crea la variable local primeros vecinos

                    numero_segundos_vecinos = 0 #crea la variable local segundos vecinos

                    for ivecino in (-1,0,1): #recorro las posibles posiciones vecinas
                        for jvecino in (-1,0,1):

                            nueva_fila2 = nueva_fila + ivecino #defino a nueva_fila2 que vendría ser la ubicación en fila del posible vecino

                            nueva_columna2 = nueva_columna + jvecino #defino a nueva_columna2 que vendría a ser la ubicación en columna del posible vecino

                            if nueva_fila2 >= numero_filas: #si la nueva fila tiene un índice mayor al número de filas (recordar que python empieza a contar desde 0)
                                nueva_fila2 = numero_filas-1 #defino a la nueva fila como la ultima fila (total nunca el frente debería llegar a recorrer todo el espacio, y debería estar ocupado por ceros)
                                print("Te quedaste corto con la cápsula!!")

                            elif nueva_fila2 < 0: #la posición no puede ser negativa
                                nueva_fila2 = 0 #asigno a la nueva fila la fila 0, total ahí  todo el frente debería estar lleno de células, es la región del bulk

                                #puede pensarse como que es el bulk y las células no tienen mucho movimiento, que las que retroceden se compensan por las que vienen de atrás
                                #o puede pensarse como que está la pared de la cápsula

                            if nueva_columna2 >= numero_columnas: # si la célula intenta irse por la derecha, debería aparecer por la izquierda, es mi condición de borde
                                nueva_columna2 = nueva_columna2-numero_columnas #con esto le hago dar la vuelta, si el frente tiene 10 columnas, el indice de la ultima será 9. Si nueva columna es 10, y le resto 10, la nueva posición será la 0

                            elif nueva_columna2 < 0:
                                nueva_columna2 = nueva_columna2 + numero_columnas #así, si cae en la posición -1 y le sumo 10, pasaría a estar en la 9

                            if not (ivecino == 0 and jvecino == 0): #realiza la condición si la posición a escanear no es la misma que ocupa la célula
                                if frente[nueva_fila2][nueva_columna2]>0 and (ivecino == 0 or jvecino == 0): #si alguno de los dos es cero es porque se trata de un primer vecino
                                    numero_primeros_vecinos = numero_primeros_vecinos + 1 #por eso le sumo uno al conteo de primeros vecinos
                                elif frente[nueva_fila2][nueva_columna2]>0 and (ivecino != 0 and jvecino != 0): #si las posiciones son distintas de 0, se trata de segundos vecinos
                                    numero_segundos_vecinos = numero_segundos_vecinos +1


                    if numero_primeros_vecinos > 0 and numero_segundos_vecinos > 0:
                        return True #si el número de primeros vecinos es mayor a 0, devuelve verdadero
                    else:
                        return False #sino, devuelve vverdadero


tiempo = time() #crea la variable tiempo y le asigna el valor del tiempo antes de comenzar los pasos de montecarlo


tiempostr = str(time())+".txt"

ultima_fila_max = filas_iniciales_del_frente + 1 #en esta variable voy a guardar la fila más grande ocupada de todo el frente, para no hacer sorteos al pedo

ultima_fila_max = ultima_fila_max - 1 #tengo que definirlo así para que me apunten a espacios de memoria distintos

fdata = open("data" + tiempostr,"a")
fdata.write("numero de columnas = ")
fdata.write(str(numero_columnas))
fdata.write("\n")
fdata.write("numero de filas = ")
fdata.write(str(numero_filas))
fdata.write("\n")
fdata.write("densidad maxima = ")
fdata.write(str(densidad_maxima))
fdata.write("\n")
fdata.write("densidad inicial = ")
fdata.write(str(densidad_inicial))
fdata.write("\n")
fdata.write("pendiente = ")
fdata.write(str(pendiente))
fdata.write("\n")
fdata.write("filas iniciales del frente = ")
fdata.write(str(filas_iniciales_del_frente))
fdata.write("\n")
fdata.write("movilidad celulas = ")
fdata.write(str(movilidad_celulas))
fdata.write("\n")
fdata.write("beta = ")
fdata.write(str(beta))
fdata.write("\n")
fdata.write("factor duplicacion = ")
fdata.write(str(factor_duplicación))
fdata.write("\n")
fdata.write("adhesion celular = ")
fdata.write(str(adhesion_celular))
fdata.write("\n")
fdata.write("factor intento por paso = ")
fdata.write(str(factor_intentos_por_paso))
fdata.write("\n")
fdata.close()

salto_paso = 1

for p in range(numero_pasos): # realiza el número de pasos de montecarlo deseados

    print(p)

    if 100 > p >= 10:
        salto_paso = 2
    elif 1000 > p >= 100:
        salto_paso = 20
    elif p >= 1000:
        salto_paso = 200


    if p % salto_paso == 0:

        suma_filas = 0 #va a guardar en esta variable la suma de todas las ultimas filas ocupadas para luego calcular la altura promedio

        hco=[0]*numero_columnas    #creo una matriz donde guardaré las alturas con overhang

        hso=[0]*numero_columnas    #creo una matriz donde guardaré las alturas sin overhang

        hco_promedio = 0 #crea la variable donde guardaré la altura con overhang promedio

        hso_promedio = 0 #crea la variable donde guardaré la altura sin overhang promedio

        sumatoria_hco = 0 #creo las variables donde irá guardando la suma de las distintas alturas del frente con

        sumatoria_hso = 0 # y sin overhang


    #OBTENCIÓN DEL CONTORNO DE LA COLONIA A PARTIR DE LA MATRIZ FRENTE

    #Leer tesis de belén moglia. Las células que están adheridas mediante primeros vecinos a la que está en la posición (0,0), son parte de la colonia o frente

    #Las células que no tienen conexión con la (0,0) son células/islas que se despegaron de la colonia

    #No nos importa si hay agujeros dentro de la colonia, solo nos interesa estudiar el frente

    #Para determinar el frente: Se recorre la primera columna desde la (0,0) hasta llegar al primer espacio vacío, luego se irá moviendo la posición buscando los primeros vecinos

    #Para asegurarnos de recorrer todo el frente se establece un órden de prioridades:

    #si el último movimiento fue hacia arriba (1), las prioridades serán: izquierda, arriba, derecha, abajo

    #si el último movimiento fue hacia la izquierda (2), las prioridades serán: abajo, izquierda, arriba, derecha

    #si el último movimiento fue hacie la derecha(3), serán arriba, derecha, abajo, izquierda

    #y si fue hacia abajo (4) serán derecha, abajo, izquierda, arriba


        recorrido = 0 #variable que me indica si se terminó de recorrer el frente o no, para obtener el contorno

        fila_inicial = 0 #indica desde que fila se comienza el recorrido

        columna_inicial = 0 #indica desde que columna se comienza el recorrido

        ultimo_movimiento = 1 #guarda el último movimiento, en este caso siempre se comienza con arriba

        frente3 = [[0] * numero_columnas for i in range(numero_filas)] #realizo una matriz igual a la del frente pero vacía, para ir guardando como 1 los espacios que se corresponden con bordes

        while recorrido == 0: #se hará un loop hasta que el valor de recorrido sea 1, es decir, se haya recorrido todo el borde

            if columna_inicial == numero_columnas -1: #si la columna inicial es la última columna del frente

                if frente[fila_inicial+1][columna_inicial] > 0: #sólo buscaremos intentar subir más filas

                    ##### NOTAR QUE DE ESTA MANERA NO SERÁ POSIBLE DETECTAR SALIENTES U OVERHANGS EN UNO DE LOS BORDES DEL FRENTE
                    ### POR AHORA NECESITO QUE SEA ASÍ PARA PODER PONER UN PUNTO FINAL AL BARRIDO

                    frente3[fila_inicial][columna_inicial] = 1 #Si encuentra un casillero ocupado, a esa posición en frente 3 le pone un 1
                    fila_inicial = fila_inicial + 1 #sube en 1 la fila inicial para evaluar luego más arriba

                else: #si la fila por encima de la que ya habíamos evaluado está vacía
                    frente3[fila_inicial][columna_inicial] = 1 #a la anterior le asignamos el valor de 1
                    recorrido = 1 #terminamos el recorrido, porque ya llegamos a la última fila ocupada de la última columna


            else: #si no estamos en la última columna
                if ultimo_movimiento == 1: #y el último movimiento fue hacia arriba
                    if columna_inicial >0 and frente[fila_inicial][columna_inicial-1]>0: #primero intentamos movernos hacia la izquierda, siempre y cuando exista una columna a la izquierda
                        frente3[fila_inicial][columna_inicial]=1 #en caso de haber, a la posición anterior le ponemos 1
                        columna_inicial = columna_inicial -1 #le restamos 1 a la variable de columna que uso para barrer
                        ultimo_movimiento = 2 #aclaro que el último movimiento que hice fue hacia la izquierda
                    elif frente[fila_inicial+1][columna_inicial]>0: #si no hay nada a la izquierda, intentará fijarse arriba
                        if columna_inicial > 0:#si estamos por encima de la columna inicial (para no llenar la columna inicial de 1)
                            frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial +1
                        ultimo_movimiento = 1
                    elif frente[fila_inicial][columna_inicial+1]>0: #si no hay nada, veremos a la derecha
                        frente3[fila_inicial][columna_inicial]=1
                        columna_inicial = columna_inicial +1
                        ultimo_movimiento = 3
                    elif frente[fila_inicial-1][columna_inicial]: #como última posibilidad, nos fijamos abajo, donde seguro hay algo porque el último movimiento había sido hacia arriba
                        frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial -1
                        ultimo_movimiento = 4

                elif ultimo_movimiento == 2: #si el último movimiento no fue hacia arriba, testea si fue hacia la izquierda.
                    if frente[fila_inicial-1][columna_inicial]: #luego realizo el mismo procedimiento anterior, pero asignando las prioridades que correspondan
                        frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial -1
                        ultimo_movimiento = 4
                    elif columna_inicial >0 and frente[fila_inicial][columna_inicial-1]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        columna_inicial = columna_inicial -1
                        ultimo_movimiento = 2
                    elif frente[fila_inicial+1][columna_inicial]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial +1
                        ultimo_movimiento = 1
                    elif frente[fila_inicial][columna_inicial+1]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        columna_inicial = columna_inicial +1
                        ultimo_movimiento = 3

                elif ultimo_movimiento == 3:
                    if frente[fila_inicial+1][columna_inicial]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial +1
                        ultimo_movimiento = 1
                    elif frente[fila_inicial][columna_inicial+1]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        columna_inicial = columna_inicial +1
                        ultimo_movimiento = 3
                    elif frente[fila_inicial-1][columna_inicial]:
                        frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial -1
                        ultimo_movimiento = 4
                    elif columna_inicial >0 and frente[fila_inicial][columna_inicial-1]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        columna_inicial = columna_inicial -1
                        ultimo_movimiento = 2

                elif ultimo_movimiento == 4:
                    if frente[fila_inicial][columna_inicial+1]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        columna_inicial = columna_inicial +1
                        ultimo_movimiento = 3
                    elif frente[fila_inicial-1][columna_inicial]:
                        frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial -1
                        ultimo_movimiento = 4
                    elif columna_inicial >0 and frente[fila_inicial][columna_inicial-1]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        columna_inicial = columna_inicial -1
                        ultimo_movimiento = 2
                    elif frente[fila_inicial+1][columna_inicial]>0:
                        frente3[fila_inicial][columna_inicial]=1
                        fila_inicial = fila_inicial +1
                        ultimo_movimiento = 1


        frente4 = [[0] * numero_columnas for i in range(numero_filas)]

        for j in range(numero_columnas): #ahora voy a recorrer el frente 3, para recortar el frente y dejar sólo los células que componen la "tierra", hasta el borde
            unos = 0 #creo la variable unos donde guardaré el número de 1 que se encuentra recoriendo el frente
            for i in range(ultima_fila_max): #en filas recorre hasta la que sabemos es la última ocupada
                if frente3[i][j] == 1: #y voy sumando el número de 1 que aparecen como parte del frente
                    unos = unos +1

            hay_uno = 0 #inicializo la variable hay_uno, que vamos a usarla como contador
            for i in range(ultima_fila_max): #recorro nuevamente las filas
                if hay_uno != unos: #si hay_uno es distinto del número de 1 que sabemos que hay voy a copiar los valores en el frente4
                #si hay_uno alcanza el número de unos, quiere decir que llegamos al borde del frente y el resto deberían ser 0
                    if hay_uno%2 == 0: #y si el resto de la división de hay_uno por 2 es igual a 0, es decir, el primer 1 va a ser el borde del frente
                        #si luego aparece otro 1 es un overhang. En ese caso entre los dos 1, tienen que haber ceros (0) y el resto entre ellos va a ser 1
                        frente4[i][j]= frente[i][j] * 1 #si el resto es 0, a frente4 que va a ser el frente recortado, le doy el valor que tiene el frente en ese punto
                    if frente3[i][j]==1: #si el valor de frente3 es igual a 1, se trata de uno borde de colonia u overhang
                        hay_uno = hay_uno + 1 #en ese caso le sumamos 1 al contador hay_uno
                        frente4[i][j] = frente[i][j] * 1 #y le asignamos a esa posición de frente4 el valor de frente

        for j in range(numero_columnas): #voy a recorrer el numero de columnas y filas para saber cuantas celulas hay y cual es la mayor fila ocupada

            filas_ocupadas = 0 #crea la variable donde voy a ir guardando la cantidad de filas ocupadas, osea, SIN OVERHANG

            ultima_fila_ocupada = 0 #crea la variable donde voy a guardar la ultima fila ocupada, independientemente de si para llegar a estas hay espacios en blanco o no, osea, CON OVERHANG Â¡Â¡Â¡CHEQUEAR!!!

            for i in range(ultima_fila_max): #recorre el numero de filas para contar las celulas del frente y ver cual es la mayor fila ocupada

                if frente4[i][j]>0: #si el casillero tiene un número mayor a 0, está ocupado, más allá de si está en contacto o no con el frente !!!!!!!!!!
                    filas_ocupadas = filas_ocupadas + 1 #va agregando uno por cada fila ocupada
                    ultima_fila_ocupada = i + 1

            hco[j] = ultima_fila_ocupada #asigna a cada numero de columna la altura en filas con overhang

            hso[j] = filas_ocupadas #asigna a cada número de columna el número de filas ocupadas o altura sin overhang

            sumatoria_hco = sumatoria_hco + hco[j] #realizao la sumatoria de las alturas de cada columna para luego dividir por el número de columnas y obtener la altura promedio

            sumatoria_hso = sumatoria_hso + hso[j]


        hco_promedio = sumatoria_hco / float(numero_columnas) #altura promedio con overhang

        hso_promedio = sumatoria_hso / float(numero_columnas) #altura promedio sin overhang

        """
        rugosidad w

        w (t) = sqrt[1/N*sumatoria(h(t)-hpromedio)cuadrado]

        """

        suma_cuadrados_hco = 0 #creamos las variables para guardar las sumatoria de las desviaciones de alturas co y so, y las seteamos a 0 todos los pasos

        suma_cuadrados_hso = 0

        rugosidad_co = 0 #creamos las variables para guardar los valores de rugosidad y los seteamos a 0 todos los pasos

        rugosidad_so = 0


        for j in range(numero_columnas): #barriendo en las columnas
            suma_cuadrados_hco = suma_cuadrados_hco + (hco[j] - hco_promedio)**2 #sumamos los cuadrados de las desviaciones de las alturas
            suma_cuadrados_hso = suma_cuadrados_hso + (hso[j] - hso_promedio)**2

        rugosidad_co = sqrt(suma_cuadrados_hco/numero_columnas) #calculamos la rugosidad utilizando la función sqrt

        rugosidad_so = sqrt(suma_cuadrados_hso/numero_columnas)


        f = open(tiempostr,"a")
        f.write(str(time()))
        f.write(" ")
        f.write(str(p))
        f.write(" ")
        f.write(str(numero_celulas))
        f.write(" ")
        f.write(str(hco_promedio))
        f.write(" ")
        f.write(str(hso_promedio))
        f.write(" ")
        f.write(str(rugosidad_co))
        f.write(" ")
        f.write(str(rugosidad_so))
        f.write("\n")
        f.close()



    if p % cada_cuanto_guarda == 0:

        """
        fr = open(str(tiempo)+" "+str(p)+".txt","a")
        fr.write("\n")
        for i in range(numero_filas):
            for j in range(numero_columnas):
                fr.write(str(frente4[i][j])+" ")
            fr.write("\n")
        fr.close()
        """

        plt.imshow(numpy.asarray(frente), cmap='viridis', aspect = "equal", interpolation = "none", vmin = 0, vmax = densidad_maxima)

        plt.savefig(str(tiempo)+" "+str(p) + "frente.pdf", bbox_inches='tight', transparent = "True", frameon = "False")


        plt.imshow(numpy.asarray(frente4), cmap='viridis', aspect = "equal", interpolation = "none", vmin = 0, vmax = densidad_maxima)

        plt.savefig(str(tiempo)+" "+str(p)+ "frente4.pdf", bbox_inches='tight', transparent = "True", frameon = "False")


    conteo_numero_celulas = numero_columnas * filas_iniciales_del_frente * densidad_inicial * factor_intentos_por_paso #int((numero_celulas - celulas_restar) * factor_intentos_por_paso)


    while conteo_numero_celulas > 0: #van a haber tantos intentos como células en el frente

        sorteo = random.uniform(0,1) #sortea de antemano si va a ocurrir un posible evento de duplicaciòn o uno de movilidad

        E1 = calcular_módulo(frente) #calcula el módulo del gradiente del frente para utilizar luego

        if sorteo < factor_duplicación : #si el sorteo cae en cierto rango de valores, quizà ocurre una duplicaciòn (habrà que ajustar el intervalo)
            i = random.randrange(empezar_por_fila, ultima_fila_max) #sortea una posible fila
            j = random.randrange(numero_columnas) #sortea una posible columna

            if frente[i][j] > 0: #and numero_vecinos(numero_filas, numero_columnas, i, j, frente): #se fija si hay una célula en el lugar para ver si realiza una duplicación, y tambión si ya ha alcanzado la densidad celular máxima

                conteo_numero_celulas = conteo_numero_celulas - 1 #le resto 1 al número de células para terminar en algún momento el loop while


                if densidad_maxima >= frente[i][j] and numero_vecinos(numero_filas, numero_columnas, i, j, frente):


                    frente[i][j] = frente[i][j]+1 #a la posición sorteada de los frentes le sumo 1

                    frente2[i][j] = frente2[i][j]+1

                    numero_celulas = numero_celulas + 1 #aumenta en 1 el número de células

                """
                if random.uniform(0,1) < (frente[i][j]/densidad_maxima) and densidad_maxima > frente[i][j]: #con esto intentamos dos cosas:

                    #primero que la densidad de un casillero nunca supere la densidad máxima

                    #segundo, que ocurran más eventos de duplicación en los casilleros más poblados y menos en los menos poblados (que pueden ser células desprendidas del frente que me generan islas)

                    frente[i][j] = frente[i][j]+1 #a la posición sorteada de los frentes le sumo 1

                    frente2[i][j] = frente2[i][j]+1

                    numero_celulas = numero_celulas + 1 #aumenta en 1 el número de células


                Esto estaba pensado para intentar disminuir la formación de islas, pero cambia el exponente del escalado de 0.33 a 0.25
                """

        else: #

            i = random.randrange(empezar_por_fila, ultima_fila_max) #sortea una posible fila, entre la primera donde quiero empezar a ver eventos y la última ocupada
            j = random.randrange(numero_columnas) #sortea una posible columna
            diferencia_ultima_fila_max_nueva_fila = 0 #variable donde voy a guardar la diferencia entre la ultima fila max y la nueva fila, para agregar el numero de filas que corresponda

            if frente[i][j] > 0: #si en ese lugar hay una célula

                conteo_numero_celulas = conteo_numero_celulas - 1 #le resta 1 al contador, independientemente de que se produzca o no movimiento

                desplazamiento_fila = (random.randrange(1+2*movilidad_celulas))-movilidad_celulas #sortea para ver a donde se movera la celula, la resta se incluye para que puedan existir numeros negativos
                desplazamiento_columna = (random.randrange(1+2*movilidad_celulas))-movilidad_celulas

                frente2[i][j] = frente2[i][j]-1 #al némero de células en el casillero sorteado le resta la célula que se mueve

                nueva_fila = i+desplazamiento_fila #defino cual será la nueva fila que ocuparía la célula si se mueve
                nueva_columna = j+desplazamiento_columna #defino cual será la nueva columna que ocupará la célula si se mueve

                if nueva_fila >= numero_filas: #si la nueva fila tiene un índice mayor al número de filas (recordar que python empieza a contar desde 0)
                    nueva_fila = numero_filas-1 #defino a la nueva fila como la ultima fila (total nunca el frente debería llegar a recorrer todo el espacio, y debería estar ocupado por ceros)
                    print("Te quedaste corto con la cápsula!!!")

                elif nueva_fila < 0: #la posición no puede ser negativa
                    nueva_fila = 0 #asigno a la nueva fila la fila 0, total ahí  todo el frente debería estar lleno de células, es la región del bulk

        #


                if nueva_columna >= numero_columnas: # si la célula intenta irse por la derecha, debería aparecer por la izquierda, es mi condición de borde
                    nueva_columna = nueva_columna-numero_columnas #con esto le hago dar la vuelta, si el frente tiene 10 columnas, el indice de la ultima será 9. Si nueva columna es 10, y le resto 10, la nueva posición será la 0

                elif nueva_columna < 0:
                    nueva_columna = nueva_columna + numero_columnas #así, si cae en la posicion -1 y le sumo 10, pasaría a estar en la 9



                frente2[nueva_fila][nueva_columna] = frente2[nueva_fila][nueva_columna]+1 #al casillero donde se movería la célula se le suma 1
                E2 = calcular_módulo(frente2) #calcula el modulo del gradiente del frente que se obtendría con la nueva configuración

                DeltaE = E2 - E1 #calculo y guardo el DeltaE

                if numero_vecinos(numero_filas, numero_columnas, nueva_fila, nueva_columna, frente):  #si la célula mantiene al menos 1 primer vecino

                    if DeltaE <= 0: #si el DeltaE es igual o menor que 0, acepto la nueva configuración
                        frente[nueva_fila][nueva_columna] = frente[nueva_fila][nueva_columna]+1 #al frente se le suma una célula en la nueva posición
                        frente[i][j] = frente[i][j]-1 #se le resta la célula que se mueve
                        if (nueva_fila+1) > ultima_fila_max: #si la nueva fila de la célula está por encima de todas las filas anteriores
                                diferencia_ultima_fila_max_nueva_fila = ultima_fila_max - nueva_fila +1 #le doy un valor a la diferencia
                                ultima_fila_max = nueva_fila + 1 #esa nueva fila pasa a ser la fila máxima
                                for i in range(diferencia_ultima_fila_max_nueva_fila):
                                    frente.append([0] * numero_columnas) #si la diferencia es uno, agrego una fila, si es dos, agrego 2
                                    frente2.append([0] * numero_columnas)
                                    numero_filas = numero_filas+1
                                    empezar_por_fila = empezar_por_fila + 1
                                    for j in range(numero_columnas): #con esto voy acumulando las células que le tengo que ir restando al contador de eventos porque las fui dejando atrás
                                        celulas_restar = celulas_restar + frente[empezar_por_fila][j]


                    elif random.uniform(0,1) < (exp(-1*beta*DeltaE)): #si el DeltaE es positivo, el movimiento no disminuye el gradiente y va a ser aceptado o rechazado con cierta probabilidad
                        #la probabilidad de aceptación depende de cuan positivo sea DeltaE y cuál sea al valor de beta
                        frente[nueva_fila][nueva_columna] = frente[nueva_fila][nueva_columna]+1 #si el movimiento es aceptado, se actualiza el frente
                        frente[i][j] = frente[i][j]-1
                        if (nueva_fila+1) > ultima_fila_max: #si la nueva fila de la célula está por encima de todas las filas anteriores
                                diferencia_ultima_fila_max_nueva_fila = ultima_fila_max - nueva_fila +1 #le doy un valor a la diferencia
                                ultima_fila_max = nueva_fila + 1 #esa nueva fila pasa a ser la fila máxima
                                for i in range(diferencia_ultima_fila_max_nueva_fila):
                                    frente.append([0] * numero_columnas) #si la diferencia es uno, agrego una fila, si es dos, agrego 2
                                    frente2.append([0] * numero_columnas)
                                    numero_filas = numero_filas+1
                                    empezar_por_fila = empezar_por_fila + 1
                                    for j in range(numero_columnas):
                                        celulas_restar = celulas_restar + frente[empezar_por_fila][j]

                    else: #si el movimiento no es aceptado, el frente2 se corrije para volver a ser igual al frente
                        frente2[nueva_fila][nueva_columna] = frente2[nueva_fila][nueva_columna]-1
                        frente2[i][j] = frente2[i][j]+1

                elif random.uniform(0,1) > adhesion_celular: #si la célula no mantiene primeros vecinos es como si se estuviera despegando del frente
                    #en ese caso tiene que haber un sorteo, y si el sorteo supera el valor de adhesión celular, se le permite a la célula moverse

                    if DeltaE <= 0: #IDEM arriba
                        frente[nueva_fila][nueva_columna] = frente[nueva_fila][nueva_columna]+1
                        frente[i][j] = frente[i][j]-1
                        if (nueva_fila+1) > ultima_fila_max: #si la nueva fila de la célula está por encima de todas las filas anteriores
                                diferencia_ultima_fila_max_nueva_fila = ultima_fila_max - nueva_fila +1 #le doy un valor a la diferencia
                                ultima_fila_max = nueva_fila + 1 #esa nueva fila pasa a ser la fila máxima
                                for i in range(diferencia_ultima_fila_max_nueva_fila):
                                    frente.append([0] * numero_columnas) #si la diferencia es uno, agrego una fila, si es dos, agrego 2
                                    frente2.append([0] * numero_columnas)
                                    numero_filas = numero_filas+1
                                    empezar_por_fila = empezar_por_fila + 1
                                    for j in range(numero_columnas):
                                        celulas_restar = celulas_restar + frente[empezar_por_fila][j]

                    elif random.uniform(0,1) < (exp(-1*beta*DeltaE)):
                        frente[nueva_fila][nueva_columna] = frente[nueva_fila][nueva_columna]+1
                        frente[i][j] = frente[i][j]-1
                        if (nueva_fila+1) > ultima_fila_max: #si la nueva fila de la célula está por encima de todas las filas anteriores
                                diferencia_ultima_fila_max_nueva_fila = ultima_fila_max - nueva_fila +1 #le doy un valor a la diferencia
                                ultima_fila_max = nueva_fila + 1 #esa nueva fila pasa a ser la fila máxima
                                for i in range(diferencia_ultima_fila_max_nueva_fila):
                                    frente.append([0] * numero_columnas) #si la diferencia es uno, agrego una fila, si es dos, agrego 2
                                    frente2.append([0] * numero_columnas)
                                    numero_filas = numero_filas+1
                                    empezar_por_fila = empezar_por_fila + 1
                                    for j in range(numero_columnas):
                                        celulas_restar = celulas_restar + frente[empezar_por_fila][j]

                    else:
                        frente2[nueva_fila][nueva_columna] = frente2[nueva_fila][nueva_columna]-1
                        frente2[i][j] = frente2[i][j]+1

                else: #si por no tener vecinos la célula no se mueve, se devuelve el frente2 a la configuración anterior
                    frente2[nueva_fila][nueva_columna] = frente2[nueva_fila][nueva_columna]-1
                    frente2[i][j] = frente2[i][j]+1


print("Terminado")
