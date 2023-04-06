import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString


#Funcion: z=2x+2y = // -2y=2x // 2y=-2x // y= (-2x)/2
#Restriccion: 2*x+y <= 100
#x+3y <= 80
#x <= 45
#y <= 100
#x,y >=0

#Paso 1: transformacion Restricciones:
#y1=(100-2*x)
#y2=(80-x)/3
#x=45+(0*y) 
#y3=100+(0*x)
#y4= 0*x
#x2= 0*y

#paso 1.5 = transformar funcion:
#y5= (-2x)/2

#Funcion z = 10000x + 6000y  (Maximizar)
#Restricciones:
#20x + 50y <= 3000 (Equipaje permitido)
#x + y <= 90 (Cantidad de asientos disponibles)
#y >= 10 (Política de asientos mínimos para no fumadores)
#y >= 0 (No negatividad)
#x >= 0 (No negatividad)

#Paso 1: transformarla en ecuaciones normales respecto a Y las restricciones
#y1 = (3000 – 20x) / 50



FuncionOb= ()
def grafico():
        x = np.arange(-100, 150, 50)
        y = np.arange(-100, 150, 50)
        print(x)
        y5= (-2*x)/2 #Funcion objetivo
        y1=100-(2*x)
        y2=(80-x)/3
        x1=45+(0*y) 
        y3=100+(0*x)
        y4= 0*x
        x2= 0*y
        #Identificadores para las líneas
        primera_linea = LineString(np.column_stack((x, y1)))
        segunda_linea = LineString(np.column_stack((x, y2)))
        tercera_linea = LineString(np.column_stack((x, y3)))
        cuarta_linea = LineString(np.column_stack((x1, y))) 
        quinta_linea = LineString(np.column_stack((x, y4)))
        sexta_linea = LineString(np.column_stack((x2, y)))
        Sectima_linea = LineString(np.column_stack((x, y5)))
        octava_linea = LineString(np.column_stack((x, y5)))

        #Graficando las líneas

        plt.plot(x, y1, '-', linewidth=2, color='b')
        plt.plot(x, y2, '-', linewidth=2, color='g')
        plt.plot(x, y3, '-', linewidth=2, color='r')
        plt.plot(x1, y, '-', linewidth=2, color='y')
        plt.plot(x, y4, '-', linewidth=2, color='k')
        plt.plot(x2, y, '-', linewidth=2, color='c') #No negatividad
        #plt.plot(x, y5, ':', linewidth=1, color='r') #Funcion OBjetivo

        #plt.plot(0, 100, ':', linewidth=2, color='aquamarine')

        plt.grid()
        plt.title('Método Gráfico, la perfumeria')

        #Generando las intersecciones (vértices)
        #si importa el orden, debe ser al sentido del reloj

        primera_interseccion = sexta_linea.intersection(segunda_linea) #cian con verde
        segunda_interseccion = segunda_linea.intersection(primera_linea) #verde con azul
        #tercera_interseccion = segunda_linea.intersection(cuarta_linea) #verde con amarillo
        cuarta_interseccion = primera_linea.intersection(cuarta_linea) #azul con amarillo
        quinta_interseccion = cuarta_linea.intersection(quinta_linea) #amarillo con negro
        sexta_interseccion = quinta_linea.intersection(sexta_linea) #negro con cian


        #Graficando los vértices
        plt.plot(*primera_interseccion.xy, 'o', color='k')
        #plt.annotate(primera_interseccion,*primera_interseccion.xy)
        plt.plot(*segunda_interseccion.xy, 'o', color='k')
        #plt.plot(*tercera_interseccion.xy, 'o', color='k')
        plt.plot(*cuarta_interseccion.xy, 'o', color='k')
        plt.plot(*quinta_interseccion.xy, 'o', color='k')
        plt.plot(*sexta_interseccion.xy, 'o', color='k')
        #plt.plot(*sectima_interseccion.xy, 'o', color='k')
        #plt.plot(*octava_interseccion.xy, 'o', color='k')

        #annotations = [ primera_interseccion, segunda_interseccion, tercera_interseccion, cuarta_interseccion, quinta_interseccion, sexta_interseccion, sectima_interseccion, octava_interseccion]
        #for i, label in enumerate(annotations):
        #        plt.annotate(label, (x[i] + 0.1, y[i]))

        
        
        #Imprimiendo las coordenadas de los vértices en la consola
        print('\n COORDENADAS DE LAS INTERSECCIONES')
        print('Coordenadas de la primera intersección: {} '.format(primera_interseccion))
        print('Coordenadas de la segunda intersección: {} '.format(segunda_interseccion))
        #print('Coordenadas de la tercera intersección: {} '.format(tercera_interseccion))
        print('Coordenadas de la cuarta intersección: {} '.format(cuarta_interseccion))
        print('Coordenadas de la quinta intersección: {} '.format(quinta_interseccion))
        print('Coordenadas de la sexta intersección: {} '.format(sexta_interseccion))
        #print('Coordenadas de la sectima intersección: {} '.format(sectima_interseccion))
        #print('Coordenadas de la octava intersección: {} '.format(octava_interseccion))

        #Identificando los valores de las coordenadas (x, y) de cada vértice
        xi1m, yi1m = primera_interseccion.xy
        xi2m, yi2m = segunda_interseccion.xy
        #xi3m, yi3m = tercera_interseccion.xy
        xi4m, yi4m = cuarta_interseccion.xy

        xi5m, yi5m = quinta_interseccion.xy
        xi6m, yi6m = sexta_interseccion.xy
        #xi7m, yi7m = sectima_interseccion.xy
        #xi8m, yi8m = octava_interseccion.xy

        #Cambiamos el formato de la variable de matriz a real
        xi1 = np.float64(np.array(xi1m))
        xi2 = np.float64(np.array(xi2m))
        #xi3 = np.float64(np.array(xi3m))
        xi4 = np.float64(np.array(xi4m))
        yi1 = np.float64(np.array(yi1m))
        yi2 = np.float64(np.array(yi2m))
        #yi3 = np.float64(np.array(yi3m))
        yi4 = np.float64(np.array(yi4m))
        #Cambiamos el formato de la variable de matriz a real
        xi5 = np.float64(np.array(xi5m))
        xi6 = np.float64(np.array(xi6m))
        #xi7 = np.float64(np.array(xi7m))
        #xi8 = np.float64(np.array(xi8m))
        yi5 = np.float64(np.array(yi5m))
        yi6 = np.float64(np.array(yi6m))
        #yi7 = np.float64(np.array(yi7m))
        #yi8 = np.float64(np.array(yi8m))

        #Evaluando la función objetivo en cada vértice
        #Funcion: z=2x+2y
        FOi1 = (2 * xi1) + (2 * yi1)
        FOi2 = (2 * xi2) + (2 * yi2)
        #FOi3 = (2 * xi3) + (2 * yi3)
        FOi4 = (2 * xi4) + (2 * yi4)
        FOi5 = (2 * xi5) + (2 * yi5)
        FOi6 = (2 * xi6) + (2 * yi6)
        #FOi2 = (xi2 * 10000) + (yi2 * 6000)
        #FOi3 = (xi3 * 10000) + (yi3 * 6000)
        #FOi4 = (xi4 * 10000) + (yi4 * 6000)
        #FOi5 = (xi5 * 10000) + (yi5 * 6000)
        #FOi6 = (xi6 * 10000) + (yi6 * 6000)
        #FOi7 = (xi7 * 10000) + (yi7 * 6000)
        #FOi8 = (xi8 * 10000) + (yi8 * 6000)

        #Imprimiendo las evaluaciones de la FO en cada vértice (Consola)
        print('\n EVALUACIÓN DE LA FO EN LOS VÉRTICES')
        print('Función objetivo en la intersección 1: {} pesos'.format(FOi1))
        print('Función objetivo en la intersección 2: {} pesos'.format(FOi2))
        #print('Función objetivo en la intersección 3: {} pesos'.format(FOi3))
        print('Función objetivo en la intersección 4: {} pesos'.format(FOi4))
        print('Función objetivo en la intersección 5: {} pesos'.format(FOi5))
        print('Función objetivo en la intersección 6: {} pesos'.format(FOi6))
        #print('Función objetivo en la intersección 3: {} pesos'.format(FOi7))
        #print('Función objetivo en la intersección 4: {} pesos'.format(FOi8))

        #Calculando el mejor resultado (Maximizar)
        #ZMAX = max(FOi1, FOi2, FOi3, FOi4, FOi5, FOi6, FOi7, FOi8)
        ZMAX = max(FOi1, FOi2, FOi4, FOi5, FOi6)

        #Imprimiendo la solución óptima en la consola
        print('\n SOLUCIÓN ÓPTIMA')
        print('Solución óptima: {} pesos'.format(ZMAX))

        #Ordenando las coordenadas de los vértices (Las coordenadas x en m y las coordenadas y en n)
        #m = [xi1, xi2, xi3, xi4, xi5, xi6, xi7, xi8]
        #n = [yi1, yi2, yi3, yi4, yi5, yi6, yi7, yi8]
        m = [xi1, xi2, xi4, xi5, xi6]
        n = [yi1, yi2, yi4, yi5, yi6]


        #Graficando el polígono solución a partir de las coordenadas de los vértices 
        plt.fill(m, n, color='silver')

        #Identificando el índice del vértice de la mejor solución
        #dict1 = {0:FOi1, 1:FOi2, 2:FOi3, 3:FOi4, 4:FOi5, 5:FOi6, 6:FOi7, 7:FOi8}
        dict1 = {0:FOi1, 1:FOi2, 3:FOi4, 4:FOi5, 5:FOi6}
        posicion = max(dict1, key=dict1.get)

        #Obteniendo las coordenadas del vértice de la mejor solución de acuerdo al índice
        XMAX = m[posicion]
        YMAX = n[posicion]

        #Imprimiendo las coordenadas del vértice de la mejor solución (variables de decisión)
        print('\n VARIABLES DE DECISIÓN')
        print('Cantidad de asientos a reservar para fumadores: {} '.format(XMAX))
        print('Cantidad de asientos a reservar para no fumadores: {} '.format(YMAX))
        plt.plot(x+XMAX, y5+YMAX, ':', linewidth=1, color='m') #Funcion OBjetivo
        #Generando las anotaciones de las coordenadas y solución óptima en el gráfico
        plt.annotate('  X: {0} / Y: {1} (Coordenadas)'.format(XMAX, YMAX), (XMAX, YMAX))
        plt.annotate('  Solución óptima: {}'.format(ZMAX), (XMAX, YMAX+3))

        plt.show()

grafico()