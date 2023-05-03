#llamado de librerias y del programa de lectura de sintaxis (parser)
import sympy
from sympy import Eq,sympify,lambdify,And
import re
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations

#Definicion variables, matrices y listas
x, y = sympy.symbols('x y')
Valor_ecuaciones = ["y=0","x=0"]
Lista_inecuaciones = [x>=0, y>=0]
x_vals = []
y_vals = []
valores_maximos_x_y = [0,0]
lista_restricciones = []

def transformar_funcion(funcion):
    funcion = re.sub(r'(\d)([xy])', r'\g<1>*\2', funcion)
    # Reemplazar ' ' sobrantes
    funcion = funcion.replace(' ', '')
    return funcion

def Puntos(Valor_ecuaciones, Lista_inecuaciones,Funcion,Bandera_Max_Min):
    #iniciamos los valores minimos de cada variable

    if Bandera_Max_Min == 1:
        sol_optima = (-np.inf, -np.inf)
        valor_optimo = -np.inf
    else:
        sol_optima = (np.inf,np.inf)
        valor_optimo = np.inf

    for i in range(len(Valor_ecuaciones)):
        # Separamos la ecuación en cada signo de desigualdad
        if '<=' in Valor_ecuaciones[i]:
            parts = Valor_ecuaciones[i].split('<=')
        elif '>=' in Valor_ecuaciones[i]:
            parts = Valor_ecuaciones[i].split('>=')
        else:
            parts = Valor_ecuaciones[i].split('=')

        # Limpiamos los elementos de la lista
        parts = [p.strip() for p in parts]

        # Simplificamos el valor
        parts[1] = sympy.simplify(parts[1])
        parts[0] = sympy.simplify(parts[0])

        # Modificamos la lista original, se transforma a una ecuacion lineal
        lista_restricciones.append(parts)
        Valor_ecuaciones[i] = parts
        print(type(parts))
        Valor_ecuaciones[i] = Eq(Valor_ecuaciones[i][0], Valor_ecuaciones[i][1])
    #combinamos todas las ecuaciones resultantes para calcular los puntos de cruce
    combinaciones = list(combinations(Valor_ecuaciones,2))
    print(lista_restricciones)
    #Iniciamos for para recorrer las restricciones y buscamos los puntos factibles
    print(combinaciones)
    for comb in combinaciones:
        sol = sympy.solve(comb, (x, y))
        if not sol:
            continue
        x_sol, y_sol = float(sol[x]), float(sol[y])
        print(sol)
        #verificamos que el valor x e y son los mas altos
        if not valores_maximos_x_y[0] > x_sol:
            valores_maximos_x_y[0] = int(x_sol)
        if not valores_maximos_x_y[1] > y_sol:
            valores_maximos_x_y[1] = int(y_sol)

        # Evaluamos las desigualdades para este punto
        if not all(desig.subs({x:x_sol, y:y_sol}) for desig in Lista_inecuaciones):
            continue
        print(x_sol, y_sol)
        print(type(x_sol))
        print(float(x_sol))

        #Transformamos la funcion objetivo a sympy y evaluamos el punto en ella
        Funcion_sympy = sympify(Funcion)
        valor = Funcion_sympy.subs({x: x_sol, y: y_sol})
        if Bandera_Max_Min == 1:
            #Verificacion de cumplimiento para valor maximo
            if not valor_optimo or valor > valor_optimo:
                sol_optima = (x_sol, y_sol)
                valor_optimo = valor
        else:
            #Verificacion de cumplimiento para valor minimo
            if valor_optimo or valor < valor_optimo:
                sol_optima = (x_sol, y_sol)
                valor_optimo = valor
        #Ingreso de todos los puntos x e y que se obtienen
        x_vals.append(x_sol)
        y_vals.append(y_sol)

    #impresion de valores para verificacion u informacion extra
    print(valores_maximos_x_y)
    print("Solución óptima:", sol_optima)
    print("Valor óptimo:", valor_optimo)
    return valor_optimo,sol_optima

def graficar(Funcion, Lista_combinada_inecuaciones, x_vals, y_vals,valor_optimo,sol_optima):

    # Convertimos la expresión booleana en una función de numpy
    f_np = lambdify((x, y), Lista_combinada_inecuaciones, 'numpy')

    # Creamos la figura y los ejes
    fig,ax = plt.subplots(figsize=(10,10))

    # Graficamos los puntos
    for i in range(len(x_vals)):
        ax.scatter(x_vals[i], y_vals[i], color = "red")
        ax.text(x_vals[i]+0.1, y_vals[i]+0.1, f'({x_vals[i]},{y_vals[i]})')

    # Graficamos la región factible
    X, Y = np.meshgrid(np.linspace(0, int(valores_maximos_x_y[0]), 1000), np.linspace(0, int(valores_maximos_x_y[1]), 1000))
    Z = f_np(X, Y)
    contourf = ax.contourf(X, Y, Z, colors=['white', 'blue'], alpha=0.5)
 
    # Graficamos la línea de nivel de la función objetivo
    X, Y = np.meshgrid(np.linspace(-10, 100, 100), np.linspace(-10, 100, 100))
    Z = eval(Funcion.upper())

    ax.contour(X, Y, Z, levels=[valor_optimo], colors='blue')
    #ax.contour(X, Y, Z2, colors='black')
    
    colors1 = plt.rcParams['axes.prop_cycle'].by_key()['color'][:30]
    print(colors1)
    for i in range(2,len(lista_restricciones)):
        ZX = eval(str(lista_restricciones[i][0]).upper()+"-"+str(lista_restricciones[i][1]))
        ax.contour(X, Y, ZX, colors=colors1[i])

    # Configuramos los ejes y mostramos la figura
  
    ax.set_xlim([-2, valores_maximos_x_y[0]])
    ax.set_ylim([-2, valores_maximos_x_y[1]])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    # remplzamos los valores x e y a los optimo.
    texto_funcion = Funcion.replace('x', str(sol_optima[0])).replace('y', str(sol_optima[1]))

    #Ocupamos proxys para mostrar losi labels y su leyenda
    proxy = [plt.Rectangle((0,0),1,1,fc=contourf.collections[1].get_facecolor()[0], alpha=0.5),
            plt.Line2D([0],[0], linestyle="none", c='r', marker = 'o', label='Puntos'),
            plt.Line2D([0],[0], linestyle="solid", c='r', label='Función objetivo')]
    ax.legend(proxy, ["Región factible\n"+str(Lista_inecuaciones), "Puntos", "Función objetivo\n"+texto_funcion+" = "+str(valor_optimo)])

    plt.show()

def funcion():
    #Funcion general, se encarga de pedir datos y correr el resto de codigo
    Bandera_Max_Min = int(input("Desear maximizar(1) o minimizar(2): "))
    Funcion = input("Ingrese la funcion en orden(N*x+N*y): ")
    Funcion = transformar_funcion(Funcion)
    print(Funcion)
    while True:
        restricciones = input("Ingrese el número de restricciones: ")
        if restricciones.isdigit():
            Restricciones = int(restricciones)
            break
        else:
            print("Por favor, ingrese un número entero válido.")

    for x in range(Restricciones):
        input_valor = input("Ingrese la "+str(x+1)+" restriccion en orden N*x + N*y <= N: ")
        input_valor = transformar_funcion(input_valor)
        #Transofrmador de funcion 3x pasa a 3*x

        #ingreso de restricciones en formato str y sympy a las lsitas correspondientes
        Lista_inecuaciones.append(sympify(input_valor))
        Valor_ecuaciones.append(input_valor)

    #union de restricciones
    Lista_combinada_inecuaciones = And(*Lista_inecuaciones)

    #obtencion de los puntos optimos para min y max, y graficacion
    valor_optimo,sol_optima = Puntos(Valor_ecuaciones, Lista_inecuaciones,Funcion,Bandera_Max_Min)
    graficar(Funcion,Lista_combinada_inecuaciones,x_vals, y_vals,valor_optimo,sol_optima) 

#llamado de funciones

funcion()





