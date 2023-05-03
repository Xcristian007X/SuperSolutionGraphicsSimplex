#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
using namespace std;

// Funcion que imprime una matriz
void printMatrix(vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

//Funcion de gauss jordan
vector<vector<double>>& GaussJordan(vector<vector<double>>& matrix, double numero_a_multiplicar, int modo, bool mixto, int acumulador){
    int m = matrix.size();    //numero de filas
    int n = matrix[0].size(); //numero de elementos en una fila o numero de columnas
    int pos_variable_decision = 0;
    if (mixto == false) {
        if (modo == 1) {             //Modo 1 aplica el metodo de gauss jordan para igualar el valor de las variables a 0
            for (int verificador = 0; verificador < n; verificador++){
                if (matrix[0][verificador] == numero_a_multiplicar) {   //Si encuentra un valor que es igual que el parametro dado (que se quiere igualar a 0)
                    //if (matrix[0][verificador] < matrix[0][verificador+1] or matrix[0][verificador] > matrix[0][verificador-1] ) {
                    //cout << matrix[0][verificador] << " " << matrix[0][verificador+1] << " " << matrix[0][verificador-1] << endl; 
                    pos_variable_decision = verificador;                //guarda su posicion
                    for (int pivot = pos_variable_decision; pivot < m; pivot++){
                        if (matrix[pivot][pos_variable_decision]==1){     //Busca en su columna el pivote
                            for (int mover = 0; mover < n; mover++){
                                double resultado = matrix[pivot][mover]*-numero_a_multiplicar+matrix[0][mover];  //hace la operacion correspondiente
                                bool sendocero= __FLT_EPSILON__>= fabs(resultado);
                                bool sendocero2= numeric_limits<double>::epsilon() >= fabs(resultado);
                                if (sendocero==1 or sendocero2==1){ 
                                    matrix[0][mover]= 0;
                                } else {
                                    matrix[0][mover] = resultado;
                                }
                            }
                        }
                    }
                }
            }
        } else if (modo == 0) {        //Modo 2 aplica el metodo de gauss jordan para que las variables de desicion no esten en 0 y sean negativos
            for (int contador = 0; contador < n; contador++){
                for (int k = 1; k < m; k++) {
                    matrix[0][contador] = matrix[k][contador]*-numero_a_multiplicar+matrix[0][contador];
                }
            }
            cout << "////////////////////////////////////////" << endl;
        }
    } else if (mixto == true){
        if (modo == 1) {             //Modo 1 aplica el metodo de gauss jordan para igualar el valor de las variables a 0
            for (int verificador = 0; verificador < n; verificador++){
                if (matrix[0][verificador] == numero_a_multiplicar) {   //Si encuentra un valor que es igual que el parametro dado (que se quiere igualar a 0)
                //if (matrix[0][verificador] < matrix[0][verificador+1] or matrix[0][verificador] > matrix[0][verificador-1] ) {
                //cout << matrix[0][verificador] << " " << matrix[0][verificador+1] << " " << matrix[0][verificador-1] << endl; 
                    pos_variable_decision = verificador;                //guarda su posicion
                    for (int pivot = 0; pivot < m; pivot++){
                        if (matrix[pivot][verificador]==1){     //Busca en su columna el pivote
                            for (int mover = 0; mover < n; mover++){
                                double resultado = matrix[pivot][mover]*-numero_a_multiplicar+matrix[0][mover];  //hace la operacion correspondiente
                                bool sendocero= __FLT_EPSILON__>= fabs(resultado);
                                bool sendocero2= numeric_limits<double>::epsilon() >= fabs(resultado);
                                if (sendocero==1 or sendocero2==1){ 
                                    matrix[0][mover]= 0;
                                } else {
                                    matrix[0][mover] = resultado;
                                }
                            }
                        }
                    }
                }
            }
        } else {        //Modo 2 aplica el metodo de gauss jordan para que las variables de desicion no esten en 0 y sean negativos
            for (int contador = 0; contador < n; contador++){
                for (int k = 1; k <= acumulador; k++) {
                    matrix[0][contador] = matrix[k][contador]*-numero_a_multiplicar+matrix[0][contador];
                }
            }
            cout << "////////////////////////////////////////" << endl;
        }
    }
    return matrix;
}

vector<vector<double>>& eliminarColumna(vector<vector<double>>& matrix, int columnaAEliminar) {
    int filas = matrix.size();
    int columnas = matrix[0].size();
    for (int i = 0; i < filas; i++) {
        matrix[i].erase(matrix[i].begin() + columnaAEliminar); // Eliminar el elemento en la columna a eliminar
    }
    return matrix;
}

vector<vector<double>> simplexaumentadotable(vector<vector<double>>& tableau, vector<double>& coeficientes_F_O, vector<double>& recursos_restricciones, int n_variables_decision_heredado, int acumulador, int acumulador_inverso , bool mixto) {
    int n_restricciones = tableau.size();  // Numero de restricciones M
    int n_variables_decision = tableau[0].size();  // Numero de variables de desicion N
    
    if (mixto == true){
        cout << "Metodo Aumentado mixto - Agregar variables de holgura" << endl;
        //Crea una tabla inicial con las dimenciones del numero de restricciones y el numero de variables de decision
        vector<vector<double>> tableau2(n_restricciones, vector<double>(n_variables_decision_heredado + n_restricciones + acumulador));

        for (int a = 0; a <= n_restricciones-1; a++ ){
            for (int b = 0; b <= n_variables_decision-1; b++){
                tableau2[a][b]= tableau[a][b];
            }
        }


        for (int d = acumulador_inverso; d < n_restricciones; d++) {
            //Inserta un 1 en la interseccion de las variables de holgura con la restriccion del comentado anterior
            tableau2[d][n_variables_decision - acumulador_inverso - 1 + d] = 1;
            //Inserta el lado derecho de las restricciones con la limitante anterior
            //tableau2[d][n_variables_decision + n_restricciones-2] = recursos_restricciones[i];
            
        }
        for (int i = 1; i < tableau2.size(); i++) {
            //Inserta el lado derecho de las restricciones con la limitante anterior
            tableau2[i][n_variables_decision + acumulador_inverso-1] = recursos_restricciones[i-1];
        }

        cout << "Tabla inicial:" << endl;
        printMatrix(tableau2);

        //simplexartificial(vector<double>& coeficientes_F_O, vector<vector<double>>& coeficientes_restricciones, vector<double>& recursos_restricciones, bool is_maximization, bool mixto, int acumulador, int acumulador_inverso);
        return tableau2;

    } else {

    //La explicacion es la misma que las anteriores, solo que se llama al afuncion eliminar columna y se aplica el modo 2 de gauss jordan
    cout << "Metodo M2F" << endl;
    cout << "Tabla inicial M2F Fase 2:" << endl;

    //eliminar columnas artificiales
    for (int g = n_variables_decision_heredado+1; g < n_variables_decision-acumulador-1; g++){
        eliminarColumna(tableau, g);
        n_variables_decision--;
    }

    //incorpora los valores de la funcion objetivo
    for (int h = 0; h < n_variables_decision_heredado; h++){
        tableau[0][h]= coeficientes_F_O[h];
    }

    printMatrix(tableau);

    cout << "Tabla procesada con Gauss Jordan:" << endl;

    //CORREGIR ESTO!! ESTO SE DEBE DETECTAR CUANTAS VARIABLES BASICAS ESTAN EN LAS DECISION (X1,X2) Y ELEGIR ESAS FILAS, NO POR EL ACUMULADOR
    if (!acumulador == 0 ){
        for (int e = 0; e < n_variables_decision_heredado-acumulador; e++){
            GaussJordan(tableau,tableau[0][e],1, false, 0);
            printMatrix(tableau);
        }
        cout << "Metodo Aumentado" << endl;
        while (true) {
            int pos_variable_decision = 0;
            for (int d = 1; d < n_variables_decision-1; d++) {
                if (tableau[0][d] < tableau[0][pos_variable_decision]) {
                    pos_variable_decision = d;
                }
            }
            if (tableau[0][pos_variable_decision] >= 0) {
                break;
            }
            int pos_fila_pivot = -1;
            double coeficiente_pivot = numeric_limits<double>::max();
            for (int l = 1; l <= n_restricciones-1; l++) {
                if (tableau[l][pos_variable_decision] > 0) {
                    double coeficiente_pivot_t = tableau[l][n_variables_decision -1] / tableau[l][pos_variable_decision];
                    if (coeficiente_pivot_t < coeficiente_pivot) {
                        coeficiente_pivot = coeficiente_pivot_t;
                        pos_fila_pivot = l;
                    }
                }
            }

            if (pos_fila_pivot == -1) {
                throw runtime_error("El problema lineal tiene infinitas soluciones");
            }

            double pivot_actual = tableau[pos_fila_pivot][pos_variable_decision];
            for (int o = 0; o <= n_variables_decision -1; o++) {
                tableau[pos_fila_pivot][o] /= pivot_actual;
            }

            for (int p = 0; p <= n_restricciones-1; p++) {
                if (p != pos_fila_pivot) {
                    double valor_multiplicar = tableau[p][pos_variable_decision];
                    for (int q = 0; q <= n_variables_decision -1; q++) {
                        double superpivot = tableau[p][q] - valor_multiplicar * tableau[pos_fila_pivot][q];
                        bool sendocero= __FLT_EPSILON__>= fabs(superpivot);
                        bool sendocero2= numeric_limits<double>::epsilon() >= fabs(superpivot);
                        if (sendocero==1 or sendocero2==1){ 
                            tableau[p][q] = 0;
                        } else {
                            tableau[p][q] =superpivot;
                        }
                        
                    }
                }
            }
            cout << "pivot actual: " << pivot_actual << endl;
            cout << "pivot actual fila: " << pos_fila_pivot << endl;
            cout << "Menor Coeficiente: " << coeficiente_pivot << endl;
            cout << "/////////////////////////////////////////////////////" << endl;
            cout << endl;
            cout << "Tabla actual:" << endl;
            printMatrix(tableau);
        }
    } else {
        for (int e = 0; e < n_variables_decision_heredado; e++){
            GaussJordan(tableau,tableau[0][e],1, false, 0);
            printMatrix(tableau);
        }
    }
    vector<double> valor_variables_decision(n_variables_decision, 0);
    vector<double> valor_variables_basicas(n_restricciones, 0);
    double optimal_value = tableau[0][n_variables_decision + n_restricciones];
    for (int j = 0; j < n_variables_decision; j++) {
        bool variable_basica = false;
        int pos_variable_basica = -1;
        for (int i = 1; i < n_restricciones; i++) {
            if (tableau[i][j] == 1 && !variable_basica) {
                variable_basica = true;
                pos_variable_basica = i;
            } else if (tableau[i][j] != 0) {
                variable_basica = false;
                break;
            }
        }
        if (variable_basica) {
            valor_variables_decision[j] = tableau[pos_variable_basica][n_variables_decision + n_restricciones];
        }
    }
    for (int i = 1; i < n_restricciones; i++) {
        valor_variables_basicas[i - 1] = tableau[i][n_variables_decision + n_restricciones];
    }
    vector<double> max_valores(n_variables_decision + n_restricciones + n_restricciones);
    int pos_max_valor = -1;
    double max_val = -numeric_limits<double>::min();
    for (int k = 0; k < n_restricciones + n_variables_decision; k++) {
        if (tableau[0][k] > max_val) {
            max_val = tableau[0][k];
            pos_max_valor = k;
        }  
    }
    cout << "Solucion Optima:" << endl;
    for (int j = 0; j < n_variables_decision; j++) {
        cout << "Variable x" << j + 1 << ": " << valor_variables_decision[j] << endl;
    }
    cout << "Valor Optimo: " << optimal_value << endl;
    cout << "Shadow Prices (Dual Variables):" << endl;
    for (int i = 0; i < n_restricciones; i++) {
        cout << "Restricciones " << i + 1 << ": " << valor_variables_basicas[i] << endl;
    }
    cout << "Verdadero Precio Sombra:" << endl;
    cout << "Variable h o e en la columna " << pos_max_valor << ": " << max_val << endl;
    cout << endl;
    return tableau;
    }
}

vector<vector<double>> simplexaumentado(vector<double>& coeficientes_F_O, vector<vector<double>>& coeficientes_restricciones, vector<double>& recursos_restricciones, bool is_maximization, bool mixto, int acumulador, int acumulador_inverso) {
    int n_restricciones = coeficientes_restricciones.size();  // Numero de restricciones M
    int n_variables_decision = coeficientes_restricciones[0].size();  // Numero de variables N

    // verifica si el problema es de maximizacion o no para convertir los coeficientes en negativo
    if (!is_maximization) {
        for (auto& val : coeficientes_F_O) {
            val = -val;
        }
    }

    if (mixto == true){
        cout << "Metodo Aumentado mixto - Agregar variables de holgura" << endl;
        //Crea una tabla inicial con las dimenciones del numero de restricciones y el numero de variables de decision
        vector<vector<double>> tableau(n_restricciones + 1, vector<double>(n_variables_decision + n_restricciones + 1));

    
        for (int i = 0; i < n_restricciones; i++) {
            for (int j = 0; j < n_variables_decision; j++) {
                //Inserta los coeficientes de las restricciones sin tocar la fila de Z (el +1)
                tableau[i+1][j] = coeficientes_restricciones[i][j];
            }
            //Inserta un 1 en la interseccion de las variables de holgura con la restriccion del comentado anterior
            tableau[i+1][n_variables_decision + i] = 1;
            //Inserta el lado derecho de las restricciones con la limitante anterior
            tableau[i+1][n_variables_decision + n_restricciones] = recursos_restricciones[i];
        }

        for (int j = 0; j < n_variables_decision; j++) {
            //Inserta los coeficientes de la funcion objetivo
            tableau[0][j] = coeficientes_F_O[j];
        }

        cout << "Tabla inicial:" << endl;
        printMatrix(tableau);

        //simplexartificial(vector<double>& coeficientes_F_O, vector<vector<double>>& coeficientes_restricciones, vector<double>& recursos_restricciones, bool is_maximization, bool mixto, int acumulador, int acumulador_inverso);
        return tableau;
    } else {

    
    cout << "Metodo Aumentado" << endl;
    //Crea una tabla inicial con las dimenciones del numero de restricciones y el numero de variables de decision
    vector<vector<double>> tableau(n_restricciones + 1, vector<double>(n_variables_decision + n_restricciones + 1));

    
    for (int i = 0; i < n_restricciones; i++) {
        for (int j = 0; j < n_variables_decision; j++) {
            //Inserta los coeficientes de las restricciones sin tocar la fila de Z (el +1)
            tableau[i+1][j] = coeficientes_restricciones[i][j];
        }
        //Inserta un 1 en la interseccion de las variables de holgura con la restriccion del comentado anterior
        tableau[i+1][n_variables_decision + i] = 1;
        //Inserta el lado derecho de las restricciones con la limitante anterior
        tableau[i+1][n_variables_decision + n_restricciones] = recursos_restricciones[i];
    }

    for (int j = 0; j < n_variables_decision; j++) {
        //Inserta los coeficientes de la funcion objetivo
        tableau[0][j] = coeficientes_F_O[j];
    }

    cout << "Tabla inicial:" << endl;
    printMatrix(tableau);

    while (true) {
        //Inicializa la posicion de la variable de decicion
        int pos_variable_decision = 0;

        //Busca y selecciona a la menor variable de decicion!
        for (int j = 1; j < n_variables_decision + n_restricciones; j++) {
            //Cuando se recorre el valor de las variables de decision, si el segundo valor es menor que el primero, guarda la posicion "j" en "pos"
            if (tableau[0][j] < tableau[0][pos_variable_decision]) {
                pos_variable_decision = j;
            }
        }

        //Y si los valores ya son 0 o mayores termina!
        if (tableau[0][pos_variable_decision] >= 0.0) {
            break;
        }

        //Se inicializa la posicion de la fila pivote
        int pos_fila_pivot = -1; 

        //inicializa con el valor maximo que puede tomar un double la variable que guardara el menor coeficiciente del pivote par ahacer comparaciones
        double coeficiente_pivot = numeric_limits<double>::max();


        //Idea!!, se puede agregar una lista en el primer if que guarde todos los coeficientes para despues motrarlos en pantalla
        //ya que divide directamente y no lo guarda despues
        
        //Busca el pivote actual seleccionando el menor de los coeficientes de los pivotes validos!
        for (int i = 1; i <= n_restricciones; i++) {
            //Si los valores de las variables básicas son mayores de 0 de la columna pivote seleccionada
            if (tableau[i][pos_variable_decision] > 0) {
                //Se calcula el coeficiente resultante de dividir el lado derecho de estas filas por los valores de la columna pivote
                double coeficiente_pivot_t = tableau[i][n_variables_decision + n_restricciones ] / tableau[i][pos_variable_decision];
                //Verifica si el coeficiente temporal sacado anteriormente es menor que el coeficiente inicializado
                if (coeficiente_pivot_t < coeficiente_pivot) {
                    //Si fuera asi lo guarda definitivamente
                    coeficiente_pivot = coeficiente_pivot_t;
                    //la posicion de la fila pivot es i
                    pos_fila_pivot = i;
                }
            }
        }

        //Si la posicion de la fila pivot es -1 finaliza el programa
        if (pos_fila_pivot == -1) {
            throw runtime_error("El problema lineal tiene infinitas soluciones");
        }

        //Se inicializa el pivot actual (es la interseccion de la fila de la v. de decicion y la fila pivot)
        double pivot_actual = tableau[pos_fila_pivot][pos_variable_decision];

        //Se divide el pivot asi mismo para que valga 1 y a toda la fila por el pivot!
        for (int j = 0; j <= n_variables_decision + n_restricciones ; j++) {
            tableau[pos_fila_pivot][j] /= pivot_actual;
        }

        //Busca que el resto de los valores de la fila de la variable de decicion sea 0!
        for (int i = 0; i <= n_restricciones; i++) {
            //Si la fila no es la que esta el pivot
            if (i != pos_fila_pivot) {
                //Se instancia el valor que no es el pivot pero que estan en la fila de la variable de decicion
                double valor_multiplicar = tableau[i][pos_variable_decision];
                for (int j = 0; j <= n_variables_decision + n_restricciones; j++) {
                    //Hace el calculo correspondiente
                    tableau[i][j] -= valor_multiplicar * tableau[pos_fila_pivot][j];
                }
            }
        }
        cout << "Tabla Actual:" << endl;
        printMatrix(tableau);
    }

    //En pruebas!!
    //Se declaran las listas de la solucion optima (solo las variables de decicion que estan como basicas) y de las variables basicas, y el precio sombra 
    //y el valor de cada variable basica
    //Hay que cambiar en que solo en una lista los guarde y que tambien lo haga el nombre.
    vector<double> valor_variables_decision(n_variables_decision, 0);
    vector<double> valor_variables_basicas(n_restricciones, 0);
    double optimal_value = tableau[0][n_variables_decision + n_restricciones];

    //Encuentra el valor de las variables basicas!
    for (int j = 0; j < n_variables_decision; j++) {
        //misma idea cuando se busca la variable de decision menor
        bool variable_basica = false;
        int pos_variable_basica = -1;
        for (int i = 1; i <= n_restricciones; i++) {
            //cuando recorre la tabla y encuentra un 1 y no es una variable basica
            if (tableau[i][j] == 1 && !variable_basica) {
                //se declara como variable basica y guarda su posicion
                variable_basica = true;
                pos_variable_basica = i;
            } else if (tableau[i][j] != 0) {
                //sino si el valor de la tabla no es 0 se declara como no es una variable basica y termina
                variable_basica = false;
                break;
            }
        }
        // y si es una variable basica se guarda en su lista correspondiente, HAY QUE PROBAR BIEN CON VARIOS EJEMPLOS SI ARROJA BIEN SOLO LAS 
        //VARIABLES DE DECISION!!
        if (variable_basica) {
            valor_variables_decision[j] = tableau[pos_variable_basica][n_variables_decision + n_restricciones];
        }
    }

    // guarda el valor de los lados derechos de las restricciones y las guarda en la lista de las variables basicas
    //falta implementar el nombre (h1,h2,etc)
    for (int i = 1; i <= n_restricciones; i++) {
        valor_variables_basicas[i - 1] = tableau[i][n_variables_decision + n_restricciones];
    }

    //Se declara la lista que guardara el valor maximo del coeficiente de las variables de la funcion objetivo
    vector<double> max_valores(n_variables_decision + n_restricciones + n_restricciones);
    int pos_max_valor = -1;
    double max_val = -numeric_limits<double>::min();

    //si se encuentra un valor en la primera fila mayor que el valor predefinido se convierte en su actual valor maximo y guarda la posicion
    for (int k = 0; k < n_restricciones + n_variables_decision; k++) {
        if (tableau[0][k] > max_val) {
            max_val = tableau[0][k];
            pos_max_valor = k;
        }  
    }

    cout << "Solucion Optima:" << endl;
    for (int j = 0; j < n_variables_decision; j++) {
        cout << "Variable x" << j + 1 << ": " << valor_variables_decision[j] << endl;
    }
    cout << "Valor Optimo de Z: " << optimal_value << endl;
    cout << "Shadow Prices (Dual Variables):" << endl;
    for (int i = 0; i < n_restricciones; i++) {
        cout << "Restricciones " << i + 1 << ": " << valor_variables_basicas[i] << endl;
    }
    cout << "Verdadero Precio Sombra:" << endl;
    cout << "Variable h o e en la columna " << pos_max_valor << ": " << max_val << endl;
    cout << endl;
    return tableau;
    }
}

vector<vector<double>> simplexartificial(vector<double>& coeficientes_F_O, vector<vector<double>>& coeficientes_restricciones, vector<double>& recursos_restricciones, bool is_maximization, bool mixto, int acumulador, int acumulador_inverso) {
    int n_restricciones = coeficientes_restricciones.size();  // Numero de restricciones
    int n_variables_decision = coeficientes_restricciones[0].size();  // Numero de variables de decision

    //El metodo es casi el mismo, solo que no se insertan los coeficientes de la f.o y antes de iniciar el pivoteo se aplica gauss jordan
    if (!is_maximization) {
        for (auto& val : coeficientes_F_O) {
            val = -val;
        }
    }

    if (mixto == true) {
        cout << "Metodo M2F mixto - Agregar variables artificiales y de exceso" << endl;
        vector<vector<double>> tableau(n_restricciones + 1, vector<double>(n_variables_decision + n_restricciones));
        for (int i = 0; i < n_restricciones; i++) {
            for (int j = 0; j < n_variables_decision; j++) {
                tableau[i + 1][j] = coeficientes_restricciones[i][j];
            }
            //tableau[i + 1][n_variables_decision+i+i] = -1.0;
            //tableau[i + 1][n_restricciones+i+i] = 1.0;
        }

        for (int i = 0; i < acumulador; i++) {
            tableau[i + 1][n_variables_decision+i+i] = -1.0;
            tableau[i + 1][n_restricciones+i+i] = 1.0;
        }

        for (int j = 0; j < acumulador; j++) {
            tableau[0][n_restricciones+j*2] = 1.0;
        }
        
        cout << "Tabla inicial:" << endl;
        printMatrix(tableau);
        return tableau;
        
    } else {
    //Diferencia principal, para rellenar la matriz correctamente se le suma 2 veces i para que coloque correctamente los 1 y -1
        cout << "Metodo M2F" << endl;
        vector<vector<double>> tableau(n_restricciones + 1, vector<double>(n_variables_decision + n_restricciones + n_restricciones + 1));
        for (int i = 0; i < n_restricciones; i++) {
            for (int j = 0; j < n_variables_decision; j++) {
                tableau[i + 1][j] = coeficientes_restricciones[i][j];
            }
            tableau[i + 1][n_variables_decision+i+i] = -1.0;
            tableau[i + 1][n_restricciones+i+i] = 1.0;
            tableau[i + 1][n_variables_decision + n_restricciones + n_restricciones] = recursos_restricciones[i];
        }

        for (int j = 0; j < n_restricciones; j++) {
            tableau[0][n_restricciones+j*2] = 1.0;
        }
        
        cout << "Tabla inicial M2F Fase 1:" << endl;
        printMatrix(tableau);

        tableau = GaussJordan(tableau,1,0, false, 0);

        cout << "Tabla inicial procesada con Gauss:" << endl;
        printMatrix(tableau);

        //Otra diferencia esque se suma 2 veces el numero de restricciones, ya que son 2 veces que se agregan variables, las artificiales y de exceso
        while (true) {
            int pos_variable_decision = 0;
            for (int d = 1; d < n_variables_decision + n_restricciones + n_restricciones; d++) {
                if (tableau[0][d] < tableau[0][pos_variable_decision]) {
                    pos_variable_decision = d;
                }
            }
            if (tableau[0][pos_variable_decision] >= 0) {
                break;
            }
            int pos_fila_pivot = -1;
            double coeficiente_pivot = numeric_limits<double>::max();
            for (int l = 1; l <= n_restricciones; l++) {
                if (tableau[l][pos_variable_decision] > 0) {
                    double coeficiente_pivot_t = tableau[l][n_variables_decision + n_restricciones + n_restricciones] / tableau[l][pos_variable_decision];
                    if (coeficiente_pivot_t < coeficiente_pivot) {
                        coeficiente_pivot = coeficiente_pivot_t;
                        pos_fila_pivot = l;
                    }
                }
            }

            if (pos_fila_pivot == -1) {
                throw runtime_error("El problema lineal tiene infinitas soluciones");
            }

            double pivot_actual = tableau[pos_fila_pivot][pos_variable_decision];
            for (int o = 0; o <= n_variables_decision + n_restricciones + n_restricciones; o++) {
                tableau[pos_fila_pivot][o] /= pivot_actual;
            }

            for (int p = 0; p <= n_restricciones; p++) {
                if (p != pos_fila_pivot) {
                    double valor_multiplicar = tableau[p][pos_variable_decision];
                    for (int q = 0; q <= n_variables_decision + n_restricciones + n_restricciones; q++) {
                        double superpivot = tableau[p][q] - valor_multiplicar * tableau[pos_fila_pivot][q];
                        bool sendocero= __FLT_EPSILON__>= fabs(superpivot);
                        bool sendocero2= numeric_limits<double>::epsilon() >= fabs(superpivot);
                        if (sendocero==1 or sendocero2==1){ 
                            tableau[p][q] = 0;
                        } else {
                            tableau[p][q] =superpivot;
                        }
                        
                    }
                }
            }
            cout << "pivot actual: " << pivot_actual << endl;
            cout << "pivot actual fila: " << pos_fila_pivot << endl;
            cout << "Menor Coeficiente: " << coeficiente_pivot << endl;
            cout << "/////////////////////////////////////////////////////" << endl;
            cout << endl;
            cout << "Tabla actual:" << endl;
            printMatrix(tableau);
        }
        simplexaumentadotable(tableau, coeficientes_F_O, recursos_restricciones, n_variables_decision, 0, 0, false);
        return tableau;
        }
    }

vector<vector<double>> simplexartificialtable(vector<vector<double>>& tableau, vector<double>& coeficientes_F_O, vector<double>& recursos_restricciones, int n_variables_decision_heredado, int acumulador, int acumulador_inverso) {
    int n_restricciones = tableau.size();  // Numero de restricciones M
    int n_variables_decision = tableau[0].size();  // Numero de variables de desicion N

    //Diferencia principal, para rellenar la matriz correctamente se le suma 2 veces i para que coloque correctamente los 1 y -1
        cout << "Metodo M2F" << endl;
        cout << "Tabla inicial M2F Fase 1:" << endl;
        printMatrix(tableau);

        tableau = GaussJordan(tableau,1,0,true, acumulador);

        cout << "Tabla inicial procesada con Gauss:" << endl;
        printMatrix(tableau);

        //Otra diferencia esque se suma 2 veces el numero de restricciones, ya que son 2 veces que se agregan variables, las artificiales y de exceso
        while (true) {
            int pos_variable_decision = 0;
            for (int d = 1; d < n_variables_decision-1; d++) {
                if (tableau[0][d] < tableau[0][pos_variable_decision]) {
                    pos_variable_decision = d;
                }
            }
            if (tableau[0][pos_variable_decision] >= 0) {
                break;
            }
            int pos_fila_pivot = -1;
            double coeficiente_pivot = numeric_limits<double>::max();
            for (int l = 1; l <= n_restricciones-1; l++) {
                if (tableau[l][pos_variable_decision] > 0) {
                    double coeficiente_pivot_t = tableau[l][n_variables_decision -1] / tableau[l][pos_variable_decision];
                    if (coeficiente_pivot_t < coeficiente_pivot) {
                        coeficiente_pivot = coeficiente_pivot_t;
                        pos_fila_pivot = l;
                    }
                }
            }

            if (pos_fila_pivot == -1) {
                throw runtime_error("El problema lineal tiene infinitas soluciones");
            }

            double pivot_actual = tableau[pos_fila_pivot][pos_variable_decision];
            for (int o = 0; o <= n_variables_decision -1; o++) {
                tableau[pos_fila_pivot][o] /= pivot_actual;
            }

            for (int p = 0; p <= n_restricciones-1; p++) {
                if (p != pos_fila_pivot) {
                    double valor_multiplicar = tableau[p][pos_variable_decision];
                    for (int q = 0; q <= n_variables_decision -1; q++) {
                        double superpivot = tableau[p][q] - valor_multiplicar * tableau[pos_fila_pivot][q];
                        bool sendocero= __FLT_EPSILON__>= fabs(superpivot);
                        bool sendocero2= numeric_limits<double>::epsilon() >= fabs(superpivot);
                        if (sendocero==1 or sendocero2==1){ 
                            tableau[p][q] = 0;
                        } else {
                            tableau[p][q] =superpivot;
                        }
                        
                    }
                }
            }
            cout << "pivot actual: " << pivot_actual << endl;
            cout << "pivot actual fila: " << pos_fila_pivot << endl;
            cout << "columna variable de decision: " << pos_variable_decision << endl;
            cout << "Menor Coeficiente: " << coeficiente_pivot << endl;
            cout << "/////////////////////////////////////////////////////" << endl;
            cout << endl;
            cout << "Tabla actual:" << endl;
            printMatrix(tableau);
        }
        simplexaumentadotable(tableau, coeficientes_F_O, recursos_restricciones, n_variables_decision_heredado, acumulador, acumulador_inverso, false);
        return tableau;
        }
    

void verificador(vector<double> coeficientes_F_O, vector<vector<double>> coeficientes_restricciones, vector<double> recursos_restricciones, vector<string> signos, bool is_maximization){
    cout << "Verificador de metodologia: " << endl;
    //Verifica el signo de las restricciones y ejecuta el programa correspondiente
    int acumulador = 0;
    bool mixto;
    //recorre toda la lista de signo
    for (int a=0; a < signos.size(); a++)
        if (signos[0] == signos[a]){
            acumulador++;
        }
    int acumulador_inverso = signos.size()-acumulador;
    cout << "Hay " << acumulador << " Restricciones de tipo " << signos[0] << " y " << signos.size()-acumulador << " del otro tipo" << endl;
    if (acumulador == signos.size()){
        mixto = false;
        if (signos[0] == "<="){
            cout << "-Programa procede a ejecutarse en modo aumentado" << endl;
            simplexaumentado(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_maximization, mixto, 0, 0);
        } else {
            cout << "-Programa procede a ejecutarse en modo M2F" << endl;
            simplexartificial(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_maximization, mixto, 0, 0);
        }
    } else {
        mixto = true;
        //No listo
        if (signos[0] == "<="){
            cout << "-Programa procede a ejecutarse en primera instancia en modo aumentado" << endl;
            //Creacion de la tabla
            vector<vector<double>> aumentadomixto = simplexaumentado(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_maximization, mixto, acumulador, acumulador_inverso);
            vector<vector<double>> artificialmixto = simplexartificial(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_maximization, mixto, acumulador, acumulador_inverso);
            //Resolucion del metodo
            //simplexartificialtable(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_maximization, mixto, acumulador, acumulador_inverso);
        //listo
        } else {
            cout << "-Programa procede a ejecutarse en primera instancia en modo M2F" << endl;
            //Creacion de la tabla
            vector<vector<double>> artificialmixto2 = simplexartificial(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_maximization, mixto, acumulador, acumulador_inverso);
            vector<vector<double>> aumentadomixto2 = simplexaumentadotable(artificialmixto2, coeficientes_F_O, recursos_restricciones, coeficientes_F_O.size(), acumulador, acumulador_inverso, mixto);

            //Resolucion del metodo
            simplexartificialtable(aumentadomixto2, coeficientes_F_O, recursos_restricciones, coeficientes_F_O.size(), acumulador, acumulador_inverso);
        }
    }
    
}

int main() {
    // Problema Fabrica de productos de vidrio
    vector<double> recursos_restricciones1 = {4, 12, 18};
    vector<double> coeficientes_F_O1 = {-30000, -50000};
    vector<vector<double>> coeficientes_restricciones1= {{1, 0}, {0, 2}, {3, 2}};
    vector<string> signos1 = {"<=", "<=", "<="};

    // Problema La Perfumeria
    //b = valor de las restricciones
    vector<double> recursos_restricciones2 = {100, 80, 45, 100};

    //c = coeficientes de la funcion objetivo
    vector<double> coeficientes_F_O2 = {-2, -2};

    //A = coeficientes de las variables de decicion (restricciones)
    vector<vector<double>> coeficientes_restricciones2 = {{2, 1}, {1, 3}, {1, 0}, {0, 1}};    
    //boleano 1 si es maximizar, boleano 2 si es artificial
    vector<string> signos2 = {"<=", "<=", "<=", "<="};

    //Problema Petroleo
    vector<double> coeficientes_F_O3 = {-0.12, -0.15};
    vector<double> recursos_restricciones3 = {300,36,90};
    vector<vector<double>> coeficientes_restricciones3 = {{60, 60}, {12,6}, {10, 30}};
    vector<string> signos3 = {">=", ">=", ">="};

    vector<double> coeficientes_F_O4 = {-3, 2};
    vector<double> recursos_restricciones4 = {18,42,5};
    vector<vector<double>> coeficientes_restricciones4 = {{2,1}, {2,3}, {3, -2}};
    vector<string> signos4 = {"<=", "<=", "<="};

    //problema multivariable
    vector<double> coeficientes_F_O5 = {-1500, -1400, -1600, -1450};
    vector<double> recursos_restricciones5 = {40, 70, 0, 180, 45};
    vector<vector<double>> coeficientes_restricciones5 = {{1, 0, 1, 0}, {0, 1, 0, 1}, {2, -1, 2, -1}, {1, 1, 0, 0}, {0, 0, 1, 1}};
    vector<string> signos5 = {">=", ">=", "<=", "<=", "<="};
    

    //Simplex Artificial solo aplica con restricciones de mayor o igual
    //verificador(coeficientes_F_O1, coeficientes_restricciones1, recursos_restricciones1, signos1, true);
    //verificador(coeficientes_F_O2, coeficientes_restricciones2, recursos_restricciones2, signos2, true);
    //verificador(coeficientes_F_O3, coeficientes_restricciones3, recursos_restricciones3, signos3, false);
    //verificador(coeficientes_F_O4, coeficientes_restricciones4, recursos_restricciones4, signos4, true);
    verificador(coeficientes_F_O5, coeficientes_restricciones5, recursos_restricciones5, signos5, true);

    system("pause");
    return 0;
}