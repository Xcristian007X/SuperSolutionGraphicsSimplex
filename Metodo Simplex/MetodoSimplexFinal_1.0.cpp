#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <regex>
#include <algorithm>
#include <string>
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

void printMatrix2(vector<vector<double>>& matrix,vector<string>& coll,vector<string>&fill) {
    int cont=0;
    for(int i=0;i<coll.size();i++){
        cout << coll[i]<<"\t";
    }
    cout << endl;
    for (const auto& row : matrix) {
        cout << fill[cont]<<"\t";
        cont=cont+1;
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

vector<vector<double>> simplexaumentadotable(vector<vector<double>>& tableau, vector<double>& coeficientes_F_O, vector<double>& recursos_restricciones, int n_variables_decision_heredado, int acumulador, int acumulador_inverso , bool mixto,vector<string>col,vector<string>fil) {
    int n_restricciones = tableau.size();  // Numero de restricciones M
    int n_variables_decision = tableau[0].size();  // Numero de variables de desicion N
    int cont=0;
    int c=2;
    int aumentador=0;
    if (mixto == true){
        cout << "Metodo Aumentado mixto - Agregar variables de holgura" << endl;
        //Crea una tabla inicial con las dimenciones del numero de restricciones y el numero de variables de decision
        vector<vector<double>> tableau2(n_restricciones, vector<double>(n_variables_decision + acumulador_inverso+1));

        for(int i=0;i<col.size();i++){
            if(i==0){col[i]="VB";}
            if(i==1){col[i]="Z";}
            if(i==col.size()-1){col[i]="LD";}
            if(i>1 && i<col.size()-3){
                if(cont <n_variables_decision_heredado){
                    col[i]="x"+to_string(cont+1);
                    cont=cont+1;
                    c=i+2;
                }else if(aumentador<acumulador){
                    col[c-1]="e"+to_string(aumentador+1);
                    col[c]="a"+to_string(aumentador+1);
                    c=c+2;
                    aumentador=aumentador+1;
                }else{
                    col[c-1]="h"+to_string(aumentador+1);
                    c=c+1;
                    aumentador=aumentador+1;
                }
            }
        }
        aumentador=0;
        //Llenado de vector para nombres de las filas
        for(int i=0;i<fil.size();i++){
            if(i==0){fil[i]="Z";}
            if(i>0){
                if(aumentador<acumulador){
                    fil[i]="a"+to_string(aumentador+1);
                    aumentador=aumentador+1;
                }else{
                    fil[i]="h"+to_string(aumentador+1);
                    aumentador=aumentador+1;
                }
            }
        }
        
        

        for (int a = 0; a <= n_restricciones-1; a++ ){
            for (int b = 0; b <= n_variables_decision-1; b++){
                tableau2[a][b+1]= tableau[a][b];
            }
        }


        for (int d = acumulador+1; d < n_restricciones; d++) {
            //Inserta un 1 en la interseccion de las variables de holgura con la restriccion del comentado anterior
            tableau2[d][n_variables_decision_heredado + acumulador -1 + d + 1] = 1;
            //Inserta el lado derecho de las restricciones con la limitante anterior
            //tableau2[d][n_variables_decision + n_restricciones-2] = recursos_restricciones[i];
            
        }
        for (int i = 1; i < tableau2.size(); i++) {
            //Inserta el lado derecho de las restricciones con la limitante anterior
            tableau2[i][n_variables_decision + acumulador_inverso] = recursos_restricciones[i-1];
        }

        cout << "Tabla inicial:" << endl;
        printMatrix2(tableau2,col,fil);

        //simplexartificial(vector<double>& coeficientes_F_O, vector<vector<double>>& coeficientes_restricciones, vector<double>& recursos_restricciones, bool is_minimize, bool mixto, int acumulador, int acumulador_inverso);
        return tableau2;

    } else {

    //La explicacion es la misma que las anteriores, solo que se llama al afuncion eliminar columna y se aplica el modo 2 de gauss jordan
    cout << "Metodo M2F" << endl;
    cout << "Tabla inicial M2F Fase 2:" << endl;

    //eliminar columnas artificiales
    for (int g = n_variables_decision_heredado+2; g < n_variables_decision-acumulador_inverso-1; g++){
        eliminarColumna(tableau, g);
        n_variables_decision--;
        col.erase(col.begin() + g+1);
    }

    //incorpora los valores de la funcion objetivo
    for (int h = 0; h < n_variables_decision_heredado; h++){
        tableau[0][h+1]= coeficientes_F_O[h];
    }

    printMatrix2(tableau,col,fil);

    cout << "Tabla procesada con Gauss Jordan:" << endl;

    //CORREGIR ESTO!! ESTO SE DEBE DETECTAR CUANTAS VARIABLES BASICAS ESTAN EN LAS DECISION (X1,X2) Y ELEGIR ESAS FILAS, NO POR EL ACUMULADOR

    // if (fil [0][0] and col [0][1] and =tableau ==1)
/*     for(int i = 0; i< fil.size();i++){
        if(fil[i+1] == col[0])
    } */

    vector<string> verificador_filas;
    for(int i = 2; i <col.size();i++){
        if(col[i][0] =='x'){
            verificador_filas.push_back(col[i]);
        }
    }
    vector<string> verificador_columnas;
    for(int i = 1; i <fil.size();i++){
        if(col[i][0] == 'x'){
            verificador_columnas.push_back(col[i]);
        }
    }
    bool son_iguales = verificador_filas.size() == verificador_columnas.size() && std::equal(verificador_filas.begin(), verificador_filas.end(), verificador_columnas.begin());

    if (son_iguales) {
        if (!acumulador == 0 ){
        for (int e = 0; e < n_variables_decision_heredado-acumulador+1; e++){
            GaussJordan(tableau,tableau[0][e],1, false, 0);
            printMatrix2(tableau,col,fil);
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
            cout << "Tabla actual:" << endl;
            fil[pos_fila_pivot]=col[pos_variable_decision+1];
            printMatrix2(tableau,col,fil);
        }
    }
        std::cout << "La primera columna coincide con la primera fila en los valores de x.\n";
    } else {
        for (int e = 0; e < n_variables_decision_heredado+1; e++){
            GaussJordan(tableau,tableau[0][e],1, false, 0);
            printMatrix2(tableau,col,fil);;
        }
        std::cout << "La primera columna no coincide con la primera fila en los valores de x1, x2, x3, ..., xn.\n";
    }

    

    //Mostrar soluciones
    double optimal_value = tableau[0][n_variables_decision_heredado + n_restricciones];
    //Se declara la lista que guardara el valor maximo del coeficiente de las variables de la funcion objetivo
    vector<double> max_valores(n_variables_decision_heredado + n_restricciones + n_restricciones);
    int pos_max_valor = -1;
    double max_val = -numeric_limits<double>::min();

    //si se encuentra un valor en la primera fila mayor que el valor predefinido se convierte en su actual valor maximo y guarda la posicion
    for (int k = 1; k < n_restricciones + n_variables_decision_heredado+1; k++) {
        if (tableau[0][k] > max_val) {
            max_val = tableau[0][k];
            pos_max_valor = k-1;
        }  
    }
    cout << "Solucion Optima:" << endl;
    cout << "Valor Optimo de Z: " << optimal_value << endl;
    bool comp=true;
    for(int x=2; x<n_restricciones+n_restricciones-2;x++){
        for(int e=1; e<n_restricciones+1;e++){
            if(col[x]==fil[e]){
                cout<<col[x]<<"="<<tableau[e][n_variables_decision_heredado+n_restricciones] <<endl;
                comp=false;
            }
        }
        if(comp){
            cout<<col[x]<<"="<<0 <<endl;
        }
        comp=true;
    }
    cout << endl;
    return tableau;

    return tableau;
    }
}

vector<vector<double>> simplexaumentado(vector<double>& coeficientes_F_O, vector<vector<double>>& coeficientes_restricciones, vector<double>& recursos_restricciones, bool is_minimize, int acumulador, int acumulador_inverso) {
    int n_restricciones = coeficientes_restricciones.size();  // Numero de restricciones M
    int n_variables_decision = coeficientes_restricciones[0].size();  // Numero de variables N
    
    // verifica si el problema es de maximizacion o no para convertir los coeficientes en negativo
    if (!is_minimize) {
        for (auto& val : coeficientes_F_O) {
            val = -val;
        }
    }
    int c=0;
    vector<string> fil(n_restricciones + 1);//Vector para almacenar nombres de las filas
    int cont=1;
    vector <string> col(n_variables_decision + n_restricciones + 3);//Vector para almacenas nombres de las columnas
    //Llenado de vector para nombres columnas.
    for(int i=0;i<col.size();i++){
        if(i==0){col[i]="VB";}
        if(i==1){col[i]="Z";}
        if(i==col.size()-1){col[i]="LD";}
        if(i>1 && i<col.size()-1){
            if(cont<=n_variables_decision){
                col[i]="x"+to_string(cont);
                cont=cont+1;
            }else{
                col[i]="h"+to_string(cont-n_variables_decision);
                cont=cont+1;
            }
        }
    }
    //Llenado de vector para nombres de las filas
    for(int i=0;i<fil.size();i++){
        if(i==0){fil[i]="Z";}
        if(i>0){
            fil[i]="h"+to_string(i);
        }
    }
    
    cout << "Metodo Aumentado" << endl;
    //Crea una tabla inicial con las dimenciones del numero de restricciones y el numero de variables de decision
    vector<vector<double>> tableau(n_restricciones + 1, vector<double>(n_variables_decision + n_restricciones + 2));

    for (int i = 0; i < n_restricciones; i++) {
        for (int j = 0; j < n_variables_decision; j++) {
            //Inserta los coeficientes de las restricciones sin tocar la fila de Z (el +1)
            tableau[i+1][j+1] = coeficientes_restricciones[i][j];
        }
        //Inserta un 1 en la interseccion de las variables de holgura con la restriccion del comentado anterior
        tableau[i+1][n_variables_decision + i+1] = 1;
        //Inserta el lado derecho de las restricciones con la limitante anterior
        tableau[i+1][n_variables_decision + n_restricciones+1] = recursos_restricciones[i];
    }

    for (int j = 0; j < n_variables_decision; j++) {
        //Inserta los coeficientes de la funcion objetivo
        tableau[0][j+1] = coeficientes_F_O[j];
    }

    cout << "Tabla inicial:" << endl;
    //printMatrix(tableau);
    printMatrix2(tableau,col,fil); 

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
        for (int i = 0; i <= n_restricciones; i++) {
            //Si los valores de las variables básicas son mayores de 0 de la columna pivote seleccionada
            if (tableau[i][pos_variable_decision] > 0) {
                //Se calcula el coeficiente resultante de dividir el lado derecho de estas filas por los valores de la columna pivote
                double coeficiente_pivot_t = tableau[i][n_variables_decision + n_restricciones+1 ] / tableau[i][pos_variable_decision];
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
        for (int j = 0; j <= n_variables_decision + n_restricciones +1; j++) {
            tableau[pos_fila_pivot][j] /= pivot_actual;
        }

        //Busca que el resto de los valores de la fila de la variable de decicion sea 0!
        for (int i = 0; i <= n_restricciones; i++) {
            //Si la fila no es la que esta el pivot
            if (i != pos_fila_pivot) {
                //Se instancia el valor que no es el pivot pero que estan en la fila de la variable de decicion
                double valor_multiplicar = tableau[i][pos_variable_decision];
                for (int j = 0; j <= n_variables_decision + n_restricciones+1; j++) {
                    //Hace el calculo correspondiente
                    tableau[i][j] -= valor_multiplicar * tableau[pos_fila_pivot][j];
                }
            }
        }
        fil[pos_fila_pivot]=col[pos_variable_decision+1];
        cout << "Tabla Actual:" << endl;
        printMatrix2(tableau,col,fil);
    }
    //Mostrar soluciones
    double optimal_value = tableau[0][n_variables_decision + n_restricciones+1];
    //Se declara la lista que guardara el valor maximo del coeficiente de las variables de la funcion objetivo
    vector<double> max_valores(n_variables_decision + n_restricciones + n_restricciones);
    int pos_max_valor = -1;
    double max_val = -numeric_limits<double>::min();

    //si se encuentra un valor en la primera fila mayor que el valor predefinido se convierte en su actual valor maximo y guarda la posicion
    for (int k = 1; k < n_restricciones + n_variables_decision+1; k++) {
        if (tableau[0][k] > max_val) {
            max_val = tableau[0][k];
            pos_max_valor = k-1;
        }  
    }
    cout << "Solucion Optima:" << endl;
    cout << "Valor Optimo de Z: " << optimal_value << endl;
    bool comp=true;
    for(int x=2; x<n_restricciones+n_restricciones ;x++){
        for(int e=1; e<n_restricciones+1;e++){
            if(col[x]==fil[e]){
                cout<<col[x]<<"="<<tableau[e][n_variables_decision+n_restricciones+1] <<endl;
                comp=false;
            }
        }
        if(comp){
            cout<<col[x]<<"="<<0 <<endl;
        }
        comp=true;
    }
    cout << "Valores precios sombra "<<endl;
    for(int i=n_variables_decision+2;i<col.size()-1;i++){
        cout<<col[i]<<"="<<tableau[0][i-1]<<endl;;
    }
    cout << endl;
    return tableau;
    }


vector<vector<double>> simplexartificial(vector<double>& coeficientes_F_O, vector<vector<double>>& coeficientes_restricciones, vector<double>& recursos_restricciones, bool is_minimize, bool mixto, int acumulador, int acumulador_inverso) {
    int n_restricciones = coeficientes_restricciones.size();  // Numero de restricciones
    int n_variables_decision = coeficientes_restricciones[0].size();  // Numero de variables de decision

    //El metodo es casi el mismo, solo que no se insertan los coeficientes de la f.o y antes de iniciar el pivoteo se aplica gauss jordan
    if (!is_minimize) {
        for (auto& val : coeficientes_F_O) {
            val = -val;
        }
    }

    if (mixto == true) {
        cout << "Metodo M2F mixto - Agregar variables artificiales y de exceso" << endl;
        vector<vector<double>> tableau(n_restricciones + 1, vector<double>(n_variables_decision + acumulador*2+1));
        for (int i = 0; i < n_restricciones; i++) {
            for (int j = 0; j < n_variables_decision; j++) {
                tableau[i + 1][j] = coeficientes_restricciones[i][j];
            }
        }

        for (int i = 0; i < acumulador; i++) {
            tableau[i + 1][n_variables_decision+i+i] = -1.0;
            tableau[i + 1][n_variables_decision+1+i+i] = 1.0;
            tableau[0][n_variables_decision+1+i*2] = 1.0;
        }

         
        cout << "Tabla inicial:" << endl;
        printMatrix(tableau);
        return tableau;
        
    } else {
        int c=1;//Variable para llenar variables artificiales y de exceso.
        int aumentador=1;//indice
        vector<string> fil(n_restricciones + 1);//Vector para almacenar nombres de las filas
        int cont=1;
        vector <string> col(n_variables_decision + n_restricciones + 6);//Vector para almacenas nombres de las columnas
        //Llenado de vector para nombres columnas.
        for(int i=0;i<col.size();i++){
            if(i==0){col[i]="VB";}
            if(i==1){col[i]="Z";}
            if(i==col.size()-1){col[i]="LD";}
            if(i>1 && i<col.size()-4){
                if(cont==n_variables_decision){c=i+2;}
                if(cont<=n_variables_decision){
                    col[i]="x"+to_string(cont);
                    cont=cont+1;
                }else{
                    col[c-1]="e"+to_string(aumentador);
                    col[c]="a"+to_string(aumentador);
                    cont=cont+1;
                    c=c+2;
                    aumentador=aumentador+1;
                }
            }
        }
        //Llenado de vector para nombres de las filas
        for(int i=0;i<fil.size();i++){
            if(i==0){fil[i]="Z";}
            if(i>0){
                fil[i]="a"+to_string(i);
            }
        }
    //Diferencia principal, para rellenar la matriz correctamente se le suma 2 veces i para que coloque correctamente los 1 y -1
        cout << "Metodo M2F" << endl;
        vector<vector<double>> tableau(n_restricciones + 1, vector<double>(n_variables_decision + n_restricciones*2 + 2));
        for (int i = 0; i < n_restricciones; i++) {
            for (int j = 0; j < n_variables_decision; j++) {
                tableau[i + 1][j+1] = coeficientes_restricciones[i][j];
            }
            tableau[i + 1][n_variables_decision+i+i+1] = -1.0;
            tableau[i + 1][n_restricciones+i+i+1] = 1.0;
            tableau[i + 1][n_variables_decision + n_restricciones + n_restricciones+1] = recursos_restricciones[i];
        }

        for (int j = 0; j < n_restricciones; j++) {
            tableau[0][n_restricciones+1+j*2] = 1.0;
        }
        
        cout << "Tabla inicial M2F Fase 1:" << endl;
        printMatrix2(tableau,col,fil);

        tableau = GaussJordan(tableau,1,0, false, 0);

        cout << "Tabla inicial procesada con Gauss:" << endl;
        printMatrix2(tableau,col,fil);

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
                    double coeficiente_pivot_t = tableau[l][n_variables_decision + n_restricciones + n_restricciones+1] / tableau[l][pos_variable_decision];
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
            for (int o = 0; o <= n_variables_decision + n_restricciones + n_restricciones+1; o++) {
                tableau[pos_fila_pivot][o] /= pivot_actual;
            }

            for (int p = 0; p <= n_restricciones; p++) {
                if (p != pos_fila_pivot) {
                    double valor_multiplicar = tableau[p][pos_variable_decision];
                    for (int q = 0; q <= n_variables_decision + n_restricciones + n_restricciones+1; q++) {
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
            cout << "Tabla actual:" << endl;
            fil[pos_fila_pivot]=col[pos_variable_decision+1];
            printMatrix2(tableau,col,fil);
        }
        simplexaumentadotable(tableau, coeficientes_F_O, recursos_restricciones, n_variables_decision, 0, 0, false,col,fil);
        return tableau;
        }
    }

vector<vector<double>> simplexartificialtable(vector<vector<double>>& tableau, vector<double>& coeficientes_F_O, vector<double>& recursos_restricciones, int n_variables_decision_heredado, int acumulador, int acumulador_inverso,vector<string>col,vector<string>fil) {
    int n_restricciones = tableau.size();  // Numero de restricciones M
    int n_variables_decision = tableau[0].size();  // Numero de variables de desicion N
    int cont=0;
    int c=2;
    int aumentador=0;
    for(int i=0;i<col.size();i++){
            if(i==0){col[i]="VB";}
            if(i==1){col[i]="Z";}
            if(i==col.size()-1){col[i]="LD";}
            if(i>1 && i<col.size()-3){
                if(cont <n_variables_decision_heredado){
                    col[i]="x"+to_string(cont+1);
                    cont=cont+1;
                    c=i+2;
                }else if(aumentador<acumulador){
                    col[c-1]="e"+to_string(aumentador+1);
                    col[c]="a"+to_string(aumentador+1);
                    c=c+2;
                    aumentador=aumentador+1;
                }else{
                    col[c-1]="h"+to_string(aumentador+1);
                    c=c+1;
                    aumentador=aumentador+1;
                }
            }
        }
        aumentador=0;
        //Llenado de vector para nombres de las filas
        for(int i=0;i<fil.size();i++){
            if(i==0){fil[i]="Z";}
            if(i>0){
                if(aumentador<acumulador){
                    fil[i]="a"+to_string(aumentador+1);
                    aumentador=aumentador+1;
                }else{
                    fil[i]="h"+to_string(aumentador+1);
                    aumentador=aumentador+1;
                }
            }
        }
    //Diferencia principal, para rellenar la matriz correctamente se le suma 2 veces i para que coloque correctamente los 1 y -1
        cout << "Metodo M2F" << endl;
        cout << "Tabla inicial M2F Fase 1:" << endl;
        printMatrix2(tableau,col,fil);

        tableau = GaussJordan(tableau,1,0,true, acumulador);

        cout << "Tabla inicial procesada con Gauss:" << endl;
        printMatrix2(tableau,col,fil);

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
            cout << "Menor Coeficiente: " << coeficiente_pivot << endl;
            cout << "/////////////////////////////////////////////////////" << endl;
            cout << endl;
            cout << "Tabla actual:" << endl;
            fil[pos_fila_pivot]=col[pos_variable_decision+1];
            printMatrix2(tableau,col,fil);
        }
        simplexaumentadotable(tableau, coeficientes_F_O, recursos_restricciones, n_variables_decision_heredado, acumulador, acumulador_inverso, false,col,fil);
        return tableau;
        }
    
vector<string>columna(vector<double> coeficientes_F_O,int acumulador,int acumulador_inverso){
    vector<string>col(coeficientes_F_O.size()+acumulador_inverso+(acumulador*2)+3);
    return col;
}
vector<string>fila(int acumulador,int acumulador_inverso){
    vector<string>fil(acumulador+acumulador_inverso+1);
    return fil;
}
void verificador(vector<double> coeficientes_F_O, vector<vector<double>> coeficientes_restricciones, vector<double> recursos_restricciones, vector<string> signos, bool is_minimize){
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
    cout << "Hay " << acumulador << " Restricciones de tipo " << signos[0] << " y " << acumulador_inverso << " del otro tipo" << endl;
    if (acumulador == signos.size()){
        mixto = false;
        if (signos[0] == "<="){
            cout << "-Programa procede a ejecutarse en modo aumentado" << endl;
            simplexaumentado(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_minimize, 0, 0);
        } else {
            cout << "-Programa procede a ejecutarse en modo M2F" << endl;
            simplexartificial(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_minimize, mixto, 0, 0);
        }
    } else {
        //signo es de mayor e igual a menor e igual, sino arrojaria valores erroneos
        mixto = true;
        cout << "-Programa procede a ejecutarse en primera instancia en modo M2F" << endl;
        //Creacion de la tabla
        vector<vector<double>> artificialmixto2 = simplexartificial(coeficientes_F_O, coeficientes_restricciones, recursos_restricciones, is_minimize, mixto, acumulador, acumulador_inverso);
        cout << "-Programa procede a ejecutarse en segunda instancia en modo aumentado" << endl;
        vector<vector<double>> aumentadomixto2 = simplexaumentadotable(artificialmixto2, coeficientes_F_O, recursos_restricciones, coeficientes_F_O.size(), acumulador, acumulador_inverso, mixto,columna(coeficientes_F_O,acumulador,acumulador_inverso),fila(acumulador,acumulador_inverso));
        
        //Resolucion del metodo
        cout << "-Resolucion del problema aplicando " << endl;
        simplexartificialtable(aumentadomixto2, coeficientes_F_O, recursos_restricciones, coeficientes_F_O.size(), acumulador, acumulador_inverso,columna(coeficientes_F_O,acumulador,acumulador_inverso),fila(acumulador,acumulador_inverso));
       
    }
    
}

//PARSER y extras

vector<double> extraer_coeficientes(string restriccion, int num_variables) {
    vector<double> coeficientes(num_variables, 0.0);
    regex patron("(-?\\d*\\.?\\d*)\\*?x(\\d+)");
    for (sregex_iterator i = sregex_iterator(restriccion.begin(), restriccion.end(), patron);
         i != sregex_iterator(); ++i) {
        smatch match = *i;
        double coeficiente = 1.0;
        if (match[1].length() > 0) {
            coeficiente = stod(match[1].str());
        }
        int indice_variable = stoi(match[2].str()) - 1;
        if (indice_variable < num_variables) {
            coeficientes[indice_variable] = coeficiente;
        }
    }
    return coeficientes;
}

double extraer_lado_derecho(string restriccion) {
    regex patron("<=|>=|<|>|=");
    sregex_iterator iterador(restriccion.begin(), restriccion.end(), patron);
    string simbolo = iterador->str();
    int pos_simbolo = iterador->position();
    string lado_derecho = restriccion.substr(pos_simbolo + simbolo.length());
    return stod(lado_derecho);
}

string extraer_simbolo(string funcion){
    std::regex desigualdad("<=|>=|<|>|=");
    std::smatch match;
    if (std::regex_search(funcion, match, desigualdad)) {
        return match.str();
    }
    return "";
}

int main() {
/*  //agregado de funciones por usuario
    int minimizar;
    cout << "deseas maximizar (1) o minimizar(2): " << endl;
    cin >> minimizar;
    
    if(minimizar == 2){
        bool minimizar = false;
    }else{
        bool minimizar = true;
    }

    int num_variables;
    cout << "Ingrese la cantidad de variables en las restricciones: ";
    //cin >> num_variables;

    int num_restricciones;
    cout << "Ingrese el número de restricciones: ";
    //cin >> num_restricciones;

    vector<string> restricciones(num_restricciones);
    vector<string> restricciones_mayor_o_igual;
    vector<string> restricciones_menor_o_igual;
    for (int i = 0; i < num_restricciones; ++i) {
        cout << "Ingrese la restricción " << i + 1 << ": ";
        //cin >> restricciones[i];
        if (restricciones[i].find(">=") != string::npos) {
            restricciones_mayor_o_igual.push_back(restricciones[i]);
        }
        else if (restricciones[i].find("<=") != string::npos) {
            restricciones_menor_o_igual.push_back(restricciones[i]);
        }
    }

    // Ordenar restricciones
    vector<string> restricciones_ordenadas;
    for (string restriccion : restricciones_mayor_o_igual) {
        restricciones_ordenadas.push_back(restriccion);
    }
    for (string restriccion : restricciones_menor_o_igual) {
        restricciones_ordenadas.push_back(restriccion);
    }

    string funcion_objetivo;
    cout << "Ingrese la función objetivo: ";
    //cin >> funcion_objetivo;

    vector<double> coeficientes_objetivo = extraer_coeficientes(funcion_objetivo, num_variables);

    vector<vector<double>> coeficientes_restricciones;
    vector<double> lado_derecho;
    vector<string> simbolos_restriciones;
    for (string restriccion : restricciones_ordenadas) {
        vector<double> coeficientes = extraer_coeficientes(restriccion, num_variables);
        string simbolo = extraer_simbolo(restriccion);
        double derecho = extraer_lado_derecho(restriccion);

        lado_derecho.push_back(derecho);
        simbolos_restriciones.push_back(simbolo);
        coeficientes_restricciones.push_back(coeficientes); // Agregar coeficientes a la lista
    } */

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

    vector<double> coeficientes_F_O6 = {-0.12, -0.15};
    vector<double> recursos_restricciones6 = {300,36,90,600};
    vector<vector<double>> coeficientes_restricciones6 = {{60, 60}, {12,6}, {10, 30}, {2,2}};
    vector<string> signos6 = {">=", ">=", ">=","<="};

    vector<double> coeficientes_F_O7 = {8, 6};
    vector<double> recursos_restricciones7 = {6,13,1,6};
    vector<vector<double>> coeficientes_restricciones7 = {{4, 3}, {1,3}, {1, 0}, {0,1}};
    vector<string> signos7 = {">=", ">=", "<=","<="};
    //for (int e = 0; e < simbolos_restriciones.size(); e++){
    //    cout << simbolos_restriciones[e] << endl;
   //}

    //Simplex Artificial solo aplica con restricciones de mayor o igual
    //verificador(coeficientes_objetivo, coeficientes_restricciones, lado_derecho, simbolos_restriciones, true);
    //verificador(coeficientes_F_O1, coeficientes_restricciones1, recursos_restricciones1, signos1, true);
    //verificador(coeficientes_F_O2, coeficientes_restricciones2, recursos_restricciones2, signos2, true);
    //verificador(coeficientes_F_O5, coeficientes_restricciones5, recursos_restricciones5, signos5, true);
    //verificador(coeficientes_F_O3, coeficientes_restricciones3, recursos_restricciones3, signos3, false);
    //verificador(coeficientes_F_O4, coeficientes_restricciones4, recursos_restricciones4, signos4, true);
    //verificador(coeficientes_F_O6, coeficientes_restricciones6, recursos_restricciones6, signos6, true);
    verificador(coeficientes_F_O7, coeficientes_restricciones7, recursos_restricciones7, signos7, false);
    //z = 44, x1 = 1, x2 = 6
    system("pause");
    return 0;
}