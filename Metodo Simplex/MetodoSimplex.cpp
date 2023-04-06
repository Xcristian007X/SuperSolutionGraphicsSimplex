#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv)
{
    int filas;
    int columnas;

    cout << "Hola, Bienvenido a la Alpha del solucionario del Metodo Simplex" << endl;
    cout << "Por favor introduzca la dimension de la matriz" << endl;
    cout << "Primero empezemos con las filas:" << endl;
    cin >> filas;
    cout << "Excelente!, ahora con las columnas:" << endl;
    cin >> columnas;
    double matriz[filas][columnas];
    int filasA= (sizeof(matriz)/sizeof(matriz[0]));
    int columnasA = (sizeof(matriz[0])/sizeof(matriz[0][0]));
    cout << "Ahora rellenemos la matriz" << endl;
    for (int x=0; x <filas ;x++){
        for (int y=0; y < columnas ;y++){
            cout << "Posicion [" << x << "][" << y << "]: " << endl;
            cin >> matriz[x][y];
        };
    };
    cout << "Esta es la matriz generada: ";
    for (int x1=0; x1 <filas ;x1++){
        cout << endl;
        for (int y1=0; y1 < columnas ;y1++){
            cout << matriz[x1][y1] << " ";
        };
    };
    cout << endl;

    int zx[filas];
    int zy[columnas];
    int nx = 0;
    int ne = 0;
    int na = 0;
    int nh = 0;
    cout << "Ahora veamos las variables" << endl;
    cout << "cuantas variables x hay?:" << endl;
    cin >> nx;
    cout << "cuantas variables e hay?:" << endl;
    cin >> ne;
    cout << "cuantas variables artificiales hay?:" << endl;
    cin >> na;
    cout << "cuantas variables de holgura hay?:" << endl;
    cin >> nh;
    cout << nh << endl;
    if (na )
    if (na == 0)
    {
        cout << "jaja1" << endl;
        for (int x2=0; x2 < filasA; x2++){
            if (matriz[0][x2] < 0){
                cout << "jaja2" << endl;
            }
        }
    }
    cout << na << endl;
        
    cout << "Listo, espero que hayas podido hacer tu tarea, que tengas un buen dia :3" << endl;
    system("pause");
    return 0;
}

double proceso(int a, int b) {
    double Matriz[a][b];
    return Matriz[a][b];
}

