//
//  matrix_double.h
//
//
//  Created by Kevin Mussi on 25/02/16.
//
//

#ifndef matrix_double_h
#define matrix_double_h

/*-------------------------------- MATRICI.H ---------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

/*---------------------------------- MACRO -----------------------------------*/

#define MATRIX matrix_double
#define ERR_VALUE_DBL -9999.9999
#define PARITÀ(n) (((n) % 2 == 0) ? 1 : -1)
#define EQUAL_TO_ZERO(X) (fabs(X) < 0.000001)
#define FREE_MATRIX(X) {free(X->data); free(X);}
#define MIN(X, Y) ((X < Y) ? X : Y)
#define MAX(X, Y) ((X > Y) ? X : Y)

/*--------------------------------- TYPEDEF ----------------------------------*/

typedef struct matrix_double_struct {
    int rows;       //n. di righe
    int cols;       //n. di colonne
    double *data;   //puntatore ai dati dinamici
    size_t size;    //dimensione della matrice
} matrix_double;

typedef struct coordinate_struct {
    int row;    //coordinata y
    int col;    //coordinata x
} coordinate;

/*------------------------------ DICHIARAZIONI -------------------------------*/

MATRIX *alloca_matrice(int R, int C);
void ottieni_matrice(MATRIX *M);
int ottieni_matrice_da_file(FILE *fp, MATRIX *M);
void stampa_matrice(MATRIX *M);
MATRIX *copia_matrice(MATRIX *M);
int ricerca(MATRIX *M, short (*confronta)(double elem), coordinate **trovati);
MATRIX *matrice_identità(int n);
MATRIX *matrice_diagonale(int DIM, ...);
MATRIX *matrice_trasposta(MATRIX *M);
MATRIX *matrice_inversa(MATRIX *M);
MATRIX *prodotto_matrici(MATRIX *A, MATRIX *B);
MATRIX *risolvi_sistema(MATRIX *A, MATRIX *b);
void trasforma_matrice(MATRIX *M, double (*trasforma)(double val));
void somma_matrici(MATRIX *A, MATRIX *B);
MATRIX *sottomatrice(MATRIX *M, int n_y, int n_x, ...);
double traccia_matrice(MATRIX *M);
int matrice_simmetrica(MATRIX *M);
double determinante_matrice(MATRIX *M);
double determinante_matrice_laplace(MATRIX *M);
MATRIX *elimina_riga_colonna(MATRIX *M, int Y, int X);
int rango_matrice(MATRIX *M);
MATRIX *matrice_a_scala(MATRIX *M);
void scambia_righe(MATRIX *M, int R1, int R2);
int colonna_nulla(MATRIX *M, int COL);
int conta_righe_nulle(MATRIX *M);

/*--------------------------------- FUNZIONI ---------------------------------*/

MATRIX *alloca_matrice(int R, int C) {
    MATRIX *M;
    
    if(R < 1 || C < 1)
        return NULL;
    
    if((M = (MATRIX *) malloc(sizeof(MATRIX))) == NULL)
        return NULL;
    
    if((M->data = (double *) calloc(R*C, sizeof(double))) == NULL) {
        free(M);
        return NULL;
    }
    
    M->rows = R;
    M->cols = C;
    M->size = (size_t) (R * C) * sizeof(double);
    
    return M;
}

void ottieni_matrice(MATRIX *M) {
    register int i, j;
    
    if(M == NULL || M->data == NULL)
        return;
    
    for(i = 0; i < M->rows; i++) {
        for(j = 0; j < M->cols; j++) {
            printf("\nInsert the element (%d, %d): ", i, j);
            scanf("%lf", &(M->data[i*M->cols + j]));
        }
    }
}

int ottieni_matrice_da_file(FILE *fp, MATRIX *M) {
    register int i;
    
    if(fp == NULL || M == NULL || M->data == NULL)
        return -1;
    
    for(i = 0; i < (M->cols)*(M->rows) && !feof(fp); i++)
        fscanf(fp, "%lf;", M->data + i);
    
    return i;
} 

void stampa_matrice(MATRIX *M) {
    register int i, j;
    
    if(M == NULL || M->data == NULL)
        return;
    
    for(i = 0; i < M->rows; i++) {
        printf("\n");
        for(j = 0; j < M->cols; j++)
            printf("%12lf", M->data[i*M->cols + j]);
    }
}

MATRIX *copia_matrice(MATRIX *M) {
    register int i, j;
    MATRIX *D;
    
    if(M == NULL || M->data == NULL)
        return NULL;
    
    if((D = alloca_matrice(M->rows, M->cols)) == NULL)
        return NULL;
    
    for(i = 0; i < M->rows; i++)
        for(j = 0; j < M->cols; j++)
            D->data[i*M->cols+j] = M->data[i*M->cols+j];
    
    return D;
}

int ricerca(MATRIX *M, short (*confronta)(double elem), coordinate **trovati) {
    register int i, j;
    coordinate temp[M->cols*M->rows];
    int cont = 0;
    
    if(M == NULL) {
        return -1;
    }
    
    for(i = 0; i < M->rows; i++) {
        for(j = 0; j < M->cols; j++) {
            if(confronta(M->data[i*M->cols+j])) {
                temp[cont].row = i;
                temp[cont].col = j;
                cont++;
            }
        }
    }
    
    if(cont != 0) {
        if((*trovati = (coordinate *) calloc(cont, sizeof(coordinate))) == NULL)
            return -1;
        
        for(i = 0; i < cont; i++)
            (*trovati)[i] = temp[i];
    }
    
    return cont;
}

MATRIX *matrice_identità(int n) {
    register int i;
    MATRIX *Id;
    
    if(n < 1)
        return NULL;
    
    if((Id = alloca_matrice(n, n)) == NULL)
        return NULL;
    
    for(i = 0; i < n; i++)
        Id->data[i*n + i] = 1.0;
        
    return Id;
}

MATRIX *matrice_diagonale(int DIM, ...) {
    register int i;
    MATRIX *D;
    va_list ap;
    
    if(DIM < 1)
        return NULL;
    
    if((D = alloca_matrice(DIM, DIM)) == NULL)
        return NULL;
    
    va_start(ap, DIM);
    
    for(i = 0; i < DIM; i++)
        D->data[i*DIM + i] = va_arg(ap, double);
    
    va_end(ap);
    
    return D;
}

MATRIX *matrice_trasposta(MATRIX *M) {
    register int i, j;
    MATRIX *T;
    
    if(M == NULL || M->data == NULL)
        return NULL;
    
    if((T = alloca_matrice(M->cols, M->rows)) == NULL)
        return NULL;
    
    for(i = 0; i < M->rows; i++)
        for(j = 0; j < M->cols; j++)
            T->data[j*T->cols+i] = M->data[i*M->cols+j];
    
    return T;
}

MATRIX *matrice_inversa(MATRIX *M) {
    register int i, j;
    MATRIX *T, *temp;
    double det;
    
    if(M == NULL || M->data == NULL || M->rows != M->cols)
        return NULL;
    
    if(EQUAL_TO_ZERO(det = determinante_matrice(M)))
        return NULL;
    
    if((T = alloca_matrice(M->rows, M->rows)) == NULL)
        return NULL;
    
    for(i = 0; i < M->rows; i++) {
        for(j = 0; j < M->rows; j++) {
            temp = elimina_riga_colonna(M, i, j);
            T->data[j*T->cols+i] = PARITÀ(i+j)*determinante_matrice(temp)/det;
            FREE_MATRIX(temp);
        }
    }
    
    return T;
}

MATRIX *prodotto_matrici(MATRIX *A, MATRIX *B) {
    register int i, j, k;
    MATRIX *M;
    
    if(A == NULL || B == NULL ||
       A->data == NULL || B->data == NULL ||
       A->cols != B->rows)
        return NULL;
    
    if((M = alloca_matrice(A->rows, B->cols)) == NULL)
        return NULL;
    
    for(i = 0; i < A->rows; i++)
        for(j = 0; j < B->cols; j++)
            for(k = 0; k < A->cols; k++)
                M->data[i*B->cols + j] += A->data[i*A->cols + k] *
                B->data[k*B->cols + j];
    
    return M;
}

MATRIX *risolvi_sistema(MATRIX *A, MATRIX *b) {
    MATRIX *X, *temp;
    
    if(A == NULL || b == NULL ||
       A->data == NULL || b->data == NULL ||
       b->cols != 1 || A->rows != b->rows)
        return NULL;
    
    X = prodotto_matrici((temp = matrice_inversa(A)), b);
    
    FREE_MATRIX(temp);
    
    return X;
}

void trasforma_matrice(MATRIX *M, double (*trasforma)(double val)) {
    register int i, j;
    
    if(M == NULL || M->data == NULL)
        return;
    
    for(i = 0; i < M->rows; i++)
        for(j = 0; j < M->cols; j++)
            M->data[i*M->cols + j] = trasforma(M->data[i*M->cols + j]);
}

void somma_matrici(MATRIX *A, MATRIX *B) {
    register int i, j;
    
    if(A == NULL || B == NULL ||
       A->data == NULL || B->data == NULL ||
       A->rows != B->rows || A->cols != B->cols)
        return;
    
    for(i = 0; i < A->rows; i++)
        for(j = 0; j < A->cols; j++)
            A->data[i*A->cols+j] += B->data[i*A->cols+j];
}

MATRIX *sottomatrice(MATRIX *M, int n_y, int n_x, ...) {
    /*
     * Il secondo parametro è il numero di righe della sottomatrice, il terzo
     * è il numero di colonne, gli n_y parametri successivi sono i numeri di
     * riga da prendere da M e gli n_x parametri successivi sono i numeri di
     * colonna da prendere da M. Questi numeri di riga e colonna sono intesi
     * come nell'algebra (a partire da 1) e non come nel C (a partire da 0).
     */
    register int i, j;
    MATRIX *S;
    int Y[n_y], X[n_x], temp;
    va_list ap;
    
    if(M == NULL || M->data == NULL ||
       n_y >= M->rows || n_y < 1 ||
       n_x >= M->cols || n_x < 1)
        return NULL;
    
    if((S = alloca_matrice(n_y, n_x)) == NULL)
        return NULL;
    
    temp = n_x;
    n_x += n_y;
    va_start(ap, n_x);
    
    for(i = 0; i < n_y; i++)
        Y[i] = va_arg(ap, int);
    for(i = 0; i < temp; i++)
        X[i] = va_arg(ap, int);
    
    va_end(ap);
    
    for(i = 0; i < temp; i++)
        for(j = 0; j < n_y; j++)
            S->data[j*temp + i] = M->data[(Y[j]-1)*M->cols + (X[i]-1)];
    
    return S;
}

double traccia_matrice(MATRIX *M) {
    register int i;
    double traccia = 0.0;
    
    if(M == NULL || M->data == NULL)
        return ERR_VALUE_DBL;
    
    for(i = 0; i < M->cols; i++)
        traccia += M->data[i*M->cols + i];
    
    return traccia;
}

int matrice_simmetrica(MATRIX *M) {
    register int i, j;
    
    if(M == NULL || M->data == NULL || M->rows != M->cols)
        return 0;
    
    for(i = 0; i < M->cols; i++)
        for(j = 0; j < M->cols && j != i; j++)
            if(M->data[i*M->cols+j] != M->data[j*M->cols+i])
                return 0;
    
    return 1;
}

double determinante_matrice(MATRIX *M) {
    /*
     * Calcolo del determinante utilizzando l'algoritmo di Gauss per ridurre
     * la matrice a scala e poi moltiplicando gli elementi sulla diagonale
     * principale. Il tempo di calcolo non aumenta esponenzialmente (come accade
     * con l'algoritmo di Laplace, anzi rimane pressoché costante (anche per
     * matrici di dimensione 10x10 o più).
     */
    register int i;
    double det = 1.0;
    MATRIX *S;
    
    if(M == NULL || M->data == NULL || M->rows != M->cols)
        return ERR_VALUE_DBL;
    
    if(M->cols == 1)
        return M->data[0];
    else if(M->cols == 2)
        return (M->data[0]*M->data[3] - M->data[1]*M->data[2]);
    
    if((S = matrice_a_scala(M)) == NULL)
        return ERR_VALUE_DBL;
    
    for(i = 0; i < S->cols; i++)
        det *= S->data[i*S->cols+i];
    
    FREE_MATRIX(S);
    return det;
}

double determinante_matrice_laplace(MATRIX *M) {
    /*
     * Calcolo del determinante utilizzando il teorema di Laplace partendo dalla
     * prima riga. Il tempo di calcolo aumenta come un fattoriale, rendendo
     * l'algoritmo adatto solo per matrici di dimensioni ridotte (massimo 5x5).
     */
    register int i;
    double det = 0.0;
    MATRIX *S;
    
    if(M == NULL || M->data == NULL || M->rows != M->cols)
        return ERR_VALUE_DBL;
    
    if(M->cols == 1) {
        det = M->data[0];
    } else {
        for(i = 0; i < M->cols; i++) {
            S = elimina_riga_colonna(M, 0, i);
            det += PARITÀ(i) * M->data[i] * determinante_matrice_laplace(S);
            FREE_MATRIX(S);
        }
    }
    
    return det;
}

MATRIX *elimina_riga_colonna(MATRIX *M, int Y, int X) {
    register int i, j;
    short k = 0, l = 0;
    MATRIX *S;
    
    if(M == NULL || M->data == NULL ||
       Y > M->rows-1 || Y < 0 ||
       X > M->cols-1 || X < 0)
        return NULL;
    
    if((S = alloca_matrice(M->rows-1, M->cols-1)) == NULL)
        return NULL;
    
    for(i = 0; i < M->rows; i++) {
        if(i == Y) {
            k = 1;
        } else {
            for(j = 0; j < M->cols; j++) {
                if(j == X)
                    l = 1;
                else
                    S->data[(i-k)*(M->cols-1)+(j-l)] = M->data[i*M->cols+j];
            }
            l = 0;
        }
    }
    
    return S;
}

int rango_matrice(MATRIX *M) {
    int rank;
    MATRIX *S;
    
    if(M == NULL || M->data == NULL)
        return -1;
    
    if((S = matrice_a_scala(M)) == NULL)
        return -2;
    
    rank = S->rows - conta_righe_nulle(S);
    FREE_MATRIX(S);
    return rank;
}

MATRIX *matrice_a_scala(MATRIX *M) {
    register int i, j, k;
    MATRIX *S;
    double pivot, num;
    int cont = 0;
    
    if(M == NULL || M->data == NULL)
        return NULL;
    
    if((S = copia_matrice(M)) == NULL)
        return NULL;
    
    //ALGORITMO DI GAUSS
    for(k = 0; k < MIN(S->rows, S->cols); k++) {
        if(colonna_nulla(S, k))
            continue;
        
        if(EQUAL_TO_ZERO(S->data[k*S->cols+k])) {
            for(i = k+1; i < S->rows; i++) {
                if(!EQUAL_TO_ZERO(S->data[i*S->cols+k])) {
                    scambia_righe(S, k, i);
                    cont++; //Conto il numero di scambi di righe
                    break;
                }
            }
        }
        
        pivot = S->data[k*S->cols+k];
        
        for(i = k+1; i < S->rows; i++) {
            if(!EQUAL_TO_ZERO(S->data[i*S->rows+k])) {
                num = (S->data[i*S->cols+k]);
                
                for(j = k; j < S->cols; j++) {
                    S->data[i*S->cols+j] -= (S->data[k*S->cols+j]) * num / pivot;
                }
            }
        }
    }
    
    /*
     * Se il numero di scambi di righe è dispari, cambio il segno dell'ultima
     * riga (arbitrariamente l'ultima). Questo perché ad ogni scambio il
     * determinante della matrice cambia, cambiando il segno di una riga la
     * matrice S ha lo stesso determinante di M.
     */
    if(cont % 2 == 1) {
        for(j = 0; j < S->cols; j++)
            if(!EQUAL_TO_ZERO(S->data[(S->rows-1)*S->cols+j]))
                S->data[(S->rows-1)*S->cols+j] *= -1;
    }
    
    return S;
}

void scambia_righe(MATRIX *M, int R1, int R2) {
    register int i;
    double temp;
    
    if(M == NULL || M->data == NULL ||
       R1 < 0 || R1 >= M->rows ||
       R2 < 0 || R2 >= M->rows)
        return;
    
    for(i = 0; i < M->cols; i++) {
        temp = M->data[R1*M->cols + i];
        M->data[R1*M->cols + i] = M->data[R2*M->cols + i];
        M->data[R2*M->cols + i] = temp;
    }
}

int colonna_nulla(MATRIX *M, int COL) {
    register int i;
    
    if(M == NULL || M->data == NULL || COL < 0 || COL >= M->rows)
        return -1;
    
    for(i = 0; i < M->rows; i++)
        if(!EQUAL_TO_ZERO(M->data[i*M->cols + COL]))
            return 0; //Caso di riga non nulla
    
    return 1; //Caso di riga nulla
}

int conta_righe_nulle(MATRIX *M) {
    register int i, j, cont = 0;
    
    for(i = 0; i < M->rows; i++) {
        for(j = 0; j < M->cols; j++) {
            if(!EQUAL_TO_ZERO(M->data[i*M->cols+j])) /* Caso di riga non nulla:
                                                      si passa alla successiva */
                j = M->cols;
            else if(j == M->cols - 1) /* Caso in cui tutta la riga sia nulla:
                                       il contatore si incrementa */
                cont++;
        }
    }
    
    return cont;
}

/*----------------------------------- FINE -----------------------------------*/

#endif /* matrix_double_h */
