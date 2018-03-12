//
//  main.c
//  GAL
//
//  Created by Kevin Mussi on 29/11/15.
//  Copyright © 2015 Kevin Mussi. All rights reserved.
//

#include "util/matrix_double.h"
#include <gsl/gsl_poly.h>
#include <time.h>

#define CENTRO_2 \
{ \
temp = alloca_matrice(2, 1); \
c = alloca_matrice(2, 1); \
c->data[0] = -B->data[2]; \
c->data[1] = -B->data[5]; \
temp = risolvi_sistema(A, c); \
c->data[0] = temp->data[0]; \
c->data[1] = temp->data[1]; \
STAMPA2("\n    Center: (%lf, %lf).", c->data[0], c->data[1]); \
FREE_MATRIX(temp); \
FREE_MATRIX(c); \
} \

#define CENTRO_3 \
{ \
temp = alloca_matrice(3, 1); \
c = alloca_matrice(3, 1); \
c->data[0] = -B->data[3]; \
c->data[1] = -B->data[7]; \
c->data[2] = -B->data[11]; \
temp = risolvi_sistema(A, c); \
c->data[0] = temp->data[0]; \
c->data[1] = temp->data[1]; \
c->data[2] = temp->data[2]; \
STAMPA2("\n    Center: (%lf, %lf,", c->data[0], c->data[1]); \
STAMPA1(" %lf).", c->data[2]); \
FREE_MATRIX(temp); \
FREE_MATRIX(c); \
} \

#define STAMPA(X) \
{ \
printf(X); \
fprintf(fp, X); \
} \

#define STAMPA1(X, Y) \
{ \
printf(X, Y); \
fprintf(fp, X, Y); \
} \

#define STAMPA2(X, Y, Z) \
{ \
printf(X, Y, Z); \
fprintf(fp, X, Y, Z); \
} \

int coniche(MATRIX *B);
int quadriche(MATRIX *B);
void rand_ins(void);
void stampa_matrice_gal(MATRIX *M);

FILE *fp;

int main(int argc, const char * argv[]) {
    char scelta;
    MATRIX *M;
    int DIM;
    time_t date = time(NULL);
    clock_t start, end;
    
    if((fp = fopen("conics_quadrics_classifier_log.txt", "at")) == NULL)
        return -1;
    
    do {
        printf("\nConic or Quadric? [c / q]:\n>> ");
        scanf("%c", &scelta);
        if(scelta != 'c' && scelta != 'q') {
            printf(" Pleasy retry.");
        }
    } while (scelta != 'c' && scelta != 'q');
    
    (scelta == 'c') ? (DIM = 3) : (DIM = 4);
    
    if((M = alloca_matrice(DIM, DIM)) == NULL)
        return -2;
    
    ottieni_matrice(M);
    
    while(!matrice_simmetrica(M)) {
        printf("\nThe matrix you entered isn't symmetric! Please retry.\n");
        ottieni_matrice(M);
    }
    
    fprintf(fp, "%s\n%s\n", ctime(&date), (scelta == 'c') ? "Conic" : "Quadric");
    STAMPA("\nMatrix:");
    stampa_matrice_gal(M);
    
    start = clock();
    
    (scelta == 'c') ? coniche(M) : quadriche(M);
    
    end = clock();
    
    STAMPA1("\n\nComputation time: %.3lfms", 1000.0*(double)(end-start) /
                                                 (double)(CLOCKS_PER_SEC));
    fprintf(fp, "\n------------------------------------------------------------\n");
    printf("\n\n");
    FREE_MATRIX(M);
    fclose(fp);
    return 0;
}

void stampa_matrice_gal(MATRIX *M) {
    register int i, j;
    
    for(i = 0; i < M->cols; i++) {
        STAMPA("\n");
        for(j = 0; j < M->cols; j++) {
            STAMPA1("%12lf", M->data[i*M->cols + j]);
            if(j == M->cols-2) {
                STAMPA(" | ");
            }
        }
        if(i == M->cols-2) {
            STAMPA("\n");
            for(j = 0; j < M->cols; j++) {
                STAMPA("------------");
                if(j == M->cols-2) {
                    STAMPA("-+-");
                }
            }
        }
    }
}

int coniche(MATRIX *B) {
    double i1, i2, i3, e[2], E, v[2], G, H;
    MATRIX *A, *temp, *c;
    int rank_A, rank_B;
    
    A = elimina_riga_colonna(B, 2, 2);
    
    i1 = traccia_matrice(A);
    i2 = determinante_matrice(A);
    i3 = determinante_matrice(B);
    rank_A = (i2 == 0) ? rango_matrice(A) : 2;
    rank_B = (i3 == 0) ? rango_matrice(B) : 3;
    
    STAMPA1("\n\nInvariants:\n    I1 = %lf.", i1);
    STAMPA2("\n    I2 = %lf --> %s.", i2, (i2 == 0) ? "conic without center" : "conic with center");
    STAMPA2("\n    I3 = %lf --> %s.", i3, (i3 == 0) ? "degenerate conic" : "non-degenerate conic");
    STAMPA2("\n    Rank of A = %d.\n    Rank of B = %d.", rank_A, rank_B);
    STAMPA("\n\nConclusions:");
    
    gsl_poly_solve_quadratic(1.0, -i1, i2, &e[0], &e[1]);
    
    if(i3 != 0.0) {
        if(i2 != 0.0) {
            if(i2 > 0.0 && i1*i3 < 0.0) {
                if(A->data[0] == A->data[3]) {
                    STAMPA("\n    Circumference.");
                } else {
                    STAMPA("\n    Real ellipsis.");
                }
            } else if(i2 < 0.0) {
                if(i1 == 0) {
                    STAMPA("\n    Equilateral hyperbola.");
                } else {
                    STAMPA("\n    Hyperbola.");
                }
            } else {
                STAMPA("\n    Immaginary ellipsis.");
            }
            
            STAMPA2("\n    Canonic form: %lfx^2 + %lfy^2", e[0], e[1]);
            STAMPA1(" = %lf.", -i3/i2);
            
            CENTRO_2;
            
        } else {
            
            STAMPA("\n    Parabola.");
            
            if(EQUAL_TO_ZERO(e[0]))
                E = e[1];
            else
                E = e[0];
            
            STAMPA2("\n    Canonic form: %lfy^2 + %lfx = 0.", E, 2*sqrt(-i3/i1));
            
            STAMPA2("\n    Axis: %lfx + %lfy", B->data[0]*B->data[0]+
                                               B->data[1]*B->data[1],
                                               B->data[1]*(B->data[0]+B->data[4]));
            STAMPA1(" + %lf = 0.", B->data[0]*B->data[2]+B->data[1]*B->data[5]);
            
            //Calcolo dei vertici
            if(EQUAL_TO_ZERO(B->data[1])) {
                v[0] = -B->data[2]/B->data[0];
                v[1] = -B->data[5]/B->data[4];
            } else {
                G = B->data[1]*(B->data[0]+B->data[4]) /
                    (B->data[0]*B->data[0]+B->data[1]*B->data[1]);
                H = (B->data[0]*B->data[2]+B->data[1]*B->data[5]) /
                    (B->data[0]*B->data[0]+B->data[1]*B->data[1]);
                gsl_poly_solve_quadratic(B->data[0]*G*G-2*B->data[1]*G+B->data[4],
                                         2*(B->data[0]*G*H - B->data[1]*H -
                                            B->data[2]*G+B->data[5]),
                                         B->data[0]*H*H-2*B->data[2]*H+B->data[8],
                                         &v[1], &v[0]);
                v[0] = -G*v[1]-H;
            }
            STAMPA2("\n    Vertex: (%lf, %lf).", v[0], v[1]);
        }
    } else {
        if(rank_A == 2) {
            STAMPA("\n    Pair of intersecting lines.");
            STAMPA2("\n    Canonic form: %lfx^2 + %lfy^2 = 0.", e[0], e[1]);
            
            CENTRO_2;
            
        } else if(rank_A == 1 && rank_B == 2) {
            STAMPA("\n    Pair of parallel lines.");
            
            if(EQUAL_TO_ZERO(e[0]))
                E = e[1];
            else
                E = e[0];
            
            STAMPA2("\n    Canonic form: %lfy^2 + %lf = 0.", E, 2*sqrt(-i3/i1));
            
        } else if(rank_A == 1 & rank_B == 1) {
            STAMPA("\n    Double line.");
            STAMPA("\n    Canonic form: y^2 = 0");
        }
    }
    
    FREE_MATRIX(A);
    
    return 0;
}

int quadriche(MATRIX *B) {
    double i1, i2, i3, i4, e[3], E1, E2;
    MATRIX *A, *c, *temp;
    MATRIX *temp1, *temp2, *temp3;
    int rank_A, rank_B;
    
    A = elimina_riga_colonna(B, 3, 3);
    temp1 = elimina_riga_colonna(A, 2, 2);
    temp2 = elimina_riga_colonna(A, 1, 1);
    temp3 = elimina_riga_colonna(A, 0, 0);
    
    i1 = traccia_matrice(A);
    i2 = determinante_matrice(temp1) +
         determinante_matrice(temp2) +
         determinante_matrice(temp3);
    i3 = determinante_matrice(A);
    i4 = determinante_matrice(B);
    rank_A = (i3 == 0) ? rango_matrice(A) : 3;
    rank_B = (i4 == 0) ? rango_matrice(B) : 4;
    
    FREE_MATRIX(temp1);
    FREE_MATRIX(temp2);
    FREE_MATRIX(temp3);
    
    STAMPA2("\n\nInvariants:\n    I1 = %lf.\n    I2 = %lf.", i1, i2);
    STAMPA2("\n    I3 = %lf --> %s.", i3, (i3 == 0) ? "quadric without center" : "quadric with center");
    STAMPA2("\n    I4 = %lf --> %s.", i4, (i4 == 0) ? "degenerate quadric" : "non-degenerate quadric");
    STAMPA2("\n    Rank of A = %d.\n    Rank of B = %d.", rank_A, rank_B);
    STAMPA("\n\nConclusions:");
    
    gsl_poly_solve_cubic(-i1, i2, -i3, &e[0], &e[1], &e[2]);
    
    if(i4 != 0.0) {
        if(i3 != 0.0) {
            if(i4 > 0.0 && i2 > 0.0 && i1*i3 > 0.0) {
                STAMPA("\n    Immaginary ellipsoid.");
            } else if(i4 < 0.0 && i2 > 0.0 && i1*i3 > 0.0) {
                STAMPA("\n    Real ellipsoid.");
            } else if(i4 > 0.0 && (i2 <= 0.0 || i1*i3 <= 0.0)) {
                STAMPA("\n    One-sheeted hyperboloid.");
            } else if(i4 < 0.0 && (i2 <= 0.0 || i1*i3 <= 0.0)) {
                STAMPA("\n    Two-sheeted hyperboloid.");
            }
            
            STAMPA2("\n    Canonic form: %lfx^2 + %lfy^2", e[0], e[1]);
            STAMPA2(" + %lfz^2 = %lf.", e[2], -i4/i3);
            
            CENTRO_3;
            
        } else {
            if(i4 < 0.0) {
                STAMPA("\n    Elliptic paraboloid.");
            } else {
                STAMPA("\n    Hyperbolic paraboloid.");
            }
            
            if(EQUAL_TO_ZERO(e[0])) {
                E1 = e[1];
                E2 = e[2];
            } else {
                E1 = e[0];
                if(EQUAL_TO_ZERO(e[1]))
                    E2 = e[2];
                else
                    E2 = e[1];
            }
            
            STAMPA2("\n    Canonic form: %lfx^2 + %lfy^2", E1, E2);
            STAMPA1(" + %lfz=0.", 2*sqrt(-i4/i2));
        }
    } else {
        if(i3 != 0.0) {
            if(i3 > 0.0) {
                STAMPA("\n    Immaginary cone.");
            } else {
                STAMPA("\n    Real cone.");
            }
            
            CENTRO_3;
            
        } else {
            if(i2 > 0.0) {
                STAMPA("\n    Elliptic cylinder (real or immaginary).");
            } else if(i2 < 0.0) {
                STAMPA("\n    Hyperbolic cylinder.");
            } else if(i2 == 0.0 && rank_B == 3) {
                STAMPA("\n    Parabolic cylinder.");
            } else if(i2 == 0.0 && rank_B < 3) {
                STAMPA("\n    Pair of planes (real or immaginary) (parallel or intersecting).");
            }
        }
    }
    
    FREE_MATRIX(A);
    
    return 0;
}
