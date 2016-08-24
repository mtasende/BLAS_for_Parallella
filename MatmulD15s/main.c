/* This code was implemented by Miguel Tasende. Copyright (C) 2016, Antel (Uruguay).*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdlib.h>
#include <stdio.h>
#include <e-hal.h>
#include <e-loader.h>
#include <math.h>
#include "common.h"
#include "bli_gemm_opt_mxn_debug.h"
#include <time.h>

#include "main.h"

#define DEBUG_DEL_MATMUL 1

#define _BufOffset (0x01000000)

void imprimirMatriz(float *matriz, int filas, int columnas, int orden);
void fimprimirMatriz(FILE *f,float *matriz, int filas, int columnas, int orden);

int main(int argc, char *argv[]){

    //Defino las constantes como variables (es para que el programa sea más flexible, si quiero sacarlas de otro lado)
    int n = N;
    int m = N_FILAS;
    int k = K;
    int nsub = NSUB;
    //Lista la definición de constantes

    float a[m*k], b[k*n];
    float *a1 = a;
    float *b1 = b;
    //const int k = CORES;
    int i;
    int j;

    clock_t inicioARM, finARM;
    double tARM;

    //Inicializar vectores de entrada
	/*for (i=0; i<N*N; i++)
	{
		a[i] = i/2.0;
		//b es la Matriz Identidad x2
		if((i%N) == i/N)
            b[i] = 2.0;
        else
            b[i] = 0.0;
	}*/
	for (j=0;j<k;j++){
        for(i=0;i<m;i++){
            //a[j*m+i] = 1.0;
            //a[j*m+i] = 0.0;
            a[j*m+i] = (j*m+i)/89.0;
            //a[j*m+i] = -j-1;
            //a[j*m+i] = -i-1;
            //a[j*m+i] = j+1;
            //a[j*m+i] = i+1;
/*
            if(j==i)
                a[j*m+i] = 1;
            else
                a[j*m+i] = 0;
*/
        }
	}
	for(i=0;i<k;i++){
        for (j=0;j<n;j++){
            //b[i*n+j] = 0.0;
            b[i*n+j] = (i+j)/34.0;
            //b[i*n+j] = i+1;
            //b[i*n+j] = j+1;
            //b[i*n+j] = i*n+j+1; //Esto numera en orden normal
            //b[i*n+j] = 1.0;
/*
            if(j==i)
                b[i*n+j] = 1;
            else
                b[i*n+j] = 0;
*/

        /*
            //b es la Matriz Identidad x2 repetida cada CORES
            if((j%CORES) == i)
                b[i*N+j] = 2.0;
            else
                b[i*N+j] = 0.0;*/
            }
        }


    float alpha = 1.0;//15.7;
    float beta = 0.0;//300.5;
    float resEpiphany[m*n];
    float *resEpiphany_ = &resEpiphany[0];

   for (j=0;j<n;j++){
        for(i=0;i<m;i++){
            resEpiphany[j*m+i] = (j*m+i)/7;
        }
	}

/*
    for (j=0;j<n;j++){
        for(i=0;i<2*m;i++){
            resEpiphany[j*2*m+i] = 0;
        }
	}
*/
    /*--------------------------o--------------------------o-------------------------------------------------------------------*/
    //Esto es para comparar
    float resARM[m*n];
    float *abij;
    float bj;
    float ai;
    int l;

    inicioARM = clock();
    //Inicializo la matriz de salida
    for(i=0; i<m*n; i++){
		resARM[i] = 0;
	}

    for (l=0;l<k;++l){
		abij = resARM;

		for(j=0;j<n;++j)
		{
			bj = *(b1 + j);

			for ( i = 0; i < m; ++i )
			{
				ai = *(a1 + i);

				*abij += ai*bj;

				//printf("l= %i  (i,j) = (%i,%i) ...... ai = %.2f   bj = %.2f   abij = %.2f\n",l,i,j,ai,bj,*abij);

				abij += 1;
			}
		}

		a1 += m;
		b1 += n;

		//printf("Matriz parcial %i:\n",l);
		//imprimirMatriz(resARM,N,N,POR_COLUMNAS);
	}

    for(i=0; i<m*n; i++){
		resARM[i] = alpha*resARM[i] + beta*resEpiphany[i];
	}
    //Fin del comparador
    finARM = clock();
    tARM = (double)(finARM - inicioARM)/CLOCKS_PER_SEC;
    printf("Tiempo de procesamiento del ARM = %f\n",tARM);
    /*--------------------------o--------------------------o-------------------------------------------------------------------*/


    //Inicializo la matriz de salida
    /*for(i=0; i<K*K; i++){
		resEpiphany[i] = 0.0;
	}*/
	a1 = a;
    b1 = b;


    //Inicialización del Epiphany<-init-><-init-><-init-><-init-><-init-><-init-><-init-><-init-><-init-><-init-><-init-><-init-><-init->
    e_platform_t platform;
    e_epiphany_t dev;
    e_mem_t memory;

    e_init(NULL);
    e_reset_system();
    e_get_platform_info(&platform);
    if (e_alloc(&memory, _BufOffset, 0x00600000)){
        printf( "\nERROR: No se pudo reservar el espacio de memoria compartida!\n\n");
        exit(1);
    }
    if (e_open(&dev,0,0,platform.rows,platform.cols)){
        printf( "\nERROR: No se pudo abrir el grupo del Epiphany!\n\n");
        exit(1);
    }
	e_load_group("/home/parallella/linpack/MatmulD15s/bin/device.srec",&dev,0,0,platform.rows,platform.cols,E_TRUE);
	//<-/init-><-/init-><-/init-><-/init-><-/init-><-/init-><-/init-><-/init-><-/init-><-/init-><-/init-><-/init-><-/init->

    //bli_sgemm_opt_mxn_debug(K,&alpha,a1,b1,&beta,resEpiphany_,1,N);
    //auxinfo_t         data;
    //bli_dgemm_opt_mxn_debug(K,&alpha,a1,b1,&beta,resEpiphany_,1,N,&data);
    bli_sgemm_opt_mxn_noInit(k,&alpha,a1,b1,&beta,resEpiphany_,1,m,&dev,&memory,&platform,N_FILAS,N,N_FILAS,N);


	//Termino con el Epiphany<-finalize-><-finalize-><-finalize-><-finalize-><-finalize-><-finalize-><-finalize-><-finalize->
	if (e_close(&dev))
	{
		printf( "\nERROR: No se pudo cerrar la conexión al Epiphany!\n\n");
		exit(1);
	}
	if (e_free(&memory))
	{
		printf( "\nERROR: No se pudo liberar la memoria compartida!\n\n");
		exit(1);
	}
    e_finalize();
    //<-/finalize-><-/finalize-><-/finalize-><-/finalize-><-/finalize-><-/finalize-><-/finalize-><-/finalize-><-/finalize->


    //Muestro los resultados


/*
    float errorMaximo = fabsf(resARM[0] - resEpiphany[0]);
    float errorPromedio = fabsf(resARM[0] - resEpiphany[0]);

    float errorMaximoR = errorMaximo/fabsf(resARM[0]);
    float errorPromedioR = errorPromedio/fabsf(resARM[0]);
  */

    float errorMaximo = 0;
    float errorPromedio = 0;

    float errorMaximoR = 0;
    float errorPromedioR = 0;

    float eTemp, eTempRelativo;
    //Calculo diferencias
    for(i=0;i<m*n;i++){
        //printf("i=%i  resArm = %f   resEpiphany = %f \n",i,resARM[i],resEpiphany[i]);
        eTemp = (fabsf(resARM[i] - resEpiphany[i]));
        eTempRelativo = eTemp/fabsf(resARM[i]);

        //errorPromedio = (errorPromedio*i+eTemp)/(i+1);
        errorPromedio += eTemp;
        if(errorMaximo< eTemp)
            errorMaximo = eTemp;

        //errorPromedioR = (errorPromedioR*i+eTempRelativo)/(i+1);
        errorPromedioR += eTempRelativo;
        if(errorMaximoR< eTempRelativo)
            errorMaximoR = eTempRelativo;
    }

    errorPromedio = errorPromedio/(m*n);
    errorPromedioR = errorPromedioR/(m*n);

    printf("Error Promedio: %e    Error Máximo: %e \n",errorPromedio,errorMaximo);
    printf("Error Promedio Relativo: %e    Error Máximo Relativo: %e \n",errorPromedioR,errorMaximoR);

/*
    printf("a = \n");
    imprimirMatriz(&a[0],N_FILAS,N,POR_COLUMNAS);
    printf("\n\nb = \n");
    imprimirMatriz(&b[0],N_FILAS,N,POR_FILAS);

    printf("Resultado del ARM:\n");
    imprimirMatriz(&resARM[0],N_FILAS,N,POR_COLUMNAS);
    printf("Resultado del Epiphany:\n");
    imprimirMatriz(&resEpiphany[0],N_FILAS,N,POR_COLUMNAS);
*/
/*
    printf("Resultado del ARM:\n");
    imprimirMatriz(&resARM[0],N_FILAS,N,POR_COLUMNAS);
    printf("\n\nResultado del Epiphany:  --------------------------------------------------------------------------------------------------------\n");
    imprimirMatriz(&resEpiphany[0],N_FILAS,N,POR_COLUMNAS);
*/
/*
    printf("\n\nEPI 2 ----- EPI 2 ---- EPI 2 ----- EPI 2 ---- EPI 2 --------- EPI 2 ---- EPI 2 --------- EPI 2 ---- EPI 2 --------- EPI 2 ---- EPI 2 ----:\n");
    imprimirMatriz(&resEpiphany[N_FILAS*N],N_FILAS,N,POR_COLUMNAS);
*/

/*
    FILE *fEpi, *fArm ,*fA, *fB;
    fEpi = fopen("resEpi.txt","w");
    fArm = fopen("resArm.txt","w");
    fA = fopen("matA.txt","w");
    fB = fopen("matB.txt","w");

    fimprimirMatriz(fArm,&resARM[0],K,K,POR_COLUMNAS);
    fimprimirMatriz(fEpi,&resEpiphany[0],K,K,POR_COLUMNAS);
    fimprimirMatriz(fA,&a[0],K,K,POR_COLUMNAS);
    fimprimirMatriz(fB,&b[0],K,K,POR_FILAS);
*/

}

void imprimirMatriz(float *matriz, int filas, int columnas, int orden){
    int i;
    int j;


    for(i=0;i<filas;i++){
        for(j=0;j<columnas;j++){
            if(orden == POR_FILAS){
                printf("%.2f ",matriz[i*columnas+j]);
            }
            else if(orden == POR_COLUMNAS){
                printf("%.2f ",matriz[i+j*filas]);
            }
        }
        printf("\n");
    }
    printf("\n");
}

void fimprimirMatriz(FILE *f,float *matriz, int filas, int columnas, int orden){
    int i;
    int j;


    for(i=0;i<filas;i++){
        for(j=0;j<columnas;j++){
            if(orden == POR_FILAS){
                fprintf(f,"%.2f ",matriz[i*columnas+j]);
            }
            else if(orden == POR_COLUMNAS){
                fprintf(f,"%.2f ",matriz[i+j*filas]);
            }
        }
        fprintf(f,"\n");
    }
    fprintf(f,"\n");
}
