/* This code was implemented by Miguel Tasende. Copyright (C) 2016, Antel (Uruguay).*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <blis.h>

#include "common.h"
#include <stdbool.h>
#include <e-hal.h>
#include <e-loader.h>

#include "shmbuf.h"

#include "bli_gemm_opt_mxn_debug.h"

#define _BufOffset (0x01000000)

//void bli_sgemm_opt_mxn_noInit(int k,float* restrict alpha,float* restrict a1,float* restrict b1,float* restrict beta,float* restrict c11,int rs_c,int cs_c);
void imprimirMatriz(float *matriz, int filas, int columnas, int orden); //Para debug
#define POR_FILAS    0
#define POR_COLUMNAS 1

int main(int argc, char* const argv[]){
    int fd;
    struct shmbuf *shmp;

    //-------Consigo los punteros a memoria compartida e inicializo los semáforos--------------------
    shm_unlink("/matmul-file");
    fd = shm_open("/matmul-file",O_RDWR|O_CREAT,S_IRUSR|S_IWUSR);
    ftruncate(fd,sizeof(struct shmbuf));
    //Crear el mapeo de memoria
    shmp = mmap(0,sizeof(struct shmbuf),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    close(fd);

    sem_init(&shmp->semIni,1,0);
    sem_init(&shmp->semFin,1,0);

    //Inicialización del Epiphany -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
    e_epiphany_t dev;
    e_platform_t platform;
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
    //-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-


    //Ahora voy al loop de espera
    printf("\nEsperando pedido de cálculo...\n");
    while(1){ //Habría que poner una condición de fin mejor
        //printf("\nEsperando pedido de cálculo...\n");
        sem_wait(&shmp->semIni);

/*        printf("a = \n");
        imprimirMatriz(&shmp->a1[0],N_FILAS,N,POR_COLUMNAS);
        printf("\n\nb = \n");
        imprimirMatriz(&shmp->b1[0],N_FILAS,N,POR_FILAS);

        printf("k=%i, alpha=%f, beta=%f\n",shmp->k,shmp->alpha,shmp->beta);*/
/*
        float *ctemp, *atemp, *btemp;
        ctemp = &shmp->c11[0];
        atemp = &shmp->a1[0];
        btemp = &shmp->b1[0];
        int i,j,k;
        for(j=0;j<N;j++){
            for(i=0;i<N_FILAS;i++){
                *ctemp = 0;
                for(k=0;k<shmp->k;k++){
                    *ctemp += (*(atemp + k*N_FILAS + i)) * (*(btemp + k*N + j));
                }
                ctemp+=1;
            }
        }*/

        bli_sgemm_opt_mxn_noInit(shmp->k,&shmp->alpha,&shmp->a1[0],&shmp->b1[0],&shmp->beta,&shmp->c11[0],shmp->rs_c,shmp->cs_c,&dev,&memory,&platform,shmp->mBlis,shmp->nBlis,shmp->cs_aBlis,shmp->rs_bBlis);
        //bli_sgemm_opt_mxn_noInit(shmp->k,&shmp->alpha,&shmp->a1[0],&shmp->b1[0],&shmp->beta,&shmp->c11[0],1,N_FILAS,&dev,&memory,&platform,N_FILAS,N,N_FILAS,N);

/*        printf("\n\nResultado del Epiphany:  --------------------------------------------------------------------------------------------------------\n");
        imprimirMatriz(&shmp->c11[0],N_FILAS,N,POR_COLUMNAS);*/

        sem_post(&shmp->semFin);
    }


    //Finalización del Epiphany -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
    if (e_close(&dev)){
		printf( "\nERROR: No se pudo cerrar la conexión al Epiphany!\n\n");
		exit(1);
	}
	if (e_free(&memory)){
		printf( "\nERROR: No se pudo liberar la memoria compartida!\n\n");
		exit(1);
	}
    e_finalize();
    //-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-

    //Libero memoria
    munmap(shmp,sizeof(struct shmbuf)); //Innecesario porque se hace automáticamente cuando el programa termina

    return 0;
}



/*
void bli_sgemm_opt_mxn_noInit(int k,float* restrict alpha,float* restrict a1,float* restrict b1,float* restrict beta,float* restrict c11,int rs_c,int cs_c){
	//dim_t        mBlis    = bli_smr;
	//dim_t        nBlis    = bli_snr;

	//inc_t        cs_aBlis  = bli_spackmr;

	//inc_t        rs_bBlis  = bli_spacknr;

    //Para debug
    int mBlis = N_FILAS;
    int nBlis = N;
    int cs_aBlis = N_FILAS;
    int rs_bBlis = N;


    //nsub es el lado del cuadrado de multiplicación por core
    int nsub = NSUB;
    //ksub es la cantidad de columnas/filas de A/B que se pasan al Epiphany entero, por vez
    int ksub = KSUB;

    //printf("bli_sgemm_opt_mxn(k = %i, alpha = %e, a1 = %e, b1 = %e, beta = %e, c11 = %e, rs_c = %i, cs_c = %i, data = ?)\n",k,*alpha,*a1,*b1,*beta,*c11,rs_c,cs_c);
    //printf("mBlis = %i, nBlis = %i, cs_aBlis = %i, rs_bBlis = %i\n",mBlis,nBlis,cs_aBlis,rs_bBlis);

    int              i;
	float*            atemp;
	float*            btemp;
	float             ctemp[mBlis*nBlis];
	float*            ctempPtr;
	int               j, p, l;


    int fila,col, sFila, sCol, sProf;

    float *resTemp;

    //Cargo los vectores de entrada A y B para la próxima iteración
    //Si k no es múltiplo de KSUB, hago padding
    float filaNula[mBlis];
    float colNula[nBlis];
    int krestantes,ktemp;
    krestantes = k;
    ktemp = k;
    if(k%KSUB != 0){
        ktemp = ((k/ksub)+1)*ksub;
        for(l=0;l<mBlis;l++){
            filaNula[l] = 0;
        }
        for(l=0;l<nBlis;l++){
            colNula[l] = 0;
        }
    }

    atemp = a1;
    btemp = b1;

    if(krestantes >= ksub){
        for(p=0;p<ksub;p++){
            e_write(&memory,0,0,OFFSET_A+(p*mBlis*sizeof(float)),(atemp+(p*cs_aBlis)),mBlis*sizeof(float));
            e_write(&memory,0,0,OFFSET_B+(p*nBlis*sizeof(float)),(btemp+(p*rs_bBlis)),nBlis*sizeof(float));
        }
        krestantes -= ksub;
    }
    else if(krestantes==0){
        //No hago nada
    }
    else{
        for(p=0;p<krestantes;p++){
            e_write(&memory,0,0,OFFSET_A+(p*mBlis*sizeof(float)),(atemp+(p*cs_aBlis)),mBlis*sizeof(float));
            e_write(&memory,0,0,OFFSET_B+(p*nBlis*sizeof(float)),(btemp+(p*rs_bBlis)),nBlis*sizeof(float));
        }
        for(p=krestantes;p<ksub;p++){ //Padding
            e_write(&memory,0,0,OFFSET_A+(p*mBlis*sizeof(float)),&filaNula[0],mBlis*sizeof(float));
            e_write(&memory,0,0,OFFSET_B+(p*nBlis*sizeof(float)),&colNula[0],nBlis*sizeof(float));
        }
        krestantes = 0;
    }
    //-----------------------------

    for(sProf=0;sProf<ktemp/ksub;sProf++){

        //Configuro el comando para el Epiphany<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        int selector;
        if(sProf%2 == 0)
            selector = 0;
        else
            selector = 1;

        int com;
        if(ktemp/ksub == 1){
            com = 3; //Iteración única
            selector = 0;
        }
        else if(sProf == 0){
            com = 0; // Inicio
            selector = 0;
        }
        else if(sProf == (ktemp/ksub)-1){
            com = 2; // Fin
        }
        else{
            com = 1; //Iteración normal
        }

        for(fila=0;fila<platform.rows;fila++){
            for(col=0;col<platform.cols;col++){
                e_write(&dev,fila,col,COMANDO,&com,sizeof(com));
                e_write(&dev,fila,col,LUGAR_DE_ENTRADA,&selector,sizeof(selector));
            }
        }
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->


        int signalLocal;
        //Levanto señal de inicio del Epiphany<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        signalLocal = 1;
        e_write(&dev,0,0,COMPARTIDA_DEV,&signalLocal,sizeof(signalLocal));
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->


        //Cargo los vectores de entrada A y B para la próxima iteración
        if(sProf+1 < ktemp/ksub){//No es la penúltima iteración

            atemp = a1 + ksub*(sProf+1)*cs_aBlis;
            btemp = b1 + ksub*rs_bBlis*(sProf+1);
            off_t dirA, dirB;

            if(sProf%2 == 0){
                dirA = OFFSET_A_2;
                dirB = OFFSET_B_2;
            }
            else{
                dirA = OFFSET_A;
                dirB = OFFSET_B;
            }

            if(krestantes >= ksub){
                for(p=0;p<ksub;p++){
                    e_write(&memory,0,0,dirA+(p*mBlis*sizeof(float)),(atemp+(p*cs_aBlis)),mBlis*sizeof(float));
                    e_write(&memory,0,0,dirB+(p*nBlis*sizeof(float)),(btemp+(p*rs_bBlis)),nBlis*sizeof(float));
                }
                krestantes -= ksub;
            }
            else if(krestantes==0){
                //No hago nada
            }
            else{
                for(p=0;p<krestantes;p++){
                    e_write(&memory,0,0,dirA+(p*mBlis*sizeof(float)),(atemp+(p*cs_aBlis)),mBlis*sizeof(float));
                    e_write(&memory,0,0,dirB+(p*nBlis*sizeof(float)),(btemp+(p*rs_bBlis)),nBlis*sizeof(float));
                }
                for(p=krestantes;p<ksub;p++){ //Padding
                    e_write(&memory,0,0,dirA+(p*mBlis*sizeof(float)),&filaNula[0],mBlis*sizeof(float));
                    e_write(&memory,0,0,dirB+(p*nBlis*sizeof(float)),&colNula[0],nBlis*sizeof(float));
                }
                krestantes = 0;
            }
        }

        //Espero a que el Epiphany termine<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        signalLocal = 1;
        while (signalLocal != 0){
            e_read(&dev, 0, 0, COMPARTIDA_DEV, &signalLocal, sizeof(signalLocal));
        }
    }

    //Multiplico por alfa y sumo beta*c11 (Postproc)----------------------------------------------------------
    ctempPtr = &ctemp[0];
    e_read(&memory,0,0,0x0,ctempPtr,mBlis*nBlis*sizeof(float));

    resTemp = c11;

    for(p=0; p<nBlis; p++) {
        for(l=0;l<mBlis;l++){
            *(resTemp+p*cs_c+l*rs_c) = (*alpha)*(*(ctempPtr+p*mBlis+l)) + (*beta)*(*(resTemp+p*cs_c+l*rs_c));
        }
    }
}*/

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


