/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/* The following functions were implemented by Miguel Tasende. Copyright (C) 2016, Antel (Uruguay):
 - bli_sgemm_opt_mxn
 - bli_dgemm_opt_mxn
 - comparador

 The headers of the functions (and the "bli_cgemm_opt_mxn", "bli_zgemm_opt_mxn" imlpementations) are
 part of the BLIS framework. See the Copyright below.
 */

/*

   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2014, The University of Texas at Austin

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name of The University of Texas at Austin nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "blis.h"
#include <stdio.h>
#include "common.h"

#include <stdlib.h>
#include <stdbool.h>
#include <e-hal.h>
#include <e-loader.h>

#include <time.h>

#define _BufOffset (0x01000000)

void comparador(dim_t k,double* alpha,double* a1,double* b1,double* beta,double* c11,inc_t rs_c,inc_t cs_c);

void bli_sgemm_opt_mxn(
                        dim_t              k,
                        float*    restrict alpha,
                        float*    restrict a1,
                        float*    restrict b1,
                        float*    restrict beta,
                        float*    restrict c11, inc_t rs_c, inc_t cs_c,
                        auxinfo_t*         data
                      )
{
    //printf("Soy un kernel cores %i  -- N %i    |--------|   (mr,nr,cs_a,rs_b) = (%i,%i,%i,%i)\n",CORES,N,mr,nr,cs_a,rs_b);
	dim_t        mBlis    = bli_smr;
	dim_t        nBlis    = bli_snr;

	inc_t        cs_aBlis  = bli_spackmr;

	inc_t        rs_bBlis  = bli_spacknr;

    mBlis = N_FILAS;
    nBlis = N;
    cs_aBlis = N_FILAS;
    rs_bBlis = N;


    //nsub es el lado del cuadrado de multiplicación por core
    int nsub = NSUB;
    //ksub es la cantidad de columnas/filas de A/B que se pasan al Epiphany entero, por vez
    int ksub = KSUB;//CORES*1024/nBlis; //Asumo N = max(M,N)

    printf("bli_sgemm_opt_mxn(k = %i, alpha = %e, a1 = %e, b1 = %e, beta = %e, c11 = %e, rs_c = %i, cs_c = %i, data = ?)\n",k,*alpha,*a1,*b1,*beta,*c11,rs_c,cs_c);
    printf("mBlis = %i, nBlis = %i, cs_aBlis = %i, rs_bBlis = %i\n",mBlis,nBlis,cs_aBlis,rs_bBlis);

    dim_t              i;

    //Descubro cual es la dimension principal
/*    int tamCEpi = 0;
    if(cs_c > rs_c){
        tamCEpi = cs_c * nBlis;
    }
    else{
        tamCEpi = rs_c * mBlis;
    }


    float              c11epi[tamCEpi];
	float             *c11epiPtr;
	c11epiPtr = &c11epi[0];

	float             *amiliPtr, *bmiliPtr;*/
	float*            atemp;
	float*            btemp;
	float             ctemp[mBlis*nBlis];
	float*            ctempPtr;
	int               j, p, l;


    int fila,col, sFila, sCol, sProf;


    clock_t inicio, fin;
    double tTotal;

    //Tiempos del Epiphany
    clock_t inicioEpi, finEpi;
    double tEpi;
    //--------------------
    //Tiempo de postprocesamiento
    clock_t inicioPP, finPP;
    double tPP;
    //--------------------
    //Tiempo de postprocesamiento
    clock_t inicioCarga, finCarga;
    double tCarga;
    //--------------------

    //inicio = clock();//Se puede poner acá para medir la inicialización

    //inicio = clock();
    //Inicialización del Epiphany
    e_platform_t platform;
    e_epiphany_t dev;

    e_init(NULL);
    //e_reset_system();
    //inicio = clock();
    e_get_platform_info(&platform);

    e_mem_t memory;


//if (e_alloc(&memory, _BufOffset, mBlis*nBlis*sizeof(float)))
    if (e_alloc(&memory, _BufOffset, 0x01000000))
	{
		printf( "\nERROR: No se pudo reservar el espacio de memoria compartida!\n\n");
		exit(1);
	}
	if (e_open(&dev,0,0,platform.rows,platform.cols))
	{
		printf( "\nERROR: No se pudo abrir el grupo del Epiphany!\n\n");
		exit(1);
	}
	//e_reset_group(&dev);
	//e_load_group("/home/parallella/linpack/Matmul15s/bin/device.srec",&dev,0,0,platform.rows,platform.cols,E_TRUE);


    //e_alloc(&memory, _BufOffset, N*N*sizeof(float));

    //e_open(&dev,0,0,platform.rows,platform.cols);

    //e_reset_group(&dev);
    //e_load_group("/home/parallella/linpack/Matmul7v2/bin/device.srec",&dev,0,0,platform.rows,platform.cols,E_TRUE);
    //---------------------------------------------------

    //Esto es una inicialización para tener las dos versiones concurrentes (debug)+++++++++
    /*
    for(p=0;p<tamCEpi;p++){
        c11epi[p] = *(c11+p); //Si beta != 0, entonces cuando sumo el c11epi queda bien... esa es la idea.
    }
    */
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    //printf("Voy a inicializar C\n");

//    inicio = clock();//Se puede mover para atras para medir la inicialización

    tEpi = 0;
    tPP = 0;
    tCarga = 0;

    inicio = clock();
    float *resTemp;
    //Preparo c11epi con 0
/*    resTemp = &c11epi[0];
    for(p=0;p<nBlis;p++){
        for(l=0;l<mBlis;l++){
            *(resTemp+(l*rs_c)) = 0;
        }
        resTemp += cs_c;
    }*/
    //printf("Ya inicialicé C\n");
    //Ahora si acumulo alfa*A*B el resultado va a quedar bien (esto no es lo más eficiente si beta=0, pero facilita la lógica)
    //----------------------------------------------------------


    //Cargo los vectores de entrada A y B para la próxima iteración
    inicioCarga = clock();


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

    //e_write(&memory,0,0,OFFSET_A,atemp,ksub*mBlis*sizeof(float));
    //e_write(&memory,0,0,OFFSET_B,btemp,ksub*nBlis*sizeof(float));
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

    finCarga = clock();
    tCarga += (double)(finCarga-inicioCarga)/CLOCKS_PER_SEC;
    //-----------------------------


    //printf("Pasé de beta");
    //-----------------------------
    for(sProf=0;sProf<ktemp/ksub;sProf++){
        inicioCarga = clock();

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
            //printf("Iter única\n");
        }
        else if(sProf == 0){
            com = 0; // Inicio
            selector = 0;
            //printf("Iter Inicio\n");
        }
        else if(sProf == (ktemp/ksub)-1){
            com = 2; // Fin
            //printf("Iter Fin\n");
        }
        else{
            com = 1; //Iteración normal
            //printf("Iter normal\n");
        }

        //selector = 0;

        for(fila=0;fila<platform.rows;fila++){
            for(col=0;col<platform.cols;col++){
                e_write(&dev,fila,col,COMANDO,&com,sizeof(com));
                e_write(&dev,fila,col,LUGAR_DE_ENTRADA,&selector,sizeof(selector));
            }
        }

        finCarga = clock();
        tCarga += (double)(finCarga-inicioCarga)/CLOCKS_PER_SEC;
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->


/*
    //Cargo los vectores de entrada A y B
    inicioCarga = clock();
    atemp = a1 + ksub*(sProf)*cs_aBlis;
    btemp = b1 + ksub*rs_bBlis*(sProf);
    for(l=0;l<ksub/CORES;l++){ //Se puede probar si hay mejora por enviar bloques más grandes por vez, pero puede perder generalidad en ksub
        for(fila=0;fila<platform.rows;fila++){
            for(col=0;col<platform.cols;col++){
                e_write(&dev,fila,col,VEC_COL_A+l*mBlis*sizeof(float),atemp,mBlis*sizeof(float));
                e_write(&dev,fila,col,VEC_FILA_B+l*nBlis*sizeof(float),btemp,nBlis*sizeof(float));
                atemp += cs_aBlis;
                btemp += rs_bBlis;
            }
        }
    }

    finCarga = clock();
    tCarga += (double)(finCarga-inicioCarga)/CLOCKS_PER_SEC;
    //-----------------------------
*/


        int signalLocal;
        //Levanto señal de inicio del Epiphany<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        signalLocal = 1;
        inicioEpi = clock();
        e_write(&dev,0,0,COMPARTIDA_DEV,&signalLocal,sizeof(signalLocal));
        //printf("Inicié al Epiphany...\n");
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->



        //Cargo los vectores de entrada A y B para la próxima iteración
        inicioCarga = clock();
        if(sProf+1 < ktemp/ksub){//No es la penúltima iteración
            //printf("Estoy en una iteración intermedia");

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

/*
            for(p=0;p<ksub;p++){
                e_write(&memory,0,0,dirA+(p*mBlis*sizeof(float)),(atemp+(p*cs_aBlis)),mBlis*sizeof(float));
                e_write(&memory,0,0,dirB+(p*nBlis*sizeof(float)),(btemp+(p*rs_bBlis)),nBlis*sizeof(float));
            }*/

            //e_write(&memory,0,0,dirA,atemp,ksub*mBlis*sizeof(float));
            //e_write(&memory,0,0,dirB,btemp,ksub*nBlis*sizeof(float));

            finCarga = clock();
            tCarga += (double)(finCarga-inicioCarga)/CLOCKS_PER_SEC;
        }

        //-----------------------------



        //Espero a que el Epiphany termine<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        signalLocal = 1;
        while (signalLocal != 0){
            e_read(&dev, 0, 0, COMPARTIDA_DEV, &signalLocal, sizeof(signalLocal));
        }
        //printf("El Epiphany terminó!!\n");
        finEpi = clock();
        tEpi += (double)(finEpi-inicioEpi)/CLOCKS_PER_SEC;
        //printf("Tiempo parcial de procesamiento = %f\n",(double)(finEpi-inicioEpi)/CLOCKS_PER_SEC);
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->

    }

    //Multiplico por alfa y sumo beta*c11 (Postproc)----------------------------------------------------------
    inicioPP = clock();

    ctempPtr = &ctemp[0];
    e_read(&memory,0,0,0x0,ctempPtr,mBlis*nBlis*sizeof(float));

    resTemp = c11;

/*
    for(p=0; p<nBlis*mBlis; p++) {//Por ahora no tiene en cuenta el stride de C
        *(resTemp+p) = (*alpha)*(*(ctempPtr+p)) + (*beta)*(*(resTemp+p));
    }
*/
    for(p=0; p<nBlis; p++) {
        for(l=0;l<mBlis;l++){
            *(resTemp+p*cs_c+l*rs_c) = (*alpha)*(*(ctempPtr+p*mBlis+l)) + (*beta)*(*(resTemp+p*cs_c+l*rs_c));
        }
    }

    finPP = clock();
    tPP += (double)(finPP-inicioPP)/CLOCKS_PER_SEC;
    //----------------------------------------------------------------------------------------------------------


    fin = clock();
    tTotal = (double)(fin-inicio)/CLOCKS_PER_SEC;

    printf("Tiempo de preprocesamiento ARM y carga = %f (%2.1f%%) -- %f GFLOPS/s\n",tCarga,tCarga*100/tTotal,N*N_FILAS*(2*K-1)/tCarga/1000000000);
    printf("Tiempo procesamiento Epiphany = %f (%2.1f%%) -- %f GFLOPS/s\n",tEpi,tEpi*100/tTotal,N*N_FILAS*(2*K-1)/tEpi/1000000000);
    printf("Tiempo postprocesamiento ARM = %f (%2.1f%%) -- %f GFLOPS/s\n",tPP,tPP*100/tTotal,N*N_FILAS*(2*K-1)/tPP/1000000000);
    printf("Tiempo Total Epi-ARM = %f -- %f GFLOPS/s\n\n",tTotal,N*N_FILAS*(2*K-1)/tTotal/1000000000);

    printf("(Nota: los GFLOPS/s calculados son suponiendo que ese tramo fuera el único que afecta el tiempo utilizado)\n\n");

	// Close connection to device
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

    //e_finalize();

	/* Just call the reference implementation. */
	/*
	BLIS_DGEMM_UKERNEL_REF( k,
	                   alpha,
	                   a1,
	                   b1,
	                   beta,
	                   c11, rs_c, cs_c,
	                   data );

    float *cPtr = c11;
    float errorMaximo = 0;
    float errorPromedio = 0;

    float errorMaximoR = 0;
    float errorPromedioR = 0;

    float eTemp, eTempRelativo;
    //Calculo diferencias
    for(j=0;j<nBlis;j++){
        cPtr = c11 + j*cs_c;
        c11epiPtr = c11epi + j*cs_c;
        for(i=0;i<mBlis;i++){
            //printf("i=%i  resArm = %f   resEpiphany = %f \n",i,resARM[i],resEpiphany[i]);
            eTemp = (fabsf(*(cPtr+(i*rs_c)) - *(c11epiPtr+(i*rs_c))));
            if(*(cPtr+(i*rs_c)) != 0)
                eTempRelativo = eTemp/fabsf(*(cPtr+(i*rs_c)));
            else //Discutible...
                eTempRelativo = 0;

            errorPromedio = (errorPromedio*(i+j*K)+eTemp)/((i+j*K)+1);
            if(errorMaximo< eTemp)
                errorMaximo = eTemp;

            errorPromedioR = (errorPromedioR*(i+j*K)+eTempRelativo)/((i+j*K)+1);
            if(errorMaximoR< eTempRelativo)
                errorMaximoR = eTempRelativo;
        }
    }
    printf("Error Promedio: %e    Error Máximo: %e \n",errorPromedio,errorMaximo);
    printf("Error Promedio Relativo: %e    Error Máximo Relativo: %e \n\n",errorPromedioR,errorMaximoR);
*/
}



void bli_dgemm_opt_mxn(
                        dim_t              k,
                        double*   restrict alpha,
                        double*   restrict a1,
                        double*   restrict b1,
                        double*   restrict beta,
                        double*   restrict c11, inc_t rs_c, inc_t cs_c,
                        auxinfo_t*         data
                      )
{

    //printf("Soy un kernel cores %i  -- N %i    |--------|   (mr,nr,cs_a,rs_b) = (%i,%i,%i,%i)\n",CORES,N,mr,nr,cs_a,rs_b);
	dim_t        mBlis    = bli_dmr;
	dim_t        nBlis    = bli_dnr;

	inc_t        cs_aBlis  = bli_dpackmr;

	inc_t        rs_bBlis  = bli_dpacknr;

    mBlis = N_FILAS;
    nBlis = N;
    cs_aBlis = N_FILAS;
    rs_bBlis = N;

    int nsub = NSUB;
    int ksub = KSUB;

    //printf("bli_sgemm_opt_mxn(k = %i, alpha = %e, a1 = %e, b1 = %e, beta = %e, c11 = %e, rs_c = %i, cs_c = %i, data = ?)\n",k,*alpha,*a1,*b1,*beta,*c11,rs_c,cs_c);
    //printf("mBlis = %i, nBlis = %i, cs_aBlis = %i, rs_bBlis = %i\n",mBlis,nBlis,cs_aBlis,rs_bBlis);

    dim_t              i;

    float             aFloat[mBlis*k];
    float             bFloat[k*nBlis];
	float*            atemp;
	float*            btemp;
	float             ctemp[mBlis*nBlis];
	float*            ctempPtr;
	int               j, p, l;


    //Hago la conversión a float
    for(p=0;p<k;p++){
        for(l=0;l<mBlis;l++){
            aFloat[p*mBlis+l] = *(a1+p*cs_aBlis+l);
        }
    }
    cs_aBlis = mBlis;

    for(p=0;p<k;p++){
        for(l=0;l<nBlis;l++){
            bFloat[p*nBlis+l] = *(b1+p*rs_bBlis+l);
        }
    }
    rs_bBlis = nBlis;
    //----------------------------------------------------------------

    int fila,col, sFila, sCol, sProf;

    //Inicialización del Epiphany
    e_platform_t platform;
    e_epiphany_t dev;
    e_init(NULL);
    e_get_platform_info(&platform);

    e_mem_t memory;

    if (e_alloc(&memory, _BufOffset, 0x00600000))
	{
		printf( "\nERROR: No se pudo reservar el espacio de memoria compartida!\n\n");
		exit(1);
	}
	if (e_open(&dev,0,0,platform.rows,platform.cols))
	{
		printf( "\nERROR: No se pudo abrir el grupo del Epiphany!\n\n");
		exit(1);
	}

    //Cargo los vectores de entrada A y B para la próxima iteración
    //inicioCarga = clock();

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

    atemp = &aFloat[0];
    btemp = &bFloat[0];


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
            //printf("Iter única\n");
        }
        else if(sProf == 0){
            com = 0; // Inicio
            selector = 0;
            //printf("Iter Inicio\n");
        }
        else if(sProf == (ktemp/ksub)-1){
            com = 2; // Fin
            //printf("Iter Fin\n");
        }
        else{
            com = 1; //Iteración normal
            //printf("Iter normal\n");
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
        //inicioEpi = clock();
        e_write(&dev,0,0,COMPARTIDA_DEV,&signalLocal,sizeof(signalLocal));
        //printf("Inicié al Epiphany...\n");
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->



        //Cargo los vectores de entrada A y B para la próxima iteración
        //inicioCarga = clock();
        if(sProf+1 < ktemp/ksub){//No es la penúltima iteración
            //printf("Estoy en una iteración intermedia");

            atemp = &aFloat[0] + ksub*(sProf+1)*cs_aBlis;
            btemp = &bFloat[0] + ksub*rs_bBlis*(sProf+1);
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
        //-----------------------------


        //Espero a que el Epiphany termine<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        signalLocal = 1;
        while (signalLocal != 0){
            e_read(&dev, 0, 0, COMPARTIDA_DEV, &signalLocal, sizeof(signalLocal));
        }
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->

    }

    //Multiplico por alfa y sumo beta*c11 (Postproc)----------------------------------------------------------

    ctempPtr = &ctemp[0];
    e_read(&memory,0,0,0x0,ctempPtr,mBlis*nBlis*sizeof(float));

    double *resTemp;
    resTemp = c11;

    for(p=0; p<nBlis; p++) {
        for(l=0;l<mBlis;l++){
            *(resTemp+p*cs_c+l*rs_c) = (*alpha)*(*(ctempPtr+p*mBlis+l)) + (*beta)*(*(resTemp+p*cs_c+l*rs_c));
        }
    }

    //----------------------------------------------------------------------------------------------------------

	// Close connection to device
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

	/* Just call the reference implementation. */

/*	BLIS_DGEMM_UKERNEL_REF( k,
	                   alpha,
	                   a1,
	                   b1,
	                   beta,
	                   c11, rs_c, cs_c,
	                   data );*/


}



void bli_cgemm_opt_mxn(
                        dim_t              k,
                        scomplex* restrict alpha,
                        scomplex* restrict a1,
                        scomplex* restrict b1,
                        scomplex* restrict beta,
                        scomplex* restrict c11, inc_t rs_c, inc_t cs_c,
                        auxinfo_t*         data
                      )
{
	/* Just call the reference implementation. */
	BLIS_CGEMM_UKERNEL_REF( k,
	                   alpha,
	                   a1,
	                   b1,
	                   beta,
	                   c11, rs_c, cs_c,
	                   data );
}



void bli_zgemm_opt_mxn(
                        dim_t              k,
                        dcomplex* restrict alpha,
                        dcomplex* restrict a1,
                        dcomplex* restrict b1,
                        dcomplex* restrict beta,
                        dcomplex* restrict c11, inc_t rs_c, inc_t cs_c,
                        auxinfo_t*         data
                      )
{
	/* Just call the reference implementation. */
	BLIS_ZGEMM_UKERNEL_REF( k,
	                   alpha,
	                   a1,
	                   b1,
	                   beta,
	                   c11, rs_c, cs_c,
	                   data );
}


void comparador(inc_t k,double* alpha,double* a1,double* b1,double* beta,double* c11,inc_t rs_c,inc_t cs_c){
    /*--------------------------o--------------------------o-------------------------------------------------------------------*/
    //Esto es para comparar
    double *abij;
    double bj;
    double ai;
    int l,i,j;

    double *ctemp;

    //Inicializo la matriz de salida
    for(j=0; j<K; j++){
        ctemp = c11 + j*cs_c;
        for(i=0; i<K; i++){
            *ctemp = (*beta)*(*ctemp);
            ctemp += 1;
        }
    }

    for (l=0;l<K;++l){

		for(j=0;j<K;++j)
		{
			bj = *(b1 + j);
			abij = c11 + j*cs_c;

			for ( i = 0; i < K; ++i )
			{
				ai = *(a1 + i);

				*abij += (*alpha)*ai*bj;


				abij += 1;
			}
		}

		a1 += K;
		b1 += K;

	}
    //Fin del comparador
    /*--------------------------o--------------------------o-------------------------------------------------------------------*/

}
