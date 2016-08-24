/* This code was implemented by Miguel Tasende. Copyright (C) 2016, Antel (Uruguay).*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <stdlib.h>
#include "e-lib.h"
#include "common.h"
#include <math.h>
#include <stdbool.h>

asm(".global __stack_start_;");
asm(".set __stack_start_,0x5FF0;");

void coreMatch(int fila,int col,int *filaM,int *colM) __attribute__ ((section (".text")));
void coreAnt(int fila,int col,int *filaM,int *colM) __attribute__ ((section (".text")));
int main(void) __attribute__ ((section (".text")));
void coords(int coreid,int *fila,int *col) __attribute__ ((section (".text")));
void subMatmul(float *aIni,float *bIni, float *resInterno, float *resExternoIni) __attribute__ ((section (".text")));
void copiarDatos(float *fuente, float *destino) __attribute__ ((section (".text")));

void matmul_assembly(float *a, float *b, float *c, float *cOut);


//Variables globales
//float resEx[N_FILAS*N*sizeof(float)*K/KSUB] __attribute__ ((section ("shared_dram")));
float resEx[N_FILAS*N*sizeof(float)] __attribute__ ((section ("shared_dram")));
e_coreid_t coreid;
int fila, col, filaM, colM;
int coreNum;
char *salida;
volatile int *estado;
volatile e_barrier_t barriers[16];
e_barrier_t *tgt_bars[16];
int *offsetResp, *offsetA, *offsetB;


int main(void){


        salida = (char *) MENSAJES;
        estado = (int *) COMPARTIDA_DEV;
        //offsetResp = (int *) OFFSET_RESP;
        //offsetA = (int *) OFFSET_A;
        //offsetB = (int *) OFFSET_B;

        coreid = e_get_coreid();
        coords(coreid,&fila,&col);
        coreNum = fila + LADO_EPI*col;
        coreMatch(fila,col,&filaM,&colM);


        *estado = 0;
    while(true){

        e_barrier_init(barriers,tgt_bars);

        //Espera de la señal de inicio<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        if (coreNum == 0){
            //Espero del ARM la señal de inicio
			while(*estado == 0) {};
		}
		//Todos los núcleos esperan al cero (que recibe del ARM)
		e_barrier(barriers, tgt_bars);
		//<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->




        //Ejecución-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
        float *a, *b;

        salida = (char *) MENSAJES;


        //-----------------------------------------------Acá empieza la gran función------------------------
/*
        float * res;
        res = &resEx[coreNum*NSUB*NSUB];
        int i;
        for(i=0;i<NSUB*NSUB;i++){
            *res = coreNum;
            res += 1;
        }
        */

        float *res1,*res2,*res1Ext,*res2Ext,*acumulador;
        res1 = (float *) MATRIZ_RESULTADO1; //res1 es un buffer de resultados intermedios chico (NSUBxMSUB)
        res1Ext = e_get_global_address(filaM,colM,res1);
        res2 = (float *) MATRIZ_RESULTADO2; //res2 funciona como acumulador grande de resultados intermedios, grande (NxM)
        res2Ext = e_get_global_address(filaM,colM,res2);

        //Inicializo los buffers-----------------
        //Acá se podría poner un "if" comando = inicio. En caso contrario acumular con lo previo.
        int *com;
        com = (int *)COMANDO;
        //e_barrier_init(barriers,tgt_bars);
        int i;
        if(*com == 0 || *com == 3){ //Iteración inicial o pasada única
            for(i=0;i<MSUB*NSUB;i++){
                *(res1+i) = 0;
            }
            for(i=0;i<N*N_FILAS/CORES;i++){
                *(res2+i) = 0;
            }
        }
        else{
            for(i=0;i<MSUB*NSUB;i++){
                *(res1+i) = 0;
            }
        }
        //e_barrier(barriers,tgt_bars);
        //---------------------------------------

        int *selector;
        selector = (int *)LUGAR_DE_ENTRADA;
        float *aBase, *bBase;

        float *aInternoBase, *bInternoBase;

        aInternoBase = (float *) VEC_COL_A;
        bInternoBase = (float *) VEC_FILA_B;

        if(*selector == 0){
            aBase = (float *) (BASE_ERAM_EPI + OFFSET_A);
            bBase = (float *) (BASE_ERAM_EPI + OFFSET_B);
        }
        else{
            aBase = (float *) (BASE_ERAM_EPI + OFFSET_A_2);
            bBase = (float *) (BASE_ERAM_EPI + OFFSET_B_2);
        }
        //Alineo a la zona de este nucleo
        aBase += coreNum*(KSUB/CORES)*N_FILAS;
        bBase += coreNum*(KSUB/CORES)*N;
        e_dma_copy(aInternoBase,aBase,(KSUB/CORES)*N_FILAS*sizeof(float));
        e_dma_copy(bInternoBase,bBase,(KSUB/CORES)*N*sizeof(float));

        e_barrier(barriers,tgt_bars);

        int filaTemp, colTemp; //Esto sirve para saber el resultado de qué nucleo estoy calculando
        int core;

int ssFila, ssCol;
float *aBaseTemp,*bBaseTemp,*res2Temp;

//bBaseTemp = bInternoBase;
res2Temp = res2;
for(ssCol=0;ssCol<(N/CORES)/NSUB;ssCol++){
    //aBaseTemp = aInternoBase;
    for(ssFila=0;ssFila<N_FILAS/MSUB;ssFila++){
        coreAnt(fila,col,&filaTemp,&colTemp); //Pongo el bloque inicial para calcular como el del núcleo anterior
//Esto hace un ciclo de respuesta--------------------------------------------------------------------------------------------------
        for(core=0;core<CORES;core++){
            //Consigo a y b
            //a = aBaseTemp;
            a = (aInternoBase + ssFila*MSUB);
            //b = bBaseTemp;
            b = (bInternoBase + ssCol*NSUB + (filaTemp+LADO_EPI*colTemp)*(N/CORES));
            //b += (N/CORES)*(filaTemp+LADO_EPI*colTemp);
            //b += (N/CORES)*coreNum;
            //b += NSUB*colTemp;

            if(core == CORES-1){ //En este caso quiero el resultado final en memoria externa
                if(*com == 2 || *com == 3){ //Iteración final o pasada única
                    //subMatmul(a,b,res1,&resEx[(*offsetResp) + coreNum*NSUB*NSUB],NSUB); //Asumo que CORES es par, entonces CORE-1 es impar
                    subMatmul(a,b,res1,res2Temp);
                }
                else{
                    subMatmul(a,b,res1,res2Ext);
                }
            }
            else{
                if(core%2==0){
                    subMatmul(a,b,res2Temp,res1Ext); //Resultados intermedios en res2, finales en res1 remoto
                }
                else{
                    subMatmul(a,b,res1,res2Ext); //Resultados intermedios en res1, finales en res2 remoto
                }

                //Ahora actualizo, hacia atrás la fila y columna temp (voy a calcular lo del núcleo anterior)
                coreAnt(filaTemp,colTemp,&filaTemp,&colTemp);
            }

        }
//-----------Fin ciclo de respuesta------------------------------------------------------------------------------------------------------------
        //aBaseTemp += MSUB;
        res2Temp += MSUB*NSUB; //Si MSUB<M entonces NSUB tiene que ser 1 (en principio, por lo menos)
        res2Ext = e_get_global_address(filaM,colM,res2Temp);
}
    //bBaseTemp += NSUB;
}
    if(*com == 2 || *com == 3){
        copiarDatos(res2,&resEx[coreNum*(N*N_FILAS/CORES)]);
    }
        //-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o


        //Sincronizo y señalo el fin<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->
        e_barrier(barriers, tgt_bars);

        if(coreNum == 0){
            *estado = 0;
        }
        //<-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><->

 //   __asm__ __volatile__("idle");
    }

}

void coreMatch(int fila,int col, int *filaM,int *colM){
    if(fila == 3){
        *filaM = 0;
        if(col == 3){
            *colM = 0;
        }
        else{
            *colM = col + 1;
        }
    }
    else{
        *colM = col;
        *filaM = fila + 1;
    }
}

void coreAnt(int fila,int col,int *filaM,int *colM){
    if(fila == 0){
        *filaM = 3;
        if(col == 0){
            *colM = 3;
        }
        else{
            *colM = col - 1;
        }
    }
    else{
        *colM = col;
        *filaM = fila - 1;
    }
}


void coords(int coreid,int *fila,int *col){
    int temp = coreid - 0x808;
    *fila = temp / 0x40;
    *col = (temp % 0x40);
}


void subMatmul(float *aIni,float *bIni, float *resInterno, float *resExternoIni){
    float *a, *b;
    int k,l,i,j;

    //volatile e_barrier_t barriers[16];
    //e_barrier_t *tgt_bars[16];

    e_barrier_init(barriers,tgt_bars);

    matmul_assembly(bIni,aIni,resInterno,resExternoIni);
/*
#define kASM 1
    for(k=0;k<(KSUB/CORES)/kASM;k++){
        a = (float *) aIni;
        a += N_FILAS*k*kASM;
        b = (float *) bIni;
        b += N*k*kASM;

        //En la última iteración mando al siguiente core
        if( k==(((KSUB/CORES)/kASM)-1) ){//Esto está mal
            //Esta version de matmul_assembly tiene en cuenta todos los stride, que A y C van por columnas, y B por filas
            for(j=0;j<NSUB;j++){
                for(l=0;l<MSUB;l++){
                    *(resExternoIni+l+j*MSUB) = (*(resInterno+l+j*MSUB)) + (*(a+l)) * (*(b+j));
                }
            }
            //matmul_assembly(b,a,resInterno,resExternoIni);
        }
        else{
            //Esta version de matmul_assembly tiene en cuenta todos los stride, que A y C van por columnas, y B por filas
            for(j=0;j<NSUB;j++){
                for(l=0;l<MSUB;l++){
                    *(resInterno+l+j*MSUB) = (*(resInterno+l+j*MSUB)) + (*(a+l)) * (*(b+j));
                }
            }
            //matmul_assembly(b,a,resInterno,resInterno);
        }

    }
*/
    e_barrier(barriers,tgt_bars);
    //--------------------------------------------------------
}


void copiarDatos(float *fuente, float *destino){
    float *fuenteTemp, *destinoTemp;
    int i,j;

    e_dma_copy(destino,fuente,(N*N_FILAS/CORES)*sizeof(float));

/*
    fuenteTemp = fuente;
    destinoTemp = destino;

    for(j=0; j<N/CORES; j++) {
        for(i=0; i<N_FILAS; i++) {
            *(destinoTemp+i) = *(fuenteTemp+i);
        }

        //e_dma_copy(destinoTemp,fuenteTemp,NSUB*sizeof(float));
        destinoTemp += N_FILAS;
        //fuenteTemp += NSUB;
        fuenteTemp += N_FILAS;
    }
*/
}
