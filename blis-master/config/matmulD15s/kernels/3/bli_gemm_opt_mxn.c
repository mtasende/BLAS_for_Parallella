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
#include <string.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <time.h>

#include "shmbuf.h"

#define _BufOffset (0x01000000)

//#define DEBUGD

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
#ifdef DEBUGD
    clock_t inicio,fin;
    inicio = clock();
#endif // DEBUGD

    int fd;
    //Inicializo el buffer de memoria compartida
    struct shmbuf *shmp;
    fd = shm_open("/matmul-file",O_RDWR,S_IRUSR|S_IWUSR);
    //Creo el mapeo de memoria
    shmp = mmap(NULL,sizeof(struct shmbuf),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    close(fd);

    //Paso los valores a memoria compartida---------------------------------------------------------
    int mBlis = bli_smr;
	int nBlis = bli_snr;
	int cs_aBlis = bli_spackmr;
	int rs_bBlis = bli_spacknr;

	/*
	mBlis = N_FILAS;
	nBlis = N;
	cs_aBlis = N_FILAS;
	rs_bBlis = N;
	*/

    shmp->k = k;
    shmp->alpha = *alpha;
    shmp->beta = *beta;
    //memcpy(&shmp->a1[0],a1,N_FILAS*k*sizeof(float));
    int l,p;
    for(l=0;l<k;l++){
        memcpy(&shmp->a1[l*mBlis],(a1+l*cs_aBlis),mBlis*sizeof(float));
    }
    for(l=0;l<k;l++){
        memcpy(&shmp->b1[l*nBlis],(b1+l*rs_bBlis),nBlis*sizeof(float));
    }
    for(l=0;l<nBlis;l++){
        for(p=0;p<mBlis;p++){
            shmp->c11[l*mBlis+p] = *(c11 + l*cs_c + p*rs_c);
        }
    }
    //shmp->rs_c = rs_c;
    //shmp->cs_c = cs_c;
    shmp->rs_c = 1;
    shmp->cs_c = mBlis;

    shmp->mBlis = mBlis;
    shmp->nBlis = nBlis;
    //shmp->cs_aBlis = cs_aBlis;
    //shmp->rs_bBlis = rs_bBlis;
    shmp->cs_aBlis = mBlis;
    shmp->rs_bBlis = nBlis;
    //----------------------------------------------------------------------------------------------

    //Falta arreglar el tema del mutex, para cuando pueda haber varios procesos llamando a matmul
    sem_post(&shmp->semIni);//Inicio el cálculo

    sem_wait(&shmp->semFin);//Espero a que termine de calcular


    //Traigo los valores de la memoria compartida---------------------------------------------------
    //memcpy(c11,&shmp->c11[0],N*N_FILAS*sizeof(float));
    for(l=0;l<nBlis;l++){
        for(p=0;p<mBlis;p++){
            *(c11 + l*cs_c + p*rs_c) = shmp->c11[l*mBlis+p];
        }
    }
    //----------------------------------------------------------------------------------------------

    munmap(shmp,sizeof(struct shmbuf)); //Innecesario porque se libera cuando el programa termina

#ifdef DEBUGD
    //Calculo el tiempo
    fin = clock();
    double tiempo;
    tiempo = (double)(fin-inicio)/CLOCKS_PER_SEC;
    printf("Tiempo = %f -- %f GFLOPS/s\n\n",tiempo,N*N_FILAS*(2*K-1)/tiempo/1000000000);
#endif // DEBUGD
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
#ifdef DEBUGD
    clock_t inicio,fin;
    inicio = clock();
#endif // DEBUGD

    int fd;
    //Inicializo el buffer de memoria compartida
    struct shmbuf *shmp;
    fd = shm_open("/matmul-file",O_RDWR,S_IRUSR|S_IWUSR);
    //Creo el mapeo de memoria
    shmp = mmap(NULL,sizeof(struct shmbuf),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    close(fd);

    //Paso los valores a memoria compartida---------------------------------------------------------
    int mBlis = bli_dmr;
	int nBlis = bli_dnr;
	int cs_aBlis = bli_dpackmr;
	int rs_bBlis = bli_dpacknr;

	/*
	mBlis = N_FILAS;
	nBlis = N;
	cs_aBlis = N_FILAS;
	rs_bBlis = N;
	*/

    shmp->k = k;
    shmp->alpha = (float)(*alpha);
    shmp->beta = (float)(*beta);
    //memcpy(&shmp->a1[0],a1,N_FILAS*k*sizeof(float));
    int l,p;
    for(l=0;l<k;l++){
        for(p=0;p<mBlis;p++){
            shmp->a1[l*mBlis+p] = (float)(*(a1+l*cs_aBlis+p));
        }
    }
    for(l=0;l<k;l++){
        for(p=0;p<nBlis;p++){
            shmp->b1[l*nBlis+p] = (float)(*(b1+l*rs_bBlis+p));
        }
    }
    for(l=0;l<nBlis;l++){
        for(p=0;p<mBlis;p++){
            shmp->c11[l*mBlis+p] = (float)(*(c11 + l*cs_c + p*rs_c));
        }
    }
    //shmp->rs_c = rs_c;
    //shmp->cs_c = cs_c;
    shmp->rs_c = 1;
    shmp->cs_c = mBlis;

    shmp->mBlis = mBlis;
    shmp->nBlis = nBlis;
    //shmp->cs_aBlis = cs_aBlis;
    //shmp->rs_bBlis = rs_bBlis;
    shmp->cs_aBlis = mBlis;
    shmp->rs_bBlis = nBlis;
    //----------------------------------------------------------------------------------------------

    //Falta arreglar el tema del mutex, para cuando pueda haber varios procesos llamando a matmul
    sem_post(&shmp->semIni);//Inicio el cálculo

    sem_wait(&shmp->semFin);//Espero a que termine de calcular


    //Traigo los valores de la memoria compartida---------------------------------------------------
    //memcpy(c11,&shmp->c11[0],N*N_FILAS*sizeof(float));
    for(l=0;l<nBlis;l++){
        for(p=0;p<mBlis;p++){
            *(c11 + l*cs_c + p*rs_c) = (double)shmp->c11[l*mBlis+p];
        }
    }
    //----------------------------------------------------------------------------------------------

    munmap(shmp,sizeof(struct shmbuf)); //Innecesario porque se libera cuando el programa termina

#ifdef DEBUGD
    //Calculo el tiempo
    fin = clock();
    double tiempo;
    tiempo = (double)(fin-inicio)/CLOCKS_PER_SEC;
    printf("Tiempo = %f -- %f GFLOPS/s\n\n",tiempo,N*N_FILAS*(2*K-1)/tiempo/1000000000);
#endif // DEBUGD
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
