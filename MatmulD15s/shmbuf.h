/* This code was implemented by Miguel Tasende. Copyright (C) 2016, Antel (Uruguay).*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef SHMBUF_H
#define SHMBUF_H

#include <sys/mman.h>
#include <fcntl.h>
#include <semaphore.h>
#include <sys/stat.h>
#include "common.h"

#define KMAX 8192

struct shmbuf{
    sem_t semIni; //El semáforo de inicio
    sem_t semFin; //El semáforo de fin
    //Parámetros de sgemm----------------
    int k;
    float alpha;
    float a1[N_FILAS*KMAX];
    float b1[N*KMAX];
    float beta;
    float c11[N*N_FILAS];
    int rs_c;
    int cs_c;
    int mBlis;
    int nBlis;
    int cs_aBlis;
    int rs_bBlis;
    //-----------------------------------
};

#endif // SHMBUF_H


