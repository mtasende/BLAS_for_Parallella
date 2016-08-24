/* This code was implemented by Miguel Tasende. Copyright (C) 2016, Antel (Uruguay).*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#include <blis.h>
#include <e-hal.h>

void bli_dgemm_opt_mxn(
                        dim_t              k,
                        double*   restrict alpha,
                        double*   restrict a1,
                        double*   restrict b1,
                        double*   restrict beta,
                        double*   restrict c11, inc_t rs_c, inc_t cs_c,
                        auxinfo_t*         data
                      );

void bli_sgemm_opt_mxn(
                        dim_t              k,
                        float*   restrict alpha,
                        float*   restrict a1,
                        float*   restrict b1,
                        float*   restrict beta,
                        float*   restrict c11, inc_t rs_c, inc_t cs_c,
                        auxinfo_t*         data
                      );

void bli_sgemm_opt_mxn_noInit(int k,
                              float* restrict alpha,
                              float* restrict a1,
                              float* restrict b1,
                              float* restrict beta,
                              float* restrict c11,
                              int rs_c,
                              int cs_c,
                              e_epiphany_t *dev,
                              e_mem_t *memory,
                              e_platform_t *platform,
                              int mBlis,
                              int nBlis,
                              int cs_aBlis,
                              int rs_bBlis
                              );

