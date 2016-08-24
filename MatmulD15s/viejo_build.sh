#!/bin/bash

set -e

ESDK=${EPIPHANY_HOME}
ELIBS="-L ${ESDK}/tools/host/lib -L /home/parallella/linpack/blisReference/lib"
EINCS="-I ${ESDK}/tools/host/include -I /home/parallella/linpack/blisReference/include/blis"
#ELDF=./matmulB.ldf
ELDF=${ESDK}/bsps/current/internal.ldf
#ELDF=${ESDK}/bsps/current/fast.ldf


SCRIPT=$(readlink -f "$0")
EXEPATH=$(dirname "$SCRIPT")
cd $EXEPATH

# Create the binaries directory
mkdir -p bin/

CROSS_PREFIX=
case $(uname -p) in
	arm*)
		# Use native arm compiler (no cross prefix)
		CROSS_PREFIX=
		;;
	   *)
		# Use cross compiler
		CROSS_PREFIX="arm-linux-gnueabihf-"
		;;
esac

# Build HOST side application
${CROSS_PREFIX}gcc main.c bli_gemm_opt_mxn.c -o bin/main.elf ${EINCS} ${ELIBS} -le-hal -le-loader -lpthread -lblis -std=c99

# Build DEVICE side program
#e-gcc -T ${ELDF} e_task.c -o bin/e_task.elf -le-lib -lm
#e-gcc -O3  -T ${ELDF} device.c -o bin/device.elf -le-lib -lm -ffast-math
e-gcc -O2  -T ${ELDF} device.c -o bin/device.elf -le-lib -lm -ffast-math
#e-gcc -Wl,--print-map  -O3  -T ${ELDF} e_task.c -o bin/e_task.elf -le-lib -lm -ffast-math
#e-gcc -c -O3  e_task.c -o e_task.o
#e-ld  ${ELIBS} -T ${ELDF} e_task.o -o bin/e_task.elf

# Convert ebinary to SREC file
e-objcopy --srec-forceS3 --output-target srec bin/device.elf bin/device.srec
