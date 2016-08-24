A Single Precision BLAS Library for Parallella, with Epiphany acceleration. Generated using BLIS.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.


---------- BUILDING AND RUNNING ---------------------------------------------

First copy this folder (BLASv1) into your Parallella host.
In the rest of the document the folders are referenced relative to the BLASv1 folder.

.................................
1) MATRIX MULTIPLICATION KERNEL
a) BUILDING

Follow the commands below:
> cd MatmulD15s
> ./build.sh
> cd ../Proc2Matmul
> ./build.sh

b) RUNNING

i) Running in the same process
> cd MatmulD15s
> ./run.sh

ii) Running on another process
Open two terminals on the Parallella host ($1 and $2)

$2> cd Proc2Matmul
$2> ./run.sh

Wait until a message appears (in Spanish; may translate it later): "Esperando pedido de cálculo..."

$1> cd MatmulD15s
$1> ./runD.sh

This will send the data to the process that was initiated before, it will send it to the Epiphany, and after all the calculations, send it back.

.................................
2) BLAS (or BLIS)

i) Running the BLIS tests
Before running the BLIS tests it is necessary to initiate the secondary process (which will be called by the kernel). Remember to wait for the "Esperando pedido de cálculo..." message ("ready" message).
In a second terminal ($2) run:
$2> cd ProcMatmul
$2> ./run.sh

In the first terminal:
$1> cd blis-master/testsuite
$1> ./test_libblis.x

That will test the sgemm micro-kernel, and the complete sgemm function. If you want you can also test the dgemm function, but it will run in Single Precision ("false dgemm"), and many "FAILURE" notes will appear (the precision should be comparable to Single in any case). Many other options can be changed for the tests. Please read the BLIS documentation to learn more about it.


ii) Using the generated BLAS Library
You can add the generated BLAS library to your linker script in one of you programs, or to other scientific software like Octave, LAPACK, etc. The location of the generated BLAS library is:

./blisMatmulD15s/lib/libblis.a  (that's a link to the real library file on the same folder)

You may also want to use some of the header files in the ./blisMatmulD15s/include/blis


iii) Compiling the BLAS (BLIS) library
To compile the BLAS library with BLIS please first read the BLIS documentation.
The new configuration folder that was created for the Epiphany accelerated kernel can be found in:

./blis-master/config/matmulD15s

In the ./MatmulD15s folder (outside blis-master) there is a script that makes it easy to send the needed source code to the configuration folder (to make "fast" changes and update the config in BLIS). You may want to use it. It is called:

./MatmulD15s/sincroint


.................................
3) High-Performance Linpack Benchmark

i) Running the benchmark
To run the benchmark with the standard options follow the commands below:
First remember to initialize the secondary process in another terminal:

In a second terminal ($2) run:
$2> cd ProcMatmul
$2> ./run.sh

In the first terminal run:
$1> cd hpl-2.1/bin/blisMatmulD15s

Now please modify the file "hostSingle" changing the name to the one of your Parallella device (you can use "nano hostSingle" or do it with any text editor)

$1> ./runSingle

It takes about 132 seconds to run with the initial configuration I am distributing the software. If you want, you can run "top" in another terminal to see the two processes working "xhpl" and 
"matmuld.elf", and how the share the CPU time (matmuld.elf corresponds to the secondary process).

ii) Configuring and compiling
For Configuration of the benchmark you can change the HPL.dat file. Please read Netlib's HPL documentation.
If you want to compile the HPL code (possibly with changes), the configuration file used is in:

./hpl-2.1/Make.blisMatmulD15s

But be aware that some of the configuration options were specific to the file system that was used in this work (you may have to adapt them to yours). Among other things it is very possible that there are many non-relative folder paths.
