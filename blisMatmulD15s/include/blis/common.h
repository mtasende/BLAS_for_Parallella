//#define DEBUG_DEL_MATMUL 1
#define CORES 16
#define LADO_EPI 4
#define N_FILAS 192
#define N 256
#define K 2048//1280//Multiplo de KSUB
#define NSUB 4
#define MSUB 192 //Múltiplo de 64 o coordinar con Matmul_assembly
#define KSUB 64 //tiene que ser multiplo de CORES, es la cantidad de vectores en la dimensión k que mando por carga de chip (limitado por el espacio de entrada)

//Entrada
#define BASE_ERAM_EPI 0x8F000000
#define OFFSET_A 0x00200000
#define OFFSET_A_2 0x00300000
#define OFFSET_B 0x00400000
#define OFFSET_B_2 0x00500000
#define OFFSET_A_INI 0x00600000
#define OFFSET_B_INI 0x00700000

#define VEC_COL_A 0x6000
#define VEC_FILA_B 0x7000
//#define VEC_COL_A_2 0x7000
//#define VEC_FILA_B_2 0x7800


//Control
/*
#define LISTO 0x3800
#define LISTO_MUL_PREVIO 0x3810
#define RESULTADO_EXTERNO 0x3820
#define RESIDUO 0x3830
#define COMPARTIDA_DEV 0x3840
*/
#define LISTO 0x5c00
#define LISTO_MUL_PREVIO 0x5c10
#define RESULTADO_EXTERNO 0x5c20
#define RESIDUO 0x5c30
#define COMPARTIDA_DEV 0x5c40
#define COMANDO 0x5c50 //Esto es para mandarle info al Epiphany (inicio/intermedio/fin/otros)
#define LUGAR_DE_ENTRADA 0x5c60

//Debug
#define MENSAJES 0x5c50
#define TAM_BUFFER 128

//Resultado Temp
#define MATRIZ_RESULTADO1 0x5000
#define MATRIZ_RESULTADO2 0x2000
