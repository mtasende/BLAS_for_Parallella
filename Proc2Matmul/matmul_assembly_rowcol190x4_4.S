; Per-core Matrix multiplication
; Author: Anish Varghese <anish.varghese@anu.edu.au>, May 2014
; License: GPL v3

; Modified by Miguel Tasende <mtasendebracco@antel.com.uy>, Feb 2016
; Change: k = 8 instead of k = 32
; A, C stored by columns. B stored by rows.

; r0: pointer to A;
; r1: pointer to B;
; r2: pointer to C_in;
; r3: pointer to C_out
; r11,r12,r14,r15 - 4 points of matrix a
; r7 - value 0
; r13: stack pointer;
; r16 - r23: 8 points of Matrix B
; r24: NSUB
; r25: MSUB/strideC
; r32 - r63: intermediate results for Matrix C;

    .file   "matmul_assembly.S"
    .section    .text
    .balign 4
    .global _matmul_assembly
    .set stride, 192
    .set strideA, 1
    .set strideB, 256
    .set strideC, 32

    .macro doMult areg,index,aprev,inc
    fmadd r32,r\areg,r16
    ldrd r22,[r1,#\index + 3]
    fmadd r33,r\areg,r17
    ldr r\aprev,[r0,#\inc]
    fmadd r34,r\areg,r18
    fmadd r35,r\areg,r19
    fmadd r36,r\areg,r20
    fmadd r37,r\areg,r21
    ldrd r16,[r1,#\index + 4]
    fmadd r38,r\areg,r22
    ldrd r18,[r1,#\index + 5]
    fmadd r39,r\areg,r23
    ldrd r20,[r1,#\index + 6]
    fmadd r40,r\areg,r16
    ldrd r22,[r1,#\index + 7]
    fmadd r41,r\areg,r17
    fmadd r42,r\areg,r18
    fmadd r43,r\areg,r19
    fmadd r44,r\areg,r20
    fmadd r45,r\areg,r21
    ldrd r16,[r1,#\index + 8]
    fmadd r46,r\areg,r22
    ldrd r18,[r1,#\index + 9]
    fmadd r47,r\areg,r23
    ldrd r20,[r1,#\index + 10]
    fmadd r48,r\areg,r16
    ldrd r22,[r1,#\index + 11]
    fmadd r49,r\areg,r17
    fmadd r50,r\areg,r18
    fmadd r51,r\areg,r19
    fmadd r52,r\areg,r20
    fmadd r53,r\areg,r21
    ldrd r16,[r1,#\index + 12]
    fmadd r54,r\areg,r22
    ldrd r18,[r1,#\index + 13]
    fmadd r55,r\areg,r23
    ldrd r20,[r1,#\index + 14]
    fmadd r56,r\areg,r16
    ldrd r22,[r1,#\index + 15]
    fmadd r57,r\areg,r17
    fmadd r58,r\areg,r18
    fmadd r59,r\areg,r19
    fmadd r60,r\areg,r20
    fmadd r61,r\areg,r21
    ldrd r16,[r1,#\index + stride/2]
    fmadd r62,r\areg,r22
    ldrd r18,[r1,#\index + stride/2 + 1]
    fmadd r63,r\areg,r23
    ldrd r20,[r1,#\index + stride/2 + 2]
    .endm

    .macro doMultend areg,index,aprev,inc
    fmadd r32,r\areg,r16
    ldrd r22,[r1,#\index + 3]
    fmadd r33,r\areg,r17
    ldr r\aprev,[r0,#\inc]
    fmadd r34,r\areg,r18
    fmadd r35,r\areg,r19
    fmadd r36,r\areg,r20
    fmadd r37,r\areg,r21
    ldrd r16,[r1,#\index + 4]
    fmadd r38,r\areg,r22
    ldrd r18,[r1,#\index + 5]
    fmadd r39,r\areg,r23
    ldrd r20,[r1,#\index + 6]
    fmadd r40,r\areg,r16
    ldrd r22,[r1,#\index + 7]
    fmadd r41,r\areg,r17
    fmadd r42,r\areg,r18
    fmadd r43,r\areg,r19
    fmadd r44,r\areg,r20
    fmadd r45,r\areg,r21
    ldrd r16,[r1,#\index + 8]
    fmadd r46,r\areg,r22
    ldrd r18,[r1,#\index + 9]
    fmadd r47,r\areg,r23
    ldrd r20,[r1,#\index + 10]
    fmadd r48,r\areg,r16
    ldrd r22,[r1,#\index + 11]
    fmadd r49,r\areg,r17
    fmadd r50,r\areg,r18
    fmadd r51,r\areg,r19
    fmadd r52,r\areg,r20
    fmadd r53,r\areg,r21
    ldrd r16,[r1,#\index + 12]
    fmadd r54,r\areg,r22
    ldrd r18,[r1,#\index + 13]
    fmadd r55,r\areg,r23
    ldrd r20,[r1,#\index + 14]
    fmadd r56,r\areg,r16
    ldrd r22,[r1,#\index + 15]
    fmadd r57,r\areg,r17
    ;Change r1 to point to next sub-column of "B" (A).
    add r1,r1,#(strideC * 4)
    fmadd r58,r\areg,r18
    fmadd r59,r\areg,r19
    fmadd r60,r\areg,r20
    fmadd r61,r\areg,r21
    ;Changes here to point back to first row
    ldrd r16,[r1,#0]
    fmadd r62,r\areg,r22
    ldrd r18,[r1,#1]
    fmadd r63,r\areg,r23
    ldrd r20,[r1,#2]
    .endm

    .macro doMultincr areg,index,aprev,inc
    fmadd r32,r\areg,r16
    ldrd r22,[r1,#\index + 3]
    fmadd r33,r\areg,r17
    ldr r\aprev,[r0,#\inc]
    fmadd r34,r\areg,r18
    ;Move r0 to point to next row of Matrix A - only change
    add r0,r0,#(strideA * 4)
    fmadd r35,r\areg,r19
    fmadd r36,r\areg,r20
    fmadd r37,r\areg,r21
    ldrd r16,[r1,#\index + 4]
    fmadd r38,r\areg,r22
    ldrd r18,[r1,#\index + 5]
    fmadd r39,r\areg,r23
    ldrd r20,[r1,#\index + 6]
    fmadd r40,r\areg,r16
    ldrd r22,[r1,#\index + 7]
    fmadd r41,r\areg,r17
    fmadd r42,r\areg,r18
    fmadd r43,r\areg,r19
    fmadd r44,r\areg,r20
    fmadd r45,r\areg,r21
    ldrd r16,[r1,#\index + 8]
    fmadd r46,r\areg,r22
    ldrd r18,[r1,#\index + 9]
    fmadd r47,r\areg,r23
    ldrd r20,[r1,#\index + 10]
    fmadd r48,r\areg,r16
    ldrd r22,[r1,#\index + 11]
    fmadd r49,r\areg,r17
    fmadd r50,r\areg,r18
    fmadd r51,r\areg,r19
    fmadd r52,r\areg,r20
    fmadd r53,r\areg,r21
    ldrd r16,[r1,#\index + 12]
    fmadd r54,r\areg,r22
    ldrd r18,[r1,#\index + 13]
    fmadd r55,r\areg,r23
    ldrd r20,[r1,#\index + 14]
    fmadd r56,r\areg,r16
    ldrd r22,[r1,#\index + 15]
    fmadd r57,r\areg,r17
    fmadd r58,r\areg,r18
    fmadd r59,r\areg,r19
    fmadd r60,r\areg,r20
    fmadd r61,r\areg,r21
    ldrd r16,[r1,#\index + stride/2]
    fmadd r62,r\areg,r22
    ldrd r18,[r1,#\index + stride/2 + 1]
    fmadd r63,r\areg,r23
    ldrd r20,[r1,#\index + stride/2 + 2]
    .endm


_matmul_assembly:
; save the 26 registers (13 dual) we are using that need callee saving
    add sp,sp,#-120
    strd r0,[sp,#0]
    strd r2,[sp,#1]
    strd r4,[sp,#2]
    strd r6,[sp,#3]
    strd r8,[sp,#4]
    strd r10,[sp,#5]
    strd r14,[sp,#6]
    strd r28,[sp,#7]
    strd r30,[sp,#8]
    strd r32,[sp,#9]
    strd r34,[sp,#10]
    strd r36,[sp,#11]
    strd r38,[sp,#12]
    strd r40,[sp,#13]
    strd r42,[sp,#14]

    ;load NSUB = 4
    mov r24,#4

    ;Matrix C - Reading first row of Intermediate result Matrix C
    ldrd r32,[r2,#0]
    ldrd r34,[r2,#1]
    ldrd r36,[r2,#2]
    ldrd r38,[r2,#3]
    ldrd r40,[r2,#4]
    ldrd r42,[r2,#5]
    ldrd r44,[r2,#6]
    ldrd r46,[r2,#7]
    ldrd r48,[r2,#8]
    ldrd r50,[r2,#9]
    ldrd r52,[r2,#10]
    ldrd r54,[r2,#11]
    ldrd r56,[r2,#12]
    ldrd r58,[r2,#13]
    ldrd r60,[r2,#14]
    ldrd r62,[r2,#15]

.Lb:
    ;load (M/MSUB)=6
    mov r25,#6

    ;preload the first row and first row of A and B
    ;Matrix A
    ldr r11,[r0,#0]
    ldr r12,[r0,#1*strideB]
    ldr r14,[r0,#2*strideB]

    ;Matrix B
    ldrd r16,[r1,#0]
    ldrd r18,[r1,#1]
    ldrd r20,[r1,#2]

.loopInt:
    doMult 11,0,15,3*strideB
    doMult 12,1*stride/2,11,0
    doMult 14,2*stride/2,12,strideB
    doMultend 15,3*stride/2,14,2*strideB

    ;sacar datos a C
    strd r32,[r3,#0]
    ldrd r32,[r2,#((strideC/2) + 0)]
    strd r34,[r3,#1]
    ldrd r34,[r2,#((strideC/2) + 1)]
    strd r36,[r3,#2]
    ldrd r36,[r2,#((strideC/2) + 2)]
    strd r38,[r3,#3]
    ldrd r38,[r2,#((strideC/2) + 3)]
    strd r40,[r3,#4]
    ldrd r40,[r2,#((strideC/2) + 4)]
    strd r42,[r3,#5]
    ldrd r42,[r2,#((strideC/2) + 5)]
    strd r44,[r3,#6]
    ldrd r44,[r2,#((strideC/2) + 6)]
    strd r46,[r3,#7]
    ldrd r46,[r2,#((strideC/2) + 7)]
    strd r48,[r3,#8]
    ldrd r48,[r2,#((strideC/2) + 8)]
    strd r50,[r3,#9]
    ldrd r50,[r2,#((strideC/2) + 9)]
    strd r52,[r3,#10]
    ldrd r52,[r2,#((strideC/2) + 10)]
    strd r54,[r3,#11]
    ldrd r54,[r2,#((strideC/2) + 11)]
    strd r56,[r3,#12]
    ldrd r56,[r2,#((strideC/2) + 12)]
    strd r58,[r3,#13]
    ldrd r58,[r2,#((strideC/2) + 13)]
    strd r60,[r3,#14]
    ldrd r60,[r2,#((strideC/2) + 14)]
    strd r62,[r3,#15]
    ldrd r62,[r2,#((strideC/2) + 15)]
    ;Point to next row of Matrix C_in and C_out
    add r2,r2,#(strideC * 4)
    add r3,r3,#(strideC * 4)

    doMult 11,0,15,3*strideB
    doMult 12,1*stride/2,11,0
    doMult 14,2*stride/2,12,strideB
    doMultend 15,3*stride/2,14,2*strideB

    ;sacar datos a C
    strd r32,[r3,#0]
    ldrd r32,[r2,#((strideC/2) + 0)]
    strd r34,[r3,#1]
    ldrd r34,[r2,#((strideC/2) + 1)]
    strd r36,[r3,#2]
    ldrd r36,[r2,#((strideC/2) + 2)]
    strd r38,[r3,#3]
    ldrd r38,[r2,#((strideC/2) + 3)]
    strd r40,[r3,#4]
    ldrd r40,[r2,#((strideC/2) + 4)]
    strd r42,[r3,#5]
    ldrd r42,[r2,#((strideC/2) + 5)]
    strd r44,[r3,#6]
    ldrd r44,[r2,#((strideC/2) + 6)]
    strd r46,[r3,#7]
    ldrd r46,[r2,#((strideC/2) + 7)]
    strd r48,[r3,#8]
    ldrd r48,[r2,#((strideC/2) + 8)]
    strd r50,[r3,#9]
    ldrd r50,[r2,#((strideC/2) + 9)]
    strd r52,[r3,#10]
    ldrd r52,[r2,#((strideC/2) + 10)]
    strd r54,[r3,#11]
    ldrd r54,[r2,#((strideC/2) + 11)]
    strd r56,[r3,#12]
    ldrd r56,[r2,#((strideC/2) + 12)]
    strd r58,[r3,#13]
    ldrd r58,[r2,#((strideC/2) + 13)]
    strd r60,[r3,#14]
    ldrd r60,[r2,#((strideC/2) + 14)]
    strd r62,[r3,#15]
    ldrd r62,[r2,#((strideC/2) + 15)]
    ;Point to next row of Matrix C_in and C_out
    add r2,r2,#(strideC * 4)
    add r3,r3,#(strideC * 4)

    ;Loop interno (hay loop unrolling = 2)
    sub r25,r25,#2
    nop
    bne .loopInt

    ;Volver al inicio de la matriz B (A)
    sub r1,r1,#(stride*4)

    ;Aumentar 1 en la matriz A (B)
    add r0,r0,#(strideA * 4)

    ;Now loop back
    sub r24,r24,#1
    nop
    bne .Lb

; restore 26 saved registers
    ldrd r0,[sp,#0]
    ldrd r2,[sp,#1]
    ldrd r4,[sp,#2]
    ldrd r6,[sp,#3]
    ldrd r8,[sp,#4]
    ldrd r10,[sp,#5]
    ldrd r14,[sp,#6]
    ldrd r28,[sp,#7]
    ldrd r30,[sp,#8]
    ldrd r32,[sp,#9]
    ldrd r34,[sp,#10]
    ldrd r36,[sp,#11]
    ldrd r38,[sp,#12]
    ldrd r40,[sp,#13]
    ldrd r42,[sp,#14]
    add sp,sp,#120
    rts
    .size   _matmul_assembly, .-_matmul_assembly
    .ident  "Anish V, v0.1"
    .end
