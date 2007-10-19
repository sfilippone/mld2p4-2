#!/bin/sh

# 
#
# NUmber of attempts for each configuration
ntry=1


date=`date +%Y%m%d%H%M%S`



for np in $* 
do

# 3rd  batch: 14bis,64bis  4 sweeps
echo "mpirun -np $np -machinefile locm df_bench >>log.part$part.ren$renum.${np}p"
#/usr/local/mpich-gcc42/bin/mpirun -np $np -machinefile locm   df_bench  <<EOF
/usr/local/mpich-gcc42/bin/mpirun -np $np -machinefile locm   df_bench  >>run.gfc.kiva.log.${np}p.$date 2>err.gfc.kiva.log.${np}p.$date <<EOF
out.${np}p.$date               Out file 1: summary
stat.${np}p.$date               Out file 2: detailed for statistics
BICGSTAB            iterative method to use
1.D-6               EPS
CSR                 Matrix format
0                   IPART: Partition method: 0: BLOCK  1: BLK2 2:GRAPH
00250               ITMAX
-1                  ITRACE
2                   ISTOPC   1: NBE Infty   2: |r|2/|b|2
30                  IRST   Restart parameter for GMRES and BiCGSTAB(L)
0                   RENUM: 0: none 1: global indices (2: GPS band reduction)
$ntry               NTRY  for each comb. print out best timings
30    NPRCS  nov rst prl fc1 fl1 mlt agg smt cm smp ft2 fl2 jsw nl  omg   th1  th2  name
none   none   0   0   0   0   0   0   0   0  0   0   0   0   0   1  -1.0 1e-4 1e-4  NOPREC
diag   none   0   0   0   1   0   2   0   1  0   2   1   0   4   1  -1.0 1e-4 1e-4  DIAG
bjac   none   0   0   0   1   0   2   0   1  0   2   1   0   4   1  -1.0 1e-4 1e-4  BJAC  
as     none   0   1   0   1   0   2   0   1  0   2   1   0   4   1  -1.0 1e-4 1e-4  RAS  
as     none   1   1   0   1   0   2   0   1  0   2   1   0   4   1  -1.0 1e-4 1e-4  RAS  
as     none   2   1   0   1   0   2   0   1  0   2   1   0   4   1  -1.0 1e-4 1e-4  RAS  
as     ml     0   1   0   1   0   2   0   1  0   2   1   0   4   2  -1.0 1e-4 1e-4  2L-M-RAS-I-D4
as     ml     1   1   0   1   0   2   0   1  0   2   1   0   4   2  -1.0 1e-4 1e-4  2L-M-RAS-I-D4
as     ml     2   1   0   1   0   2   0   1  0   2   1   0   4   2  -1.0 1e-4 1e-4  2L-M-RAS-I-D4
as     ml     0   1   0   1   0   2   0   1  0   2   5   0   4   2  -1.0 1e-4 1e-4  2L-M-RAS-U-D4
as     ml     1   1   0   1   0   2   0   1  0   2   5   0   4   2  -1.0 1e-4 1e-4  2L-M-RAS-U-D4
as     ml     2   1   0   1   0   2   0   1  0   2   5   0   4   2  -1.0 1e-4 1e-4  2L-M-RAS-U-D4
as     ml     0   1   0   1   0   2   0   1  0   2   1   0   4   3  -1.0 1e-4 1e-4  3L-M-RAS-I-D4
as     ml     1   1   0   1   0   2   0   1  0   2   1   0   4   3  -1.0 1e-4 1e-4  3L-M-RAS-I-D4
as     ml     2   1   0   1   0   2   0   1  0   2   1   0   4   3  -1.0 1e-4 1e-4  3L-M-RAS-I-D4
as     ml     0   1   0   1   0   2   0   1  0   2   5   0   4   3  -1.0 1e-4 1e-4  3L-M-RAS-U-D4
as     ml     1   1   0   1   0   2   0   1  0   2   5   0   4   3  -1.0 1e-4 1e-4  3L-M-RAS-U-D4
as     ml     2   1   0   1   0   2   0   1  0   2   5   0   4   3  -1.0 1e-4 1e-4  3L-M-RAS-U-D4
as     ml     0   1   0   1   0   2   0   1  1   2   1   0   1   2  -1.0 1e-4 1e-4  2L-M-RAS-I-R
as     ml     1   1   0   1   0   2   0   1  1   2   1   0   1   2  -1.0 1e-4 1e-4  2L-M-RAS-I-R
as     ml     2   1   0   1   0   2   0   1  1   2   1   0   1   2  -1.0 1e-4 1e-4  2L-M-RAS-I-R
as     ml     0   1   0   1   0   2   0   1  1   2   5   0   1   2  -1.0 1e-4 1e-4  2L-M-RAS-U-R
as     ml     1   1   0   1   0   2   0   1  1   2   5   0   1   2  -1.0 1e-4 1e-4  2L-M-RAS-U-R
as     ml     2   1   0   1   0   2   0   1  1   2   5   0   1   2  -1.0 1e-4 1e-4  2L-M-RAS-U-R
as     ml     0   1   0   1   0   2   0   1  1   2   1   0   1   3  -1.0 1e-4 1e-4  3L-M-RAS-I-R
as     ml     1   1   0   1   0   2   0   1  1   2   1   0   1   3  -1.0 1e-4 1e-4  3L-M-RAS-I-R
as     ml     2   1   0   1   0   2   0   1  1   2   1   0   1   3  -1.0 1e-4 1e-4  3L-M-RAS-I-R
as     ml     0   1   0   1   0   2   0   1  1   2   5   0   1   3  -1.0 1e-4 1e-4  3L-M-RAS-U-R
as     ml     1   1   0   1   0   2   0   1  1   2   5   0   1   3  -1.0 1e-4 1e-4  3L-M-RAS-U-R
as     ml     2   1   0   1   0   2   0   1  1   2   5   0   1   3  -1.0 1e-4 1e-4  3L-M-RAS-U-R
2                  Number of matrices
kivap004.mtx            none                                     
kivap001.mtx            none                                     
thm50x30.mtx            none                                     
thm200x120.mtx            none                                     
thm1000x600.mtx            none                                     
a400x400.mtx     b400.mtx
kivap007.mtx            none                                     
!!! preconditioner templates
bja       none      0  0  0  0  0  0  0  0  0  0  0  0.0     Block Jacobi
none      none      0  0  0  0  0  0  0  0  0  0  0  0.1     No preconditioner
diagsc    none      0  0  0  0  0  0  0  0  0  0  0  0.0     Diagonal scaling
as        none      1  4  1  0  0  0  0  0  0  0  0  0.0     Additive Schwarz 1 overlap
as        none      1  4  0  0  0  0  0  0  0  0  0  0.0     Restricted Additive Schwarz 1 overlap
EOF

cat out.${np}p.$date >>dat.out.kiva.$date
cat stat.${np}p.$date >>dat.stat.kiva.$date

done


