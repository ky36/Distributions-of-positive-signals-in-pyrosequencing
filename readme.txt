On this site you can find the source code of the simulation program 
sbs_ffcm_simul.c to simulate the distribution of sequence length for 
sequencing by synthesis (SBS), and compare the simulation results 
with the analytical results.  The program is written in C programming 
language. To compile the program, use something like this 
(different platforms/compilers differ slightly):

gcc -Wall -o sbs_ffcm_simul sbs_ffcm_simul.c 

I also include two compiled binary (executable) files, one for 32-bit computers
(sbs_ffcm_simul_32), one for 64-bit computers (sbs_ffcm_simul_64).

To use the program:

./sbs_ffcm_simul -f <flow cycles> -r <repeats> -a <prob of a> -b <prob of b> -c <prob of c> -d <prob of g> -A a0,a1,... -B b0,b1,...  -C c0,c1,... -D d0,d1,...

The parameters are:

-f: number of flow cycles
-r: number of simulations
-a: nucleotide context probability of "a"
-b: nucleotide context probability of "b"
-c: nucleotide context probability of "c"
-d: nucleotide context probability of "d"
-A: a list of nucleotide incorporation probabilities of nucleotide "a"
-B: a list of nucleotide incorporation probabilities of nucleotide "b"
-C: a list of nucleotide incorporation probabilities of nucleotide "c"
-D: a list of nucleotide incorporation probabilities of nucleotide "d"

For example,

./sbs_ffcm_simul -f 100 -r 200000 -a 0.3333 -b 0.0909 -c 0.4329 -d 0.1429 -A .1090909091,.5,.3,.09090909091 -B .3166666667,.25,.3333333333,.1 -C .646031
7460,.1428571429,.1,.1111111111 -D .425,.2,.25,.125


Yong Kong
Yale University
yong.kong@yale.edu


