# Distributions of positive signals in pyrosequencing

This site contains information for the simulation program for
"Distributions of positive signals in pyrosequencing" (see reference [4]).

The name of the simulation program is `sbs_ffcm_r.c`.  The structure and
usage of the program is very similar to `sbs_ffcm_simul.c`, the program
to simulate sequence length distribution.  Actually `sbs_ffcm_r.c` includes
the functions of `sbs_ffcm_simul.c`, and more.  So please also refer to
the [] for more information.
 
To compile the program, use something like this:

gcc -o sbs_ffcm_r -W -Wall -O3 sbs_ffcm_r.c -lm

Two binary executables are also included: `sbs_ffcm_r_32` for 32-bit computers
and `sbs_ffcm_r_64` for 64-bit computers.

To run the program, 

`./sbs_ffcm_r -f 10 -r 2000000 -a 0.333333333 -b 0.090909091 -c 0.432900433 -d 0.142857143 -A 1 -B 1 -C 1 -D 1`


The file `sbs_ffcm_r_100_repeat_20000000_A0.333_1.00_B0.091_1.00_C0.433_1.00_D0.143_1.00` contains the output of
a sample simulation run.


## References
40
1. Kong Y. Statistical distributions of pyrosequencing. Journal of computational biology: a journal of computational molecular cell biology. 2009;16(1):31-42. Epub 2008/12/17. [doi: 10.1089/cmb.2008.0106](https://www.liebertpub.com/doi/10.1089/cmb.2008.0106) PubMed PMID: 19072582.
41
2. Kong Y. Statistical distributions of sequencing by synthesis with probabilistic nucleotide incorporation. Journal of computational biology: a journal of computational molecular cell biology. 2009;16(6):817-27. Epub 2009/06/16. [doi: 10.1089/cmb.2008.0215](https://www.liebertpub.com/doi/10.1089/cmb.2008.0215). PubMed PMID: 19522665.
42
3. Kong Y. Length distribution of sequencing by synthesis: fixed flow cycle model. J Math Biol. 2013;67(2):389-410. Epub 2012/06/13. [doi: 10.1007/s00285-012-0556-3](https://link.springer.com/article/10.1007/s00285-012-0556-3). PubMed PMID: 22689207.
43
4. Kong Y. Distributions of positive signals in pyrosequencing. J Math Biol. 2014;69(1):39-54. Epub 2013/06/01. [doi: 10.1007/s00285-013-0691-5](https://link.springer.com/article/10.1007/s00285-013-0691-5). PubMed PMID: 23722629; PMCID: PMC3795870.
