/* 
 * simulation of sequence-by-synthesis with cycle numer fixed, and 
 * compare the simulation with theoretical results
 * 
 * The "I" in the name refers to the def I of book 4, p27
 *
 * Yong Kong
 * 11/7/2008
 *
 * $Id: sbs_ffcm_r.c,v 1.2 2012/08/31 21:36:09 yk336 Exp $
 */

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/time.h>
#include <ctype.h>
#include <time.h>
#include <dirent.h>

#include <assert.h>


#define MAX_RJ      20            /* max of rj */
#define MAX_C       256           /* upper limit of incorporation cycle */
#define MAX_N       1024          /* max of sequence length */


long long baseCnt[4];             /* base cnt */
double nSum = 0.0, nSum2 = 0.0;   /* seq length */
double rSum = 0.0, rSum2 = 0.0;   /* cond seq length */
double rjSum[MAX_RJ];             /* how many runs with length 1, 2, etc */
double rjSum2[MAX_RJ];

double ndist[MAX_N];              /* length counts */
double ndist4[MAX_N][4];          /* length counts of the 4 nt flows */

double rdist[MAX_N];              /* r length counts */
double rdist4[MAX_N][4];          /* r length counts of the 4 nt flows */

double dist[MAX_RJ][MAX_N];       
    /* dist[i][j]: count of runs of length i appearing j times */



void usage(const char *name) {
  fprintf(stderr, "%s: -a <prob of A> -b <prob of C> -c <prob of G> -d <prob of T> -A a0,a1,... -B b0,b1,...  -C c0,c1,... -D d0,d1,...\n", name);
  exit(-1);
}

void seq_gen(int cycle, double *cumProb, double (*cumPr)[MAX_C], int *pr_cnt) {
  int f = 1;   /* cycle count */
  int n = 0;   /* sequence length */
  int r = 0;   /* cond sequence length */
  int i = 0, j = 0, ii, ok, lastj = 0;
  int rj = 0;  /* track the number of bases at one cycle (not quad cycle): 
		  the length of one run */
  //  int lastbase = (unsigned) (drand48() * 4);
  double rand;
  int prev_i = -1, prev_f = -1;
  int rjOne[MAX_RJ];  /* for one sequence, how many 1's, 2's, etc.:
		       rjOne[0]: total 0's;
		       rjOne[1]: total 1's;
		       rjOne[2]: total 2's;
		      */

  for (i=0; i<MAX_RJ; i++)
    rjOne[i] = 0;


  i = 0;  /* flow; it not a good idea to use i, j, k to do real things! */
  while (1) { /* cycle */
    /* create a base of the template */
    ok = 0;
    rand = drand48();
    //    lastj = j;  /* the last base that ends up in cycle f */
    for (j = 0; j < 4; j++) { /* j: base of the template */
      if (rand < cumProb[j]) {
	ok = 1;
	baseCnt[j]++;
	break;
      }
    }
    // printf("J: %d %d %d\n", j, i, f);
    //printf("[%1d]", j);
    /* now j is the base to be incorporated */

    if (ok == 0) {
      printf("The base is not incorporated! (rand = %f)\n", rand);
      exit(-1);
    } 

    /* 
     * flush the flows, until the base is incorporated, or the cycle
     *  limit is reached 
     */
    while (f <= cycle) {
      // printf("<%1d:%d>", i, f);
      if (i != j) { /* flow and base do not match */
	i++;        /* i: flow */
	
	if (rj < MAX_RJ)  /* reset rj */
          rjOne[rj]++;	
        rj = 0;
	if (i > 3) {
	  i = 0;
	  f++;
	  if (f > cycle) {
	    break;
	  }
	}
	// printf("I: %d %d %d\n", i, j, f);
      } else { /* incorporate the base based on probs */
	ok = 0;
	rand = drand48();
	// printf("II: %d %d %d %f\n", i, j, f, rand);
	/* ii: index of gi  */
	for (ii = 0; ii < pr_cnt[i]; ii++) {
	  if (rand < cumPr[i][ii]) { /* the base will be incorporated */ 
	    ok = 1;
	    break;
	  } else {
	    f++;
	    // printf("H %d %d %d %f %f\n", i, ii, f, cumPr[i][ii], rand); 
	    if (f > cycle) {
	      break;
	    }
	  }
	}

	//	printf("[%1d;%1d:%1d:%1d]", j,f,i,ii);

	if (ok == 0 && f <= cycle) {
	  printf("The base is not incorporated! (f: %d n: %d r: %d base: %d flow: %d g_cnt: %d ii: %d cumPr: %f rand: %f)\n", 
		 f, n, r, j, i, pr_cnt[i], ii, cumPr[i][ii], rand);
	  exit(-1);
	}      
	 
	if (f <= cycle) {
	  /* base got incorporated */
	  //	  printf("%1d(%1d:%d) ", j, rj, f);
	  rj++;
	  n++;    /* seq length */
	  lastj = j;

	  if ((i != prev_i) || 
	      (i ==  prev_i && f != prev_f)) {
	    r++;

	    prev_f = f;
	    prev_i = i;
	  }
	  break;
	}
      }
    } /* while (f <= cycle) */
    // printf("F %d\n", f);
    if (f > cycle) {
      break;
    }
  }

  //printf(" (%1d:%1d) :: %d %d %d (%d %d %d)\n", 
  //	 lastj, j, 
  //	 n, r, f, 
  //	 rjOne[0], rjOne[1], rjOne[2]);  /* j: the last one out of cycle */
  nSum += n;
  nSum2 += (n * n);
  rSum += r;
  rSum2 += (r * r);

  if (n < MAX_N) {
    ndist[n]++;
    ndist4[n][lastj]++;
  }

  if (r < MAX_N) {
    rdist[r]++;
    rdist4[r][lastj]++;
  }

  for (i=0; i<MAX_RJ; i++) {
    rjSum[i] += rjOne[i];
    rjSum2[i] += rjOne[i] * rjOne[i];

    if (rjOne[i] < MAX_N) {
      dist[i][rjOne[i]]++;
    }
  }

}




/* s2 */
double get_s2(double *prob) {
  double s23, s2;

  s23 = prob[2] + prob[3];
  s2 = prob[0] * (prob[1] + s23)
    +  prob[1] * s23
    +  prob[2] * prob[3];

  return s2;
}

/* s3 */
double get_s3(double *prob) {
  double s3;

  s3 = prob[0] * prob[1] * prob[2]
    +  prob[0] * prob[1] * prob[3]
    +  prob[0] * prob[2] * prob[3]
    +  prob[1] * prob[2] * prob[3];

  return s3;
}

/* s4 */
double get_s4(double *prob) {
  double s4;

  s4 = prob[0] * prob[1] * prob[2] * prob[3];

  return s4;
}

/* rfac */
double get_rfac(double *prob) {
  double rfac;

  rfac = prob[0]*prob[0]*(1.0 - prob[0]) +
    prob[1]*prob[1]*(prob[2] + prob[3]) +
    prob[2]*prob[2]*(prob[3]);
    
  return rfac;
}


/* t1'(1) and "var" of t1 */
double get_t1p(double *prob, double (*pr)[MAX_C], int *pr_cnt, 
	       double *var_t1) {
  int i, j;
  double g, g2, g2sum, t1p;
  
  t1p = 0.0;
  g2sum = 0.0;
  for (i = 0; i < 4; i++) {
    g = 0.0;
    g2 = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g += j * pr[i][j];
      g2 += j*(j-1) * pr[i][j];
    }
    g *= prob[i];
    g2 *= prob[i];

    t1p += g;
    g2sum += g2;
  }

  *var_t1 = g2sum + t1p - t1p*t1p;

  return t1p;
}

/* t1''(1)  */
double get_t1p2(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g, t1p2;
  
  t1p2 = 0.0;
  for (i = 0; i < 4; i++) {
    g = 0.0;
    for (j = 2; j < pr_cnt[i]; j++) {
      g += j * (j-1) * pr[i][j];
    }
    g *= prob[i];

    t1p2 += g;
  }

  return t1p2;
}


/* t1'''(1)  */
double get_t1p3(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g, t1p3;
  
  t1p3 = 0.0;
  for (i = 0; i < 4; i++) {
    g = 0.0;
    for (j = 3; j < pr_cnt[i]; j++) {
      g += j * (j-1) * (j-2) * pr[i][j];
    }
    g *= prob[i];

    t1p3 += g;
  }

  return t1p3;
}


/* t2'(1) -- see yale notebook 3 p12 */
double get_t2p(double *prob, double (*pr)[MAX_C], int *pr_cnt, double t1p) {
  int i, j;
  double g[4], t, t2p;
  
  for (i = 0; i < 4; i++) {
    g[i] = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g[i] += j * pr[i][j];
    }
    g[i] *= prob[i]; /* g_i'(1) */
  }

  t = 0.0;
  for (i = 0; i < 4; i++) 
    t += prob[i] * g[i];

  t2p = t1p - t;

  return t2p;
}

/* t2''(1) */
double get_t2p2(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g1[4], g2[4], t2p2;
  
  for (i = 0; i < 4; i++) {
    g1[i] = 0.0;
    g2[i] = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g1[i] += j * pr[i][j];
      g2[i] += j * (j-1) * pr[i][j];
    }
    g1[i] *= prob[i]; /* g_i'(1) */
    g2[i] *= prob[i]; /* g_i''(1) */
  }

  t2p2 = 0.0;
  for (i = 0; i < 4; i++) {
    t2p2 += (g2[i] * (1.0 - prob[i]));
  }

  t2p2 += 2.0 * (g1[0] * g1[1] + g1[0] * g1[2] + g1[0] * g1[3]
		 + g1[1] * g1[2] + g1[1] * g1[3]
		 + g1[2] * g1[3]);

  return t2p2;
}


/* t3'(1) */
double get_t3p(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g[4], t3p;
  
  for (i = 0; i < 4; i++) {
    g[i] = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g[i] += j * pr[i][j];
    }
    g[i] *= prob[i]; /* g_i'(1) */
  }

  t3p = g[0] * (prob[1] * prob[2] + prob[1] * prob[3] + prob[2] * prob[3]) +
        g[1] * (prob[0] * prob[2] + prob[0] * prob[3] + prob[2] * prob[3]) +
        g[2] * (prob[0] * prob[1] + prob[0] * prob[3] + prob[1] * prob[3]) +
        g[3] * (prob[0] * prob[1] + prob[0] * prob[2] + prob[1] * prob[2]);

  return t3p;
}


double pi (double *prob, int p) {
  return pow(prob[0], (double)p) + pow(prob[1], (double)p)
    + pow(prob[2], (double)p) + pow(prob[3], (double)p);
}

/* 
 * this one is before I found the theoretical results;
 * it uses the avg of n to calculate avg of count of runs of length j;
 * remarkably, the avg is not that far;
 * the var, on the other hand, seems more than doubled;
 */
double rj_avg_var(double *prob, int n, int j, double *var) {
  double rv;
  double pr0, pr1, pr2, pr20, pr21, pr22, pr23;
  double vA, vB, vC, vD, vE, vF1, vF2;

  if (j == n) {
    rv = pi(prob, n);
    *var = 0.0;  /* not exact right */
  } else {
    pr0 = pi(prob, j);
    pr1 = pi(prob, j + 1);
    pr2 = pi(prob, j + 2);

    pr20 = pi(prob, 2*j);
    pr21 = pi(prob, 2*j+1);
    pr22 = pi(prob, 2*j+2);
    pr23 = pi(prob, 2*j+3);

    vA = pr0 - 2.0 * pr1 + pr2;
    vB = 1 - ((2*j - 1)*pr0 - (4*j+2)*pr1 + (2*j+3)*pr2);
    vC = -2*pr20 + 6*pr21 - 6*pr22 + 2*pr23;

    vD = (j-1.0)*pr0 - 2.0*j*pr1 + (j+1.0)*pr2;
    vE = (4*j-2)*pr20 -12*j*pr21 + (12*j+6)*pr22 - 4*(j+1)*pr23;
    vF1 = (j-1)*(3*j-1)*pr0*pr0 + 4*j*(3*j+2)*pr1*pr1 + (j+1)*(3*j+5)*pr2*pr2;
    vF2 = (-4*j)*(3*j-1)*pr0*pr1 + (6*j*j+4*j+2)*pr0*pr2 - 4*(j+1)*(3*j+2)*pr1*pr2;
    
    rv = vA * n - vD;
    *var = (vA*vB + vC)*n - vD + vE + vF1 + vF2;
  } 
  
  return rv;
}

double rj_avg_var_new(double *prob, int f, int j, double *var) {
  double rv;
  double pr0, pr1, pr2, pr3, pr4, e2, e3, e4, e22, e23, e24, e32;
  double pr2_0, pr2_1, pr2_2, pr2_3, pr2_4;
  double vA, vB, vC, vD1, vD2, vD, vE, vF, vA2;
  double var1, var0;

  pr0 = pi(prob, j);
  pr1 = pi(prob, j + 1);
  pr2 = pi(prob, j + 2);
  pr3 = pi(prob, j + 3);
  pr4 = pi(prob, j + 4);

  pr2_0 = pi(prob, 2*j);
  pr2_1 = pi(prob, 2*j + 1);
  pr2_2 = pi(prob, 2*j + 2);
  pr2_3 = pi(prob, 2*j + 3);
  pr2_4 = pi(prob, 2*j + 4);


  e2 = get_s2(prob);
  e3 = get_s3(prob);
  e4 = get_s4(prob);

  e22 = e2 * e2;
  e23 = e22 * e2;
  e24 = e23 * e2;

  e32 = e3 * e3;

  vA = pr0 - 2.0 * pr1 + pr2;
  vB = pr0 - 3.0 * pr1 + 3.0 * pr2 - pr3;
  vC = pr0 - 4.0 * pr1 + 5.0 * pr2 - 2.0*pr3;

  vD1 = pr1 - 2.0 * pr2 + pr3;
  vD2 = pr0 - 8.0 * pr1 + 17.0 * pr2 - 14.0*pr3 + 4.0 * pr4;
  vD = vD1*vD1 + vA * vD2;

  vE = pr2_0 - 3.0 * pr2_1 + 3.0 * pr2_2 - pr2_3;
  vF = pr2_0 - 4.0 * pr2_1 + 6.0 * pr2_2 - 4.0 * pr2_3 + pr2_4;
    
  rv = vA * f/e2 - (e2 * vB - 2.0 * e3 * vA)/(e22);

  vA2 = vA * vA;
  var1 = (2.0*vA2*e3 + (vA - 2.0*vE)*e22 - vA*vC*e2) / e23;
  var0 = (vD*e22 - 6.0*vA2*e2*e4 + 8.0*vA2*e32 + (2.0*vA - 4.0*vE)*e22*e3 - 8.0*vA*vB*e2*e3 + (3.0*vA2 - vB + 3.0*vF)*e23) / e24;
  *var = var1 * f + var0;
  
  return rv;
}


double get_t01 (double *prob, double (*pr)[MAX_C]) {
  int i;
  double rv = 0.0;

  for (i=0; i<4; i++) {
    rv += prob[i]*pr[i][0];
  }
  return rv;
}

double P2 (double *prob, double (*pr)[MAX_C]) {
  int i;
  double rv = 0.0;

  for (i=0; i<4; i++) {
    rv += prob[i]*prob[i]*pr[i][0];
  }
  return rv;
}

double P3 (double *prob, double (*pr)[MAX_C]) {
  int i;
  double rv = 0.0;

  for (i=0; i<4; i++) {
    rv += prob[i]*prob[i]*prob[i]*pr[i][0];
  }
  return rv;
}


double get_h11p (double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g, h11p, pp1 = 0.0;
  
  /* e1 of m0 */
  for (i = 0; i < 4; i++) {
    pp1 += prob[i]*pr[i][0];
  }

  h11p = 0.0;
  for (i = 0; i < 4; i++) {
    g = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g += j * pr[i][j];
    }
    g *= prob[i]*(pp1 - prob[i]*pr[i][0]);

    h11p += g;
  }

  return h11p;
}

/*
 * yale book 4, p14
 */
double r_avg(double *prob, double (*pr)[MAX_C], int *pr_cnt, int f) {
  double pp2, pp3, e2, e3, t1p, t1p2, t2p, t01, h11p, var_t1;
  double u, v, avg, term1, term2;
  
  e2 = get_s2(prob);
  e3 = get_s3(prob);

  t01 = get_t01(prob, pr);
  pp2 = P2(prob, pr);
  pp3 = P3(prob, pr);
  h11p = get_h11p(prob, pr, pr_cnt);
  //  printf("h11p: %f\n", h11p);


  t1p  = get_t1p(prob, pr, pr_cnt, &var_t1);
  t1p2 = get_t1p2(prob, pr, pr_cnt);
  t2p  = get_t2p(prob, pr, pr_cnt, t1p);

  u = 1.0 / ( 1.0 - pp2 );
  v = 1.0 / (e2 + t1p);

  //term1 = (-2.0*e2 + pp2 - pp3) / v;
  //term2 = (2.0*e3) / u;

  term1 = (t1p*t01 -t1p -h11p -2.0*e2 + pp2 - pp3) / v;
  term2 = (t1p2 + 2.0*t2p + 2.0*e3) / u;

  avg = v / u * f + v*v*(term1 + term2);

  return avg;
}

int main(int argc, char *argv[]) { 
  int i, j, k;
  int cycle = 10, repeat = 1000, ch, rsum = 0, peak, zerocnt;
  int fixedSeed = 0, pad = 0;
  double prob[] = {0.25, 0.25, 0.25, 0.25};  /* prob of 4 nucleotides */
  double cumProb[4];                         /* cumulative probs */
  double pr[4][MAX_C];                       /* incorporation probs */
  double cumPr[4][MAX_C];

  double sum = 0.0;
  char *t;
  int pr_cnt[4];                  /* upper index of incorporation probs */

  double nAvg, nVar;
  double rAvg, rVar;
  double nAvgth, nVarth;
  double rAvgth, rVarth;
  double nfac, nfac2, nfac3, nfac4, vfac, tmp2, u, v;
  double s2, s3, s4, t1p, var_t1, t2p, t1p2, t1p3, t2p2, t3p;
  double rjAvg, rjVar;
  double rjAvgth, rjVarth;
  double rjAvgthNew, rjVarthNew;
  double sumL[MAX_RJ], sumL2[MAX_RJ];
  double ffc_factor, mean_ffc, var_ffc;

  char fn[1024];
  FILE *fp;

  while ((ch = 
	  getopt(argc, argv, 
		 "a:b:c:d:sf:r:p:A:B:C:D:"
		 )) != -1)
    switch (ch) {
    case 'a':
      prob[0] = atof(optarg);             break;
    case 'b':
      prob[1] = atof(optarg);             break;      
    case 'c':
      prob[2] = atof(optarg);             break;
    case 'd':
      prob[3] = atof(optarg);             break;      
    case 's':
      fixedSeed = 1;                      break;   
    case 'p':
      pad       = 1;                      break;   
    case 'f':
      cycle = atoi(optarg);              break;   
    case 'r':
      repeat = atoi(optarg);              break;   
      

    case 'A':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[0][i] = atof(t);
      
      pr_cnt[0] = i;
      break;   

    case 'B':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[1][i] = atof(t);
 
      pr_cnt[1] = i;
      break;   

    case 'C':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[2][i] = atof(t);
 
      pr_cnt[2] = i;
      break;   

    case 'D':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[3][i] = atof(t);

      pr_cnt[3] = i;
      break;   


    default:
      usage(argv[0]);
    }
  
  /* make sure the probs add up to 1 */
  for (i = 0; i < 4; i++) {
    if (prob[i] < 0.0 || prob[i] > 1.0) {
      usage(argv[0]);
    }
    sum += prob[i];
  }

  if (fabs(sum - 1.0) > 1.0e-10) {
    printf("prob sum: %.12f\n", sum);
    usage(argv[0]);
  }

  for (i = 0; i < 4; i++) {
    sum = 0.0;
    for (j = 0; j < pr_cnt[i]; j++) {

      if (pr[i][j] < 0.0 || pr[i][j] > 1.0) {
	usage(argv[0]);
      }
      
      sum += pr[i][j];
    }
    if (fabs(sum - 1.0) > 1.0e-10) {
      printf("%c\t%f\n", i + 65, sum);
      usage(argv[0]);
    }
  }

  /* set up a unique file name */
  sprintf(fn, "sbs_ffcm_r_%d_repeat_%d", cycle, repeat);
  for (i=0; i<4; i++) {
    sprintf(fn, "%s_%c%4.3f", fn, i+65, prob[i]);
    for (j=0; j<pr_cnt[i]; j++)
      sprintf(fn, "%s_%3.2f", fn, pr[i][j]);
  }
  printf("%s\n", fn);
  fp = fopen(fn, "w");
  if (fp == NULL) {
    fprintf(stderr, "cannot open file %s to write!\n", fn);
    exit(-1);
  }


  fprintf(fp, "######## Input parameters ###########################################\n");
  if (fixedSeed) {
    srand48(12345);
    fprintf(fp, "#\tseed: fixed as 12345\n");
  } else {
    srand48(time(NULL));
    fprintf(fp, "#\tseed: use time\n");
  }

  fprintf(fp, "#\ttflowcycle: %d\n", cycle);
  fprintf(fp, "#\trepeat: %d\n", repeat);

  /* set up cumulative probs */
  cumProb[0] = prob[0];
  fprintf(fp, "\n#\tnucleotide composition probabilities\tcumulative\n");  
  fprintf(fp, "#\t%d\t%f\t%f\n", 0, prob[0], cumProb[0]);
  for (i = 1; i < 4; i++) {
    cumProb[i] = cumProb[i-1] + prob[i];
    fprintf(fp, "#\t%c\t%f\t%f\n", i+97, prob[i], cumProb[i]);
  }

  fprintf(fp, "\n#\tnucleotide incorporation probabilities\tcumulative\n"); 
  for (i = 0; i < 4; i++) {
    cumPr[i][0] = pr[i][0];
    fprintf(fp, "#\t%c\t%d\t%f\t%f\n", i+97, 0, pr[i][0], cumPr[i][0]);
    for (j = 1; j < pr_cnt[i]; j++) {
      cumPr[i][j] = cumPr[i][j-1] + pr[i][j];
      fprintf(fp, "#\t%c\t%d\t%f\t%f\n", i+97, j, pr[i][j], cumPr[i][j]);
   }
    fprintf(fp, "\n");
  }

  for (i = 0; i < 4; i++) {
    baseCnt[i] = 0;
  }

  for (i = 0; i < MAX_RJ; i++) {
    rjSum[i] = 0.0; 
    rjSum2[i] = 0.0; 
  }

  for (i = 0; i < MAX_RJ; i++) 
    for (j = 0; j < MAX_N; j++) 
      dist[i][j] = 0.0; 
  
  for (j = 0; j < MAX_N; j++) {
    ndist[j] = 0.0;
    rdist[j] = 0.0;
    for (i = 0; i < 4; i++)
      ndist4[j][i] = 0.0; 
      rdist4[j][i] = 0.0; 
  }



  /*** run the simulation ***/
  for (k = 0; k < repeat; k++) {
    seq_gen(cycle, cumProb, cumPr, pr_cnt);
  }


  fprintf(fp, "#####################################################################\n");
  fprintf(fp, "######## Theoretical results of sequence length distribution ########\n");


  s2 = get_s2(prob);
  s3 = get_s3(prob);
  s4 = prob[0] * prob[1] * prob[2] * prob[3]; 
  t1p  = get_t1p(prob, pr, pr_cnt, &var_t1);
  t2p  = get_t2p(prob, pr, pr_cnt, t1p);
  t1p2 = get_t1p2(prob, pr, pr_cnt);
  t1p3 = get_t1p3(prob, pr, pr_cnt);
  t2p2 = get_t2p2(prob, pr, pr_cnt);
  t3p  = get_t3p(prob, pr, pr_cnt);

  nfac = s2 + t1p;
  vfac = 2.0*s3 + s2 - 3.0*s2*s2 + t1p2 + t1p - t1p*t1p + 2.0*t2p 
    - 4.0*t1p*s2;  /* w in jcb p824 */

  nfac2 = nfac  * nfac;
  nfac3 = nfac2 * nfac;
  nfac4 = nfac3 * nfac;

  // 4/15/2011
  u = nfac;
  v = 2*t2p + t1p2 + 2*s3;

  tmp2 = v*v;


  /* 4/25/2010 */
  nAvgth = 1.0/nfac * cycle 
    + 0.5/(nfac2) * (-2.0*s2*s2 + 2.0*s3 
		     - 4.0*t1p*s2 - 2.0*t1p*t1p + t1p2 + 2.0*t2p); 
  
  // 4/15/2011 a newer version
  nVarth = 1.0/(nfac3) * (v - (3*s2 + t1p -1)*u) * cycle 
    + 1.0/12.0/(nfac4) * (
			  6.0*u*v*(3*u-4*s2) + 15.0*v*v
			  -8.0*u*(
				  6.0*t3p + t1p3 + 3.0*t2p2 + 6.0*s4
				  + 3.0*u*(t1p2 + t2p)
				  )
			  );



  fprintf(fp, "\tAvg\t%f\tVar\t%f\n", 
	  nAvgth, 
	  nVarth);

  if (nAvgth >= 0.5*MAX_N) {
    fprintf(stderr, "MAX_N should be increased!\n");
    exit(-1);
  }

  fprintf(fp, "######## Simulation results of sequence length distribution #########\n");
  /* avg and var from simulation */
  nAvg = nSum / (double)repeat;
  nVar = nSum2 / (double)repeat -  nAvg * nAvg;
  fprintf(fp, "\tAvg\t%f\tVar\t%f\n", 
	  nAvg, 
	  nVar);

  fprintf(fp, "#####################################################################\n");


  fprintf(fp, "\n\n");
  fprintf(fp, "#####################################################################\n");
  fprintf(fp, "############### Theoretical results of r distribution ###############\n");

  rAvgth = r_avg(prob, pr, pr_cnt, cycle); /* incomplete incorporation */
  rVarth = 2.0*s3/s2 * cycle + (2.0*s2*s2*s2 - 3*s2*s3 + 5*s3*s3 - 4*s2*s4) 
    /(s2 * s2);  // do we have a var for incomplete incorporation ??? 2012/8/22



  // 2012/8/31 this factor appears in the pyro version of mean and var of r(f)
  ffc_factor = prob[1]*prob[2]*prob[3]
    + prob[1]*prob[1]*(prob[2]+prob[3])
    + prob[2]*prob[2]*prob[3];

  // 2012/8/31 (see pyro_run_relation.mw)
  mean_ffc = 2.0 * cycle + (s3-s2 + prob[0]*s2 + ffc_factor)/s2;
  var_ffc  = 2.0 * s3 * cycle /s2 + (3.0*s3*s3 -s2*s3 -2.0*s2*s4   -prob[0]*(-1.0+prob[0])*s2*s2 -  ffc_factor*ffc_factor + (2.0*prob[1]*prob[3]*prob[2]*prob[2] + ffc_factor)*s2)/(s2*s2);
  

  //printf("ffc_factor: %.10f  e3: %.10f\n", ffc_factor, s3);

  fprintf(fp, "\tAvg\t%f\tVar\t%f\n", rAvgth, rVarth);
  fprintf(fp, "\tAvg\t%f\tVar\t%f (FFCM)\n", mean_ffc, var_ffc);

  fprintf(fp, "############### Simulation results of r distribution ################\n");
  rAvg = rSum / (double)repeat;
  rVar = rSum2 / (double)repeat -  rAvg * rAvg;
  fprintf(fp, "\tAvg\t%f\tVar\t%f\n", rAvg, rVar);
  fprintf(fp, "#####################################################################\n");




  fprintf(fp, "\n");
  fprintf(fp, "######### base counts \n");
  for (i = 0; i < 4; i++) {
    fprintf(fp, "#\t%c\t%lld\n", i+97, baseCnt[i]);
  }

  /* use a simple peak finder to truncate the unneccessary printout */
  fprintf(fp, "\n");
  fprintf(fp, "######### length distribution########################################\n");
  fprintf(fp, "n\tabcd_count\tabcd_frac\ta_count\ta_frac\tb_count\tb_frac\tc_count\tc_frac\td_count\td_frac\n");
  peak = 0;
  zerocnt = 0;
  for (j=0; j<MAX_N; j++) {
    fprintf(fp, "%d\t%f\t%f\t", j, ndist[j], ndist[j]/repeat);
    for (i = 0; i < 4; i++)
      fprintf(fp, "%f\t%f\t", ndist4[j][i], ndist4[j][i]/repeat);
    fprintf(fp, "\n");
    rsum += ndist[j];

    if (ndist[j] > 0) 
      peak = 1;

    if (peak && ndist[j] == 0.0)
      zerocnt++;

    if (zerocnt > 100)
      break;
  }




  
  fprintf(fp, "#\n");
  fprintf(fp, 
	  "#\tj\trjAvg\trjVar\trjAvgth\trjVarth\trjAvgthNew\trjVarthNew\n");
  for (i = 0; i < MAX_RJ; i++) {
    rjAvg = rjSum[i]/(double)repeat;
    rjVar = rjSum2[i]/(double)repeat - rjAvg * rjAvg;

    rjAvgth = rj_avg_var(prob, nAvgth, i, &rjVarth);

    rjAvgthNew = rj_avg_var_new(prob, cycle, i, &rjVarthNew);

    fprintf(fp, "#\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", 
    i, rjAvg, rjVar, rjAvgth, rjVarth, rjAvgthNew, rjVarthNew); 
  }

  for (i=0; i<MAX_RJ; i++) {
    //    fprintf(fp, "\n\n");
    sumL[i] = sumL2[i] = 0.0;
    for (j=0; j<MAX_N; j++) {
      //  fprintf(fp, "%d\t%d\t%f\n", i, j, dist[i][j]);
      sumL[i]  += j*  dist[i][j];
      sumL2[i] += j*j*dist[i][j];
    }
    
    sumL[i]  = sumL[i] / (double)repeat;
    sumL2[i] = sumL2[i] /(double)repeat - sumL[i]*sumL[i];
    //fprintf(fp, "\n\n%f\t%f\n", sumL[i], sumL2[i]);
  }

  fprintf(fp, "#\tthe next line: (avg and var) from simulation\n");
  fprintf(fp, "#");
  for (i=0; i<MAX_RJ; i++) {
    fprintf(fp, "(%f %f)\t", sumL[i], sumL2[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "Total sequences: %d\n", rsum);


  fprintf(fp, "\n");
  fprintf(fp, "######### r distribution########################################\n");
  fprintf(fp, "n\tabcd_count\tabcd_frac\ta_count\ta_frac\tb_count\tb_frac\tc_count\tc_frac\td_count\td_frac\n");
  rsum = 0;
  peak = 0;
  zerocnt = 0;
  for (j=0; j<MAX_N; j++) {
    fprintf(fp, "%d\t%f\t%f\t", j, rdist[j], rdist[j]/repeat);
    for (i = 0; i < 4; i++)
      fprintf(fp, "%f\t%f\t", rdist4[j][i], rdist4[j][i]/repeat);
    fprintf(fp, "\n");
    rsum += rdist[j];

    if (rdist[j] > 0) 
      peak = 1;

    if (peak && rdist[j] == 0.0)
      zerocnt++;

    if (zerocnt > 10)
      break;
  }
  fprintf(fp, "Total sequences: %d\n", rsum);

  fprintf(fp, "# rj distribution\n");
  for (j=0; j<MAX_N; j++) {
    fprintf(fp, "%d\t", j);
    for (i=0; i<MAX_RJ; i++) {
      fprintf(fp, "%f\t", dist[i][j]);
    }
    fprintf(fp, "\n");
  }


  fclose(fp);
  exit (0);
}
