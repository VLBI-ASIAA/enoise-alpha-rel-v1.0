/*****************************************************************************
 * *    <enoise: none-zero baseline VLBI data simulator>                        * 
 * *    Copyright (C) <2013> <Zheng Meyer-Zhao>                                 *
 * *                                                                            *
 * *    This program is free software: you can redistribute it and/or modify    *
 * *    it under the terms of the GNU General Public License as published by    *
 * *    the Free Software Foundation, either version 3 of the License, or       *
 * *    (at your option) any later version.                                     *
 * *                                                                            *
 * *    This program is distributed in the hope that it will be useful,         *
 * *    but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * *    GNU General Public License for more details.                            *
 * *                                                                            *
 * *    You should have received a copy of the GNU General Public License       *
 * *    along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 * *****************************************************************************/

/*
 *  random/genran.c
 *
 *  Author: Zheng Meyer-Zhao
 *  2013/05/24
 *
 *  generate random numbers with standard deviation 1.0
 *  users can define the random generator seed using option -s
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <gsl/gsl_randist.h>

#define BLOCK 10000
#define MAX_NAME_LEN 20
#define MAX_NUM_FILES 10

/* user controlled */
unsigned long int seed = 25786;
int honest = 4;
int nfiles = 0;

char names[MAX_NUM_FILES][MAX_NAME_LEN];
double durs[MAX_NUM_FILES];
/* random number generator data */
gsl_rng         *rng_inst;
const gsl_rng_type * T;

static int usage(char *name)
{
  printf("Usage: %s [options]\n\n", name);
  printf("where the options are:\n\n"
      "   -s <float>    seed to use\n"
      "   -h <int>      quality of random number 1 - 5, where 5 is the best\n"
      "   -n <string>   file name and duration, e.g. common.dat:4.4\n"
  );
  return(0);
}

void set_seed(unsigned long int aseed)
{
  seed = aseed;
}

void set_honest(int ahonest)
{
  honest = ahonest;
}

void set_datainfo(char *datainfo)
{
  char name[MAX_NAME_LEN-1];
  double dur;

  char *ptrs;
  int numStr = 0;
  char *temp[2];

  assert(nfiles < MAX_NUM_FILES);

  ptrs = strtok(datainfo, ":");
  while (ptrs != NULL && numStr < 2)

  {
    temp[numStr] = ptrs;
    numStr++;
    ptrs = strtok(NULL, " ");
  }

  strcpy(name, temp[0]);
  dur = atof(temp[1]);
  printf("%s, %.6f\n", name, dur);

  strcpy(names[nfiles], name);
  durs[nfiles] = dur;
  nfiles++;
}

static int parsecmd(int argc, char **argv)
{
  int c;
  while ((c = getopt(argc, argv, "s:h:n:")) != -1)
    switch(c) {
      case 's': set_seed(atof(optarg));   break;
      case 'h': set_honest(atoi(optarg)); break;
      case 'n': set_datainfo(optarg);     break;
      default :                           return(1);
    }
  return(0);
}

static void common_gaussian_setup()
{
  gsl_rng_env_setup();
  switch (honest)
  {
    case 1:  rng_inst = gsl_rng_alloc(gsl_rng_rand);    break;
    case 2:  rng_inst = gsl_rng_alloc(gsl_rng_taus2);   break;
    case 3:  rng_inst = gsl_rng_alloc(gsl_rng_mt19937); break;
    case 4:  rng_inst = gsl_rng_alloc(gsl_rng_ranlux);  break;
    case 5:  rng_inst = gsl_rng_alloc(gsl_rng_ranlux389);  break;
    default:         rng_inst = gsl_rng_alloc(gsl_rng_default); break;
  }
  gsl_rng_set(rng_inst, seed);
  printf("GSL type: %s (honest=%d), ", gsl_rng_name (rng_inst), honest);
  printf("seed=%u: ", seed);
  printf("%u, ", gsl_rng_get (rng_inst));
  printf("%u ...\n", gsl_rng_get (rng_inst));
}

void gen_random(int ii)
{
  printf("Generating random data ...\n");

  int us, ss;
  float *cd = (float *)malloc(BLOCK * sizeof(float));
  float *ptr = cd;

  FILE *fp;
  fp = fopen(names[ii], "w");
  
  printf("durs is %.6f\n",durs[ii]);
  for (us = 0; us < durs[ii] * 1000000; us++)
  {
    for (ss = 0; ss < BLOCK; ss++, cd++)
    {
      *cd = gsl_ran_gaussian_ziggurat(rng_inst, 1.0);
    }
    cd = ptr;
    fwrite(cd, sizeof(float), BLOCK, fp);
  }

  free(cd);

  fclose(fp);
  printf("Data writen to %s ...\n", names[ii]);
  printf("size of float is %d\n", sizeof(float));

}

int main(int argc, char **argv)
{
  int ii, err = 0;
  if (err = parsecmd(argc, argv) || argc == 1) 
  {
    return(usage(argv[0]));
  }
  common_gaussian_setup();

  printf("nfiles is %d\n", nfiles);
  for(ii = 0; ii < nfiles; ii++)
    gen_random(ii);

  gsl_rng_free(rng_inst);
  return(0);
}
