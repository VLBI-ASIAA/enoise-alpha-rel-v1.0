/*****************************************************************************
*    <enoise: none-zero baseline VLBI data simulator>                        * 
*    Copyright (C) <2014> <Zheng Meyer-Zhao>                                 *
*                                                                            *
*    This program is free software: you can redistribute it and/or modify    *
*    it under the terms of the GNU General Public License as published by    *
*    the Free Software Foundation, either version 3 of the License, or       *
*    (at your option) any later version.                                     *
*                                                                            *
*    This program is distributed in the hope that it will be useful,         *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of          *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
*    GNU General Public License for more details.                            *
*                                                                            *
*    You should have received a copy of the GNU General Public License       *
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
*****************************************************************************/

/*
 *  cmdparser.c
 *
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "cmdparser.h"
#include "enoise.h"

int cmdline(int *argc, char ***argv, SetupInfo *setup)
{
  int x = options(*argc, *argv, setup);
  *argc -= optind;
  *argv += optind;
  return(x);
}

static int options(int argc, char **argv, SetupInfo *setup)
{
  int c;
  while ((c = getopt(argc, argv, "vb:d:n:r:t:q:")) != -1)
    switch(c) {
      case 'v': setup->verb++;                                           break;

      case 'b': set_bsmx(atoi(optarg), setup);                           break;
      case 'd': set_dura(atof(optarg), setup);                           break;
      case 'n': if ((c = set_comm(optarg, setup)))    return(c<0?-1:c);
                else                                                     break;
      case 'r': set_reps(atof(optarg), setup);                           break;
      case 't': set_thrs(atof(optarg), setup);                           break;
      default: return(-1);
    }
  return(0);
}

static void set_bsmx(int b, SetupInfo *setup)
{
  setup->bsmx = b;
  if (setup->verb>0) fprintf(stdout, "Bits states with %d counts\n", setup->bsmx);
}

static void set_dura(double d, SetupInfo *setup)
{
  setup->dura = d;
  if (setup->verb>0) fprintf(stdout, "Running for %g secs\n", setup->dura);
}

static void set_reps(double r, SetupInfo *setup)
{
  setup->nusr = (int)r;
  r = (double)setup->nusr / 1e6;
  if (setup->verb>0) fprintf(stdout, "Reporting every %g secs\n", r);
}

static void set_thrs(double t, SetupInfo *setup)
{
  setup->thrs = t;
  if (setup->verb>0) fprintf(stdout, "2-bit threshold %g sigma\n", setup->thrs);
}

static void set_corr(double c, SetupInfo *setup)
{
  setup->corr = c;
  if (setup->verb>0) fprintf(stdout, "Correlation coefficient now %g\n", setup->corr);
}

/* could add phase */
static int add_tone(char *frqamp, SetupInfo *setup)
{
  double  frq, amp, phs = 0;
  if ((3 != sscanf(frqamp, "%lf,%lf,%lf", &frq, &amp, &phs)) &&
      (2 != sscanf(frqamp, "%lf,%lf", &frq, &amp))) return(-1);
  if (setup->verb>0) fprintf(stdout,
      "Added tone %g at %g MHz (%g deg)\n", amp, frq, phs);
  setup->tphs[setup->ntones] = phs * M_PI / 180.0;
  setup->tfrq[setup->ntones] = frq;
  setup->tamp[setup->ntones] = amp;
  setup->ntones++;
}

static void set_fftn(int a, SetupInfo *setup)
{
  if (a < 0) a = 0;
  setup->fftsps = a;
  if (setup->verb>0) fprintf(stdout, "Requesting FFT %d points/us\n", setup->fftsps);
}

static void set_slices(int a, SetupInfo *setup)
{
  if (a < 1) a = 1;
  setup->slices = a;
  if (setup->verb>0) fprintf(stdout, "Requesting FFT every %d us\n", setup->slices);
}

static void set_pathd(char *path, SetupInfo *setup)
{
  strcpy(setup->pathd, path);
  if (setup->verb>0) fprintf(stdout, "Path to ndata directory is %s\n", setup->pathd); 
}

static void set_pathr(char *path, SetupInfo *setup)
{
  strcpy(setup->pathr, path);
  if (setup->verb>0) fprintf(stdout, "Path to drate directory is %s\n", setup->pathr);
}

static void set_skyfreq(double freq, SetupInfo *setup)
{
  setup->skyfreq = freq;
  if (setup->verb>0) fprintf(stdout, "Sky frequency is %f MHz\n", setup->skyfreq);
}

static void set_limit(double limit, SetupInfo *setup)
{
  setup->limit = limit;
  if (setup->verb > 0) fprintf(stdout, "Limit is %f\n", setup->limit);
}

static void set_spsmul(int spsmul, SetupInfo *setup)
{
  setup->spsmul = spsmul;
  if (setup->verb > 0) fprintf(stdout, "sps multiplier is %d\n", setup->spsmul);
}

static void set_stnoise(int stnoise, SetupInfo *setup)
{
  setup->stnoise = stnoise;
  if (setup->verb > 0) fprintf(stdout, "station noise ratio is %f\n", setup->stnoise);
}

static int set_comm(char *c, SetupInfo *setup)
{
  if (!strcmp(c, "help")) {
    printf( "  corr:<float>         amplitude of common signal relative\n"
            "                       to the station receiver noise (0.01)\n"
            "  tone:<freq,amp>      adds a tone at the <float> freq,amp\n"
            "  tone:<freq,amp,phs>  specifies phase as well\n"
            "  fftn:<int>           force underlying fft to have this many\n"
            "                       samples/us (0 -- means work it out)\n"
            "  slices:<int>         number of 1-us sample groups per fft (1)\n"
            "  pathd:<string>       path to simulation data, maximum 50 characters\n"
            "  pathr:<string>       path to delay and rate info, maximum 50 characters\n"
            "  skyfreq:<float>      sky frequency in MHz\n"
            "  limit:<float>        threshold to calculate delay error\n"
            "  spsmul:<int>         multiplier of sps (samples per microsecond)\n"
            "  stnoise:<float>      ratio of station noise [0..1]\n"
    );
    return(-1);
  } else if (!strncmp(c, "corr:", 5)) {
    if (atof(c+5) < 0 || atof(c+5) >= 1.0) {
      printf("corr value should be between 0.0 and 1.0\n");
      return(-1);
    }
    else
      set_corr(atof(c+5), setup);
  } else if (!strncmp(c, "tone:", 5)) {
    add_tone(c+5, setup);
  } else if (!strncmp(c, "fftn:", 5)) {
    set_fftn(atoi(c+5), setup);
  } else if (!strncmp(c, "slices:", 7)) {
    set_slices(atoi(c+7), setup);
  } else if (!strncmp(c, "pathd:", 6)) {
    set_pathd(c+6, setup);
  } else if (!strncmp(c, "pathr:", 6)) {
    set_pathr(c+6, setup);
  } else if (!strncmp(c, "skyfreq:", 8)) {
    set_skyfreq(atof(c+8), setup);
  } else if (!strncmp(c, "limit:", 6)) {
    set_limit(atof(c+6), setup);
  } else if (!strncmp(c, "spsmul:", 7)) {
    set_spsmul(atoi(c+7), setup);
  } else if (!strncmp(c, "stnoise:", 8)) {
    set_stnoise(atoi(c+8), setup);  
  } else {
    fprintf(stderr, "What is %s ?\n", c);
    return(-1);
  }
  return(0);
}

void set_csigma(SetupInfo *setup)
{
  setup->csigma = sqrt (setup->corr / (1.0 - setup->corr));
  if (setup->verb>0) fprintf(stdout, "Csigma is %g\n", setup->csigma);
}

/*
 * eof
 */
