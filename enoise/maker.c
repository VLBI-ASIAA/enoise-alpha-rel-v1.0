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
 *  maker.c
 *  
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 *
 *  This file contains functions to simulate the data.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "maker.h"
#include "enoise.h"
#include "station.h"
#include "stdata.h"
#include "tones.h"

/*
 *  generate fake data
 */
int fake_the_data(SetupInfo *setup)
{
  int ns;
  
  if (setup->nstn == 0) return(0);   /* no work */
  
  // calculate samples per us
  prep_the_data(setup);
  if (setup->sps == 0) return(1);    /* setup problem */
  if (setup->verb > 0) fprintf(stdout, "setup->sps is %d\n", setup->sps);
  //setup->slot = 1.0 / setup->sps;
  
  // set the loop limits based on number of slices
  if (loop_numerology(setup)) return(2);
  
  // setup tones and fill in tone array with
  // accumulated value of all the tones
  setup_tones(setup);
  setup->tone_arr = (double *)malloc(setup->sps * sizeof(double));
  fill_tonearray(setup);
  
  // allocate memory for data array
  setup->fdata = (float *)malloc(setup->sps * sizeof(float));
  setup->data = (double *)malloc(setup->sps * sizeof(double));
  
  for (ns = 0; ns < setup->nstn; ns++) {
  	if (fabricate_stdata(ns, setup)) return(1);
  }
  
  // free the memory allocated to tone and data array
  free(setup->tone_arr);
  free(setup->fdata);
  free(setup->data);
  
  /* report bit states and exit */
  return(close_and_report(setup));
}

/*
 * This is smart enough, perhaps.  The samples/us of each
 * station must evenly divide the value we return.
 *
 * In addition to working out the total samples per us slice,
 * we also set the oversampling and channel bw which are used often.
 */
static int work_out_samples_per_us(SetupInfo *setup)
{
  int tries[] = { 8, 16, 32, 64, 125, 128, 256, 500, 512, 1000, 1024,
                  2000, 2048, 4000, 4096, 8000, 8192, 16000, 16384,
                  32000, 32768, 64000, 65536, 128000, 131072 };
  int nt, rv, ns;
  VdifConf *vdcp;
  if (setup->verb>1) fprintf(stdout, "Working out samples per (us) slice\n");
  for (nt = 0; nt < sizeof(tries)/sizeof(int); nt++) {
    rv = tries[nt];                             /* SlicesFixme */
    for (ns = 0, vdcp = setup->vcnf; ns < setup->nstn; ns++, vdcp++)
        if (rv != (vdcp->srate * (rv / vdcp->srate)) ||
            rv < 2*vdcp->smpus) break;          /* fail */
        else if (setup->verb>1) fprintf(stdout, "  St[%d] %d/%d = %d, pass\n",
            ns, rv, vdcp->srate, rv/vdcp->srate);
    if (ns == setup->nstn) break;
    else if (setup->verb>1) fprintf(stdout, "  St[%d] %d/%d = %d, fail\n",
            ns, rv, vdcp->srate, rv/vdcp->srate);
  }
  rv *= setup->slices;
  rv *= setup->spsmul;

  if (setup->verb>0) fprintf(stdout, "Picked %d samples per us\n", rv);
  for (ns = 0, vdcp = setup->vcnf; ns < setup->nstn; ns++, vdcp++) {
    vdcp->over = rv / vdcp->srate / setup->slices;     /* SlicesFixme */
    vdcp->chbw = setup->slices * vdcp->srate / 2.0;    /* SlicesFixme */
  }
  if (setup->fftsps > rv) {
    rv = setup->fftsps;
    if (setup->verb>0) fprintf(stdout,
        "Using %d samples per %d-us per input\n", rv, setup->slices);
  } else if (setup->fftsps > 0) {
    fprintf(stderr, "Asserting sampling %d < %d, good luck!\n", setup->fftsps,rv);
    rv = setup->fftsps;
  }
  return(rv);
}

/*
 * Prepare to fake all the data:
 *   Initialize each station (and channels)
 */
static void prep_the_data(SetupInfo *setup)
{
  int ss;
  setup->sps = work_out_samples_per_us(setup);

  for (ss = 0; ss < setup->nstn; ss++)
    assert(initialize_station(ss, setup) == 0);
}

/*
 * Based on number of slices, set the loop limits
 */
static int loop_numerology(SetupInfo *setup)
{
  double runt;
  /* loop numerology SlicesFixme */
  setup->one_sec = ONE_SEC_IN_USEC / setup->slices;
  setup->nusr /= setup->slices;
  setup->idur = (int)floor(setup->dura);
  runt = (double)(setup->dura - (setup->idur));
  setup->fdus = (int)rint(runt*(double)(setup->one_sec));
  if (setup->slices * (setup->one_sec) != ONE_SEC_IN_USEC) {
    fprintf(stderr,
        "Processing slices must evenly divide a second\n"
        "You have %d * %d != %d\n", setup->slices, setup->one_sec, ONE_SEC_IN_USEC);
    return(2);
  }
  if (setup->verb>0) fprintf(stdout, setup->verb>1
        ?  "Running %d.%06d secs, reporting every %d %d-us slices\n"
        :  "Running %d.%06d secs\n", setup->idur, (setup->fdus)*setup->slices,
           setup->nusr, setup->slices);
  return(0);
}

/*
 * Close all the file descriptors and report on the bit states.
 */
static int close_and_report(SetupInfo *setup)
{
  int ii, ch;
  uint64_t *bsp;
  double bsch, pc, div = 1e6, rms, mag;
  char *lab = "Ms";
  if (setup->verb>0) for (ii = 0; ii < setup->nstn; ii++)
    for (ch = 0; ch < setup->vcnf[ii].chans; ch++) {
      bsch = (double)(setup->vcnf[ii].bc[ch]);
      lab = (bsch>1.e4) ? "Ms" : "ks";
      div = (bsch>1.e4) ? 1.e6 : 1.e3;
      bsp = setup->vcnf[ii].bs[ch];
      pc = ( bsp[1] + bsp[2] ) * 100.0 / bsch;
      rms = sqrt(setup->vcnf[ii].sq[ch] / bsch);
      mag = setup->vcnf[ii].ab[ch] / bsch;
      fprintf(stdout,
          "BS[%d][%d] %.3f %.3f %.3f %.3f (%.1f%% %.3f %s) %.6f %.6f\n",
          ii, ch, bsp[0] / bsch, bsp[1] / bsch,
                  bsp[2] / bsch, bsp[3] / bsch, pc, bsch / div,
          lab, rms, mag);
    }
  if (setup->verb>0) fprintf(stdout, "Flushing data...");
  for (ii = 0; ii < setup->nstn; ii++) destroy_station(ii, setup);
  if (setup->verb>0) fprintf(stdout, " done.\n");
  return(0);
}

/*
 * eof
 */
