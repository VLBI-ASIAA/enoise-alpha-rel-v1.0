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
 *  stdata.c
 *
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 *
 *  This file contains functions to fabricate station data.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "stdata.h"
#include "station.h"
#include "enoise.h"

/*
 * fabricate data for every station
 */
int fabricate_stdata(int ns, SetupInfo *setup)
{
  Input *input = malloc(sizeof(Input));
  double dzero;
  int num_d_slots, isec;
  char cname[70], dname[70], nname[70], rname[70];
  StInfo *stinfo = malloc(sizeof(StInfo));

  // initialize struct stinfo
  stinfo->ns = ns;
  stinfo->delay = 0.0;
  stinfo->delay_rate = 0.0;
  stinfo->accel = 0.0;
  stinfo->numsamples = 0;
  stinfo->rem_n_samps = 0;
  stinfo->insertpos = 0;
  stinfo->icount = 0;
  stinfo->insertstat = 0;
  stinfo->remtones = 0;
  stinfo->prev_tone = 0.0;
  stinfo->sampcount = 0;
  stinfo->current_rate = 0.0;


  // construct filename for station ns
  sprintf(cname, "%s/common.dat", setup->pathd);
  sprintf(dname, "%s/delay.dat", setup->pathd);
  sprintf(nname, "%s/noise%d.dat", setup->pathd, ns);
  sprintf(rname, "%s/rate%d.dat", setup->pathr, ns);

  input->fc = fopen(cname, "r");
  input->fd = fopen(dname, "r");
  input->fn = fopen(nname, "r");
  input->fr = fopen(rname, "r");

  if (input->fc == NULL || input->fd == NULL || 
       input->fn == NULL || input->fr == NULL) {
    fprintf(stdout, "Failed to open one or more data file(s).\n");
    return(1);
  }

  if (setup->verb > 0) fprintf(stdout, "Station %d\n", ns);
  // get delay, delay rate and accel for second 0
  //get_delay_and_rate(&(stinfo->delay), &(stinfo->delay_rate), input->fr);
  get_delay_and_rate(stinfo, input->fr);

  if ((stinfo->delay == 0) || (stinfo->delay_rate == 0) || (fabs(stinfo->accel) < 1e-12)) {
    fprintf(stdout, "Delay, delay_rate or accel is not properly calculated ...\n");
    return(1);
  } else {
    fprintf(stdout, "Sec 0: Delay %.9f, Delay_rate %.9f, Accel %.20f\n", stinfo->delay, stinfo->delay_rate, stinfo->accel);
  }

  dzero = stinfo->delay;
  // current_rate at sec 0 is delay_rate
  stinfo->current_rate = stinfo->delay_rate;
  // calculate number_of_samples using the absolute value of current_rate
  stinfo->numsamples = (unsigned long long int)rint((setup->limit / fabs(stinfo->current_rate)) * 1e6) - 1;
  if (setup->verb > 2) fprintf(stdout, "Number of samples to insert at sec 0 is %llu\n", stinfo->numsamples);
  

  // calculate the number of signal samples to skip from the common signal
  //num_d_slots = (int)rint(stinfo->delay/setup->slot);
  num_d_slots = (int)rint(stinfo->delay * setup->sps);

  if (setup->verb > 0) fprintf(stdout, "Number of signal samples to skip from common.dat is %d.\n", num_d_slots);
  if (setup->verb > 2) fprintf(stdout,
    "Common.dat Original file pointer address is %d\n", ftell(input->fc));

  // move the file pointer to skip signal samples due to delay
  fseek(input->fc, num_d_slots * sizeof(float), SEEK_SET);

  if (setup->verb > 2) fprintf(stdout, 
    "Common.dat File pointer address after skipping samples is %d\n", ftell(input->fc));

  // move the tone pointer
  // it has to match the skipped number of signal samples
  stinfo->remtones = setup->sps - (num_d_slots % setup->sps);
  
  for (isec = 0; isec <= setup->idur; isec++) {
    //reset sample counter at the beginning of every second
    stinfo->sampcount = 0;

    if (fabricate_stdata_sec(isec, input, dzero, stinfo, setup)) return(1);
  }

  fclose(input->fc);
  fclose(input->fd);
  fclose(input->fn);
  fclose(input->fr);
  free(input);
  free(stinfo);
  if (setup->verb > 0) fprintf(stdout, "Station %d DONE!\n", ns);
  return(0);
}

/*
 * fabricate data for second isec
 */
static int fabricate_stdata_sec(int isec, Input *input, double dzero, StInfo *stinfo, SetupInfo *setup)
{
  int cus;

  if (isec != 0) {
    //get_delay_and_rate(&(stinfo->delay), &(stinfo->delay_rate), input->fr);
    get_delay_and_rate(stinfo, input->fr);
    if (stinfo->delay == 0 || stinfo->delay_rate == 0) {
      fprintf(stdout, "Delay or delay_rate is not properly calculated ...\n");
      return(1);
    } else {
      fprintf(stdout, "Sec %d: Delay %.9f, Delay_rate %.9f, Accel %.20f\n",
        isec, stinfo->delay, stinfo->delay_rate, stinfo->accel);
    }
  }

  // for every us
  for (cus = 0; cus < setup->one_sec; cus++) {
    if ((isec == setup->idur) && (cus > setup->fdus)) break;

    fabricate_stdata_us(cus, input, stinfo, setup);
  }

  return(0);
}

/*
 * fabricate station data for microsecond cus
 */
static void fabricate_stdata_us(int cus, Input *input,
             StInfo *stinfo, SetupInfo *setup)
{
  int ss;
  fill_fdata(input, stinfo, setup->fdata, setup->sps, setup->verb, setup->limit);
  
  // convert float to double
  for (ss = 0; ss < setup->sps; ss++) {
    *(setup->data + ss) = (double)*(setup->fdata + ss);
  }
  
  // add amplitude to data array
  add_amp(setup->data, setup->sps, setup->csigma, setup->verb);

  // add tones to data array
  add_tones(stinfo, setup);
  fabricate_station(cus, setup, input->fn, stinfo);
  if (setup->verb > 1) fprintf(stdout, "Station %d at us %d is fabricated!\n", stinfo->ns, cus); 
}

/* support functions to fabricate station data per us */

/*
 * fill in data array with random numbers from input files
 */
static void fill_fdata(Input *input, StInfo *stinfo, 
             float *fdata, int sps, int verb, double limit)
{
  int nfilled = 0;

  while (nfilled != sps) {
    if (stinfo->insertstat) {
      // if delay_rate is positive
      // skip one signal sample
      if (stinfo->delay_rate > 0)
        skip_signal(&(stinfo->remtones), sps, input->fc, verb);
      // if delay_rate is negative
      // insert one delay sample
      else
        insert_delay(fdata, &nfilled, sps, stinfo, input->fd, verb);
      
      // set insertstat back to zero after insertion
      stinfo->insertstat = 0;
    }

    insert_signal(fdata, &nfilled, sps, stinfo, input->fc, verb, limit);
  }
}

/*
 * add amplitude to data array
 */
static void add_amp(double *data, int sps, double csigma, int verb)
{
  int ss;

  for (ss = 0; ss < sps; ss++)
  {
    *(data + ss) *= csigma;
  }
  if(verb > 2)
    fprintf(stdout, "Added amplitude csigma value %f\n", csigma);
}

/*
 * add tones to data array
 */
static void add_tones(StInfo *stinfo, SetupInfo *setup)
{
  int ss;
  int count = setup->sps - stinfo->remtones;

  for (ss = 0; ss < setup->sps; ss++) {
    if (ss == (stinfo->insertpos - stinfo->icount + 1)
         && stinfo->icount != 0) {
      *(setup->data + ss) += stinfo->prev_tone;
      stinfo->icount -= 1;
      if(setup->verb > 1)
        fprintf(stdout, "icount is %d\n", stinfo->icount);
    } else {
      *(setup->data + ss) += *(setup->tone_arr + count);
      stinfo->prev_tone = *(setup->tone_arr + count);
      if(setup->verb > 2)
        fprintf(stdout, "Tone value %f at position %d is added to datavalue at position %d\n",
                  stinfo->prev_tone, count, ss);
      stinfo->remtones -= 1;

      if(stinfo->remtones == 0) {
        stinfo->remtones = setup->sps;
        count = 0;
      } else {
        count++;
      }
    }
  }

  if(setup->verb > 1) fprintf(stdout, "Remaining number of tones is %d\n", stinfo->remtones);
  if(setup->verb > 1) fprintf(stdout, "The last tone added is %f\n", stinfo->prev_tone);
}

/*
 * skip one common signal
 */
static void skip_signal(int *remtones, int sps, FILE *fc, int verb)
{
  float temp[1];
  int rc;

  // read data from input file stream
  rc = read_data(temp, 1, fc);
  assert(rc != -1);
  if(rc == 1)
    fprintf(stdout, "Common data file has been reset!\n");

  if (verb > 1) fprintf(stdout, "Skip 1 data sample from common data\n");
  // shift one tone in the tone array
  *remtones -= 1;
  if (*remtones == 0) *remtones = sps;
}

/*
 * insert one delay signal
 */
static void insert_delay(float *fdata, int *nfilled,
              int sps, StInfo *stinfo, FILE *fd, int verb)
{
  int rc;
  // read data from input file stream
  rc = read_data(fdata + *nfilled, 1, fd);
  assert(rc != -1);
  if(rc == 1)
    fprintf(stdout, "Delay data file has been reset!\n");
  if (verb > 1) fprintf(stdout, "Read 1 delay sample to data array\n");

  stinfo->icount++;
  stinfo->insertpos = *nfilled;
  *nfilled += 1;
  stinfo->sampcount += 1;
}

/*
 * insert common signals
 */
static void insert_signal(float *fdata, int *nfilled,
              int sps, StInfo *stinfo, FILE *fc, int verb, double limit)
{
  int rc;

  if (verb > 1) fprintf(stdout, 
    "Remaining number of samples is %d\n", stinfo->rem_n_samps);
  if (stinfo->rem_n_samps == 0)
  {
    stinfo->current_rate = stinfo->delay_rate + 2 * stinfo->accel * (double)stinfo->sampcount * 1e-6 / sps;
    stinfo->numsamples = (unsigned long long int)rint((limit / fabs(stinfo->current_rate)) * 1e6) - 1;
    stinfo->rem_n_samps = stinfo->numsamples;
  }

  if (stinfo->rem_n_samps > (sps - *nfilled)) {
    // read data from input file stream
    rc = read_data(fdata + *nfilled, sps - *nfilled, fc);
    assert(rc != -1);
    if(rc == 1)
      fprintf(stdout, "Common data file has been reset!\n");

    if (verb > 1) fprintf(stdout, 
      "Read %d data samples to data array\n", sps - *nfilled);
    // update remaining number of samples and sample count
    stinfo->rem_n_samps -= (sps - *nfilled);
    stinfo->sampcount += (sps - *nfilled);
    *nfilled = sps;
  } else {
    // read data from input file stream
    rc = read_data(fdata + *nfilled, stinfo->rem_n_samps, fc);
    assert(rc != -1);
    if(rc == 1)
      fprintf(stdout, "Common data file has been reset!\n");
    if (verb > 1) 
      fprintf(stdout, "Read %d data samples to data array\n", stinfo->rem_n_samps);
    // update number of samples filled and sample count
    *nfilled += stinfo->rem_n_samps;
    stinfo->sampcount += stinfo->rem_n_samps;
    stinfo->rem_n_samps = 0;
    stinfo->insertstat = 1;
  }
}

/* 
 * retrieve delay and delay rate of the station
 */
static void get_delay_and_rate(StInfo *stinfo, FILE *fr)
{
  char buff[200];
  char *ptrs;
  char *temp[ORDER];
  size_t idx;

  for(idx = 0; idx < ORDER; idx++)
    temp[idx] = "0.0";

  assert(fgets(buff, sizeof(buff), fr) != NULL);
  ptrs = strtok(buff, " ");

  int num = 0;
  while (ptrs != NULL && num < ORDER)
  {
    temp[num] = ptrs;
    num++;
    ptrs = strtok(NULL, " ");
  }

  stinfo->delay = atof(temp[0]);
  stinfo->delay_rate = atof(temp[1]);
  stinfo->accel = atof(temp[2]);
  fprintf(stdout, "delay is %.9f, delay_rate is %.9f, accel is %.12f\n", stinfo->delay, stinfo->delay_rate, stinfo->accel);
}

/*
 * eof
 */
