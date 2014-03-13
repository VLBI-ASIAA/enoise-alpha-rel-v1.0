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
 *  stdata.h
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

#ifndef __STDATA_H__
#define __STDATA_H__

#include "enoise.h"

/* A representation of station info at certain second*/
typedef struct stinfo_sec {
	int ns;
	double delay;
	double delay_rate;
  double accel;
  double current_rate;
	/* remaining number of samples to insert */	
	unsigned long long int numsamples;
	unsigned long long int rem_n_samps;
  int sampcount;
	/* info of delay insertion*/
	int insertpos;
	int icount;
	int insertstat;
	int remtones;
	double prev_tone;
} StInfo;

/* A representation of input data files*/
typedef struct input {
	FILE *fc;
	FILE *fd;
	FILE *fn;
	FILE *fr;
} Input;

/* fabricate station data*/
int fabricate_stdata(int ns, SetupInfo *setup);

/* fabricate station data for second isec */
static int fabricate_stdata_sec(int isec, Input *input,
             double dzero, StInfo *stinfo, SetupInfo *setup);

/* fabricate station data for microsecond cus */
static void fabricate_stdata_us(int cus, Input *input,
             StInfo *stinfo, SetupInfo *setup);

/* support functions to fabricate station data per us */
static void fill_fdata(Input *input, StInfo *stinfo, 
             float *fdata, int sps, int verb, double limit);
static void add_amp(double *data, int sps, double csigma, int verb);
static void add_tones(StInfo *stinfo, SetupInfo *setup);
static void skip_signal(int *remtones, int sps, FILE *fc, int verb);
static void insert_delay(float *fdata, int *nfilled,
              int sps, StInfo *stinfo, FILE *fd, int verb);
static void insert_signal(float *fdata, int *nfilled,
              int sps, StInfo *stinfo, FILE *fc, int verb, double limit);

/* retrieve station delay and delay rate*/
//static void get_delay_and_rate(double *delay, double *delay_rate, FILE *fr);
static void get_delay_and_rate(StInfo *stinfo, FILE *fr);

#endif /* __STDATA_H__ */

/*
 * eof
 */

