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
 *  station.h
 *
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 *
 *	This file contains functions to create and output stations.
 */

#ifndef __STATION_H__
#define __STATION_H__

#include "enoise.h"
#include "stdata.h"

/* station functions */
int new_station(char *station, SetupInfo *setup);
int initialize_station(int ns, SetupInfo *setup);
void destroy_station(int ns, SetupInfo *setup);
void fabricate_station(int cus, SetupInfo *setup, FILE *fn, StInfo *stinfo);
/*
 * Read data from input file streams
 * check whether there are enough numbers left in the file
 * if not, reset the file position back to the beginning of the file
 * @param arr Start position of the array to read data in
 * @param numsamps Number of samples to read
 * @param fs File stream to read data from
 */
int read_data(float *arr, int numsamps, FILE *fs);

/* support functions to create new stations */
static int mylog2(int nc);
static char *compute_derived(char *type, VdifConf *vdcp);
static char *update_type(char *type, VdifConf *vdcp);
static void update_start(char *start, VdifConf *vdcp);
static VdifConf *parse_station(char *station, SetupInfo *setup);
static void starting_header(VdifExt *vh, VdifConf *vdcp, int nstn);

/* support functions to initialize stations */
static void describe_station(VdifConf *vdcp, SetupInfo *setup);

/* support functions to fabricate station data per us */
//static void receive_data(double *data, int sps, double *work, FILE *fn, int verb);
static void receive_data(SetupInfo *setup, FILE *fn, int ns);
//static int packetize_station(double ratediff, SetupInfo *setup, StInfo *stinfo);
static int packetize_station(int cus, SetupInfo *setup, StInfo *stinfo);
static void output_station(int ns, SetupInfo *setup);
//static void channelize(int ch, VdifConf *vc, int sps, double theta, double diff);
//static void channelize(int ch, VdifConf *vc, int sps, int cus, double delay_rate, double lofreq, double theta);
static void channelize(int ch, VdifConf *vc, int sps, int cus, StInfo *stinfo, double lofreq, double theta);
static void load_chan(int ch, VdifConf *vc, int sps, uint64_t bsmx, int verb);
static void update_thresh_mult(int ch, VdifConf *vc);
static void spec_out(int wh, int ch,
              double *dd, fftw_complex *cd, int nn, FILE *sp);
//static void change_phase(fftw_complex *cnum, double theta, double diff);
static void change_phase(fftw_complex *cnum, int idx, int sps, double theta, int cus, double delay_rate, double lofreq);
static double fraction_of(double val);
static void real_to_complex(double* work, fftw_complex* workC, int sps);
static void complex_to_real(fftw_complex* chanC, double* chan, int sps);
static void symmetry_check(fftw_complex* freq, int sps);
static void hilbert_transform(fftw_complex* pSrc, int sps);
static void applyfringerot(fftw_complex *chanC, double *chan, int sps, double theta, int cus, StInfo *stinfo, double lofreq);
static void applyinversefringerot(fftw_complex *chanC, double *chan, int sps, double theta, int cus, StInfo *stinfo, double lofreq);

#endif /* __STATION_H__ */

/*
 * eof
 */

