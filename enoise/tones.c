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
 *  tones.c
 *
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 *
 *  This file contains functions to setup and use tones.
 */

#include <stdio.h>
#include <math.h>

#include "tones.h"

/*
 *  setup tones
 */
void setup_tones(SetupInfo *setup)
{
  int tt;
  double MHz_mult;

  MHz_mult = 2 * M_PI / setup->sps;
  for (tt = 0; tt < setup->ntones; tt++) setup->tfrq[tt] *= MHz_mult;
}

/*
 *  create tone array
 */
void fill_tonearray(SetupInfo *setup)
{
	int ss, tt;

  if(setup->verb > 1) fprintf(stdout,
                       "Fill in tone array with acumulated value of all tones!\n");

  for (ss = 0; ss < setup->sps; ss++) {
    for (tt = 0; tt < setup->ntones; tt++)
      *(setup->tone_arr + ss) += setup->tamp[tt] * sin(setup->tfrq[tt]*ss);
    if(setup->verb > 2) fprintf(stdout,
                         "Accumulated tone at position %d is %f\n", ss, *(setup->tone_arr + ss));
  }
}

/*
 * eof
 */
