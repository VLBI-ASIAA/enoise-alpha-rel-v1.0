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
 *  cmdparser.h
 *
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 *
 *	This file contains functions to parse the command line arguments
 *  and set the value of user controlled variables.
 */

#ifndef __CMDPARSER_H__
#define __CMDPARSER_H__

#include "enoise.h"

/* parse command line arguments */
int cmdline(int *argc, char ***argv, SetupInfo *setup);
static int options(int argc, char **argv, SetupInfo *setup);

/* support functions to set user controlled variables*/
static int add_tone(char *frqamp, SetupInfo *setup);
static void set_bsmx(int b, SetupInfo *setup);
static void set_dura(double d, SetupInfo *setup);
static void set_reps(double r, SetupInfo *setup);
static void set_thrs(double t, SetupInfo *setup);
static void set_corr(double c, SetupInfo *setup);
static void set_fftn(int a, SetupInfo *setup);
static void set_slices(int a, SetupInfo *setup);
static void set_pathd(char *path, SetupInfo *setup);
static void set_pathr(char *path, SetupInfo *setup);
static void set_skyfreq(double freq, SetupInfo *setup);
static void set_limit(double limit, SetupInfo *setup);
static void set_spsmul(int spsmul, SetupInfo *setup);

static int set_comm(char *c, SetupInfo *setup);
void set_csigma(SetupInfo *setup);

#endif /* __CMDPARSER_H__ */

/*
 * eof
 */
