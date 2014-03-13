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
 *  maker.h
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

#ifndef __MAKER_H__
#define __MAKER_H__

#define ONE_SEC_IN_USEC 1000000

#include "enoise.h"

/* generate fake data */
int fake_the_data(SetupInfo *setup);

/* support functions to fake the data */
static int work_out_samples_per_us(SetupInfo *setup);
static void prep_the_data(SetupInfo *setup);
static int loop_numerology(SetupInfo *setup);
static int close_and_report(SetupInfo *setup);

#endif /* __MAKER_H__ */

/*
 * eof
 */

