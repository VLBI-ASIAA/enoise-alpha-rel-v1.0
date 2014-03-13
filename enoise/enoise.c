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

/**
 *  enoise.c
 *
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 *
 *  This file contains main() with basic usage.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "enoise.h"
#include "cmdparser.h"
#include "station.h"
#include "maker.h"

static int usage(char *name)
{
  printf("Usage: %s [options] stn1 [stn2 ...]\n\n", name);
  printf("where the options are:\n\n"
      "  -v            verbose, may be repeated for more\n"
      "  -b <int>      samples used for bit statistics (1024^2)\n"
      "  -d <float>    duration of observation (2.0)\n"
      "  -n <string>   make adjustments to common noise;\n"
      "                use \"help\" for details\n"
      "  -r <float>    report processing interval (usecs)\n"
      "  -t <float>    2-bit threshold in sigma units (1.00)\n"
      "\n"
      "The -r flag produces a progress report for the \"duration\".\n"
      "The -t option is not implemented.\n"
      "\n"
      "Use \"help\" as a station name for station configuration options.\n"
      "\n"
      "A sample invocation giving a (very) short sample with 2 channels\n"
      "of ALMA data (station AL) and 4 channels of typical VLBI data\n"
      "(station SP) is:\n"
      "\n"
      "enoisez -v -d 0.00064 -n corr:0.05 \\\n"
      " -n tone:5,0.01 -n tone:40,0.01 -n tone:75,0.01 -n tone:110,0.01 \\\n"
      " -n pathd:/path/to/ndata -n pathr:/path/to/drate \\\n"
      " AL:AL.vdif:alma62.5x2,1.5:24@6415080 SP:SP.vdif:vlbi32.0x4:24@6415080\n"
      "\n"
      "If -d is 0.0, the first packet is generated, but no\n"
      "data is created; which allows you to check what you will eventually\n"
      "get without waiting a long time for it....\n"
  );
  return(0);
}

/*
 * Help message for station choices below
 */
static void print_type_help(void)
{
  printf(
    "The station specification is a string of the form:\n"
    "\n"
    "  ID:file[:type[choff][:time]]\n"
    "\n"
    "where the type specifies the VDIF packet format for 2-letter\n"
    "station ID, written to file with one of these types:\n"
    "\n"
    "  vlbi512      one single-channel of 512 MHz\n"
    "  alma500      one single-channel of 500 MHz\n"
    "  trad8.0xN    N (256) channels of 8.0 MHz\n"
    "  vlbi32.0xN   N (64) channels of 32.0 MHz\n"
    "  vlbi64.0xN   N (32) channels of 64.0 MHz\n"
    "  alma62.5xN   N (32) channels of 62.5 MHz\n"
    "  sma32.0xN    N (64) channels of 32.0 MHz\n"
    "  carma32.0xN  N (64) channels of 32.0 MHz\n"
    "  smto32.0xN   N (64) channels of 32.0 MHz\n"
    "  iram3032.0xN N (64) channels of 32.0 MHz\n"
    "\n"
    "where in the multi-channel cases the number of channels\n"
    "needed for ~ 2 GHz sampling is indicated in parentheses.\n"
    "The number of channels is required to be a power of 2 and\n"
    "at most %d; and in all cases, 2-bit sampling is assumed.\n",
    MAX_NUM_CHANNELS
  );
  printf(
    "\n"
    "A channel specification may follow.  It can be either or both of\n"
    "    +C or -C    to shift all channels in FFT space\n"
    "    ,G          to insert a gap of G spectral points\n"
    "where the units of C and G are MHz (most likely)\n"
    "In the -C case, the width of the zeroth channel will be reduced\n"
    "so that only one side band is used.  If both are given, the shift\n"
    "must come first.  Both are floating point quantities, but rounding\n"
    "in FFT space may not give you exactly what you want.\n"
    "\n"
    "Finally, the (starting) time specification is epoch@seconds.\n"
    "\n"
  );
}


/*
 * Main Entry.
 */
int main(int argc, char **argv)
{
  int errs = 0;
  SetupInfo *setup = malloc(sizeof(SetupInfo));  

  /* default value for user controlled variables */
  setup->sps = 0;
  setup->idur = 0;
  setup->fdus = 0;
  setup->one_sec = 0;
  setup->verb = 0;
  setup->skyfreq = 230000.0;
  setup->limit = 1.0;
  setup->spsmul = 1;
  setup->stnoise = 1.0;
  setup->dura = 2.00;
  setup->bsmx = 1024 * 1024;
  setup->nusr = 1000;
  setup->thrs = 1.00;
  setup->corr = 0.01;
  setup->fftsps = 0;
  setup->slices = 1;
  setup->ntones = 0; 
  strcpy(setup->pathd, "../ndata");
  strcpy(setup->pathr, "../drate");
  /* initialize tone-related arrays */
  int ii = 0;
  for (ii = 0; ii < MAX_NUM_TONES; ii++) {
    setup->tphs[ii] = 0;
    setup->tfrq[ii] = 0;
    setup->tamp[ii] = 0;
  }

  if ((argc == 1) || !strcmp(argv[1], "--help")) return(usage(argv[0]));
  if (cmdline(&argc, &argv, setup)) return(1);

  if (argc && !strcmp(*argv, "help"))
  {
    print_type_help();
    return(0);
  }
 
  set_csigma(setup);

  while (argc-- && !errs)
    errs += new_station(*argv++, setup);

  if (!errs) errs += fake_the_data(setup);
  
  free(setup);
  return(errs);
}

/*
 * eof
 */ 
