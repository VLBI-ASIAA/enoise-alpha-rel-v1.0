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
 *	enoise.h
 *
 *	Based on the anoise simulator developed by Geoff Crew
 *	enoise enhances anoise by simulating non-zero-baseline data
 *
 * 	anoise:
 * 	A deviant form of bnoise/vnoise which makes vdif formed data simulating
 * 	VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 * 	data and correlate it against Mark5b-like or other flavors of data.
 */

#ifndef __ENOISE_H__
#define __ENOISE_H__

#include <stdint.h>
#include <fftw3.h>

#define VDIF_HDR_OCTETS         (32/8)
#define MAX_OCTETS_PER_FRAME    1024
#define MAX_OCTETS_PER_PACKET   (VDIF_HDR_OCTETS + MAX_OCTETS_PER_FRAME)

/* original packet packing code required 32 here */
#define MAX_NUM_CHANNELS        256

#define MAX_NUM_TONES 		32
#define MAX_NUM_STATIONS        20
#define EPSILON                 1e-5

#define ORDER                   3

/*
 * A representation for a generic VDIF header
 * Bit representations are considered non-portable,
 * but this is fine for our needs.
 *
 * Note this differs slightly from the vnoise.h version.
 */
typedef struct vdif_ext {
    /* word0 */
    uint32_t secs_inre:30;
    uint32_t legacy:1;
    uint32_t invalid:1;
    /* word1 */
    uint32_t df_num_insec:24;
    uint32_t ref_epoch:6;
    uint32_t UA:2;
    /* word2 */
    uint32_t df_len_octets:24;
    uint32_t num_chan_log2:5 ;
    uint32_t ver:3;
    /* word3 */
    uint32_t stationID:16 ;
    uint32_t threadID:10 ;
    uint32_t bpsm1:5 ;
    uint32_t dt:1 ;
    /* word4 */
    uint32_t magic:24;  /* formerly seql */
    uint32_t edv:8;     /* formerly edv5 */
    /* word5 */
    uint32_t status;    /* formerly seqh */
    /* word6 */
    uint32_t seql;      /* formerly psn */
    /* word7 */
    uint32_t seqh;      /* formerly eud8 */
} VdifExt;

/*
 * Just the bits from above that a human might actually want to supply.
 * Plus all the other information needed to do the deed.
 *
 * The fundamental numerology is per-us (MHz), but the global "slices"
 * variable can provide finer spectral resolution.
 */
typedef struct vdif_cnf {
    /* input information */
    char        *file;                      /* name of file to hold outputdata */
    short       id;                         /* station ID */
    int         epoch;                      /* ref epoch */
    int         esecs;                      /* secs in epoch */
    int         bytes;                      /* data bytes per packet */
    int         chans;                      /* data channels */
    int         srate;                      /* samples per channel per us */
    double      choff;                      /* channel offset (FFT units) */
    double      chgap[MAX_NUM_CHANNELS];    /* channel gap (FFT units) */
    /* working data */
    int         smpus;                      /* samples (all channels) per us */
    int         opp;                        /* octets per packet */
    int         pps;                        /* packets per second */
    FILE        *ofp;                       /* output packet FILE */
    FILE        *osp;                       /* output spectral FILE */
    char        *name;                      /* output spectral file name */
    double      *work;                      /* input time-series work array */
    fftw_complex*workC;                     /* input time-series work array in Complex */
    fftw_complex*spec;                      /* intermediate spectral work array */
    fftw_complex*chop;                      /* intermediate spectral work array */
    fftw_complex*chanC;                     /* output time-series work array in Complex */
    fftw_plan   plan_ws;                    /* FFTW3 plan work->spec */
    fftw_plan   plan_cc;                    /* FFTW3 plan chop->chan */
    double      *chan;                      /* output time-series work array */
    double      sigma;                      /* station sigma */
    double      thresh;                     /* station 2-bit threshold */
    uint64_t    *load;                      /* load point into packet */
    int         shift;                      /* working current shift */
    double      opus;                       /* octets per us */
    int         over;                       /* oversampling factor */
    double      chbw;                       /* channel bandwidth in fft space */
    uint64_t    vhd[MAX_OCTETS_PER_PACKET]; /* packet */
    uint64_t    bs[MAX_NUM_CHANNELS][4];    /* bit states */
    uint64_t    bc[MAX_NUM_CHANNELS];       /* sample count */
    double      sq[MAX_NUM_CHANNELS];       /* signal^2 for RMS */
    double      ab[MAX_NUM_CHANNELS];       /* |signal| for |S| */
    double      tm[MAX_NUM_CHANNELS];       /* threshold multiplier */
} VdifConf;

/* A representation of the setup used by all stations */
typedef struct setup_info {
  double *tone_arr;
  float *fdata;
  double *data;
  //double slot;
  int sps;
  int idur;
  int fdus;
  int one_sec;
  /* user controlled variables */
  int verb;
  double skyfreq;        /* sky frequency*/
  double limit;          /* threashold to calculate delay error */
  int spsmul;            /* sps multiplier */
  float stnoise;         /* ratio of station noise */
  double dura;           /* target duration of data generation */
  uint64_t bsmx;         /* number of states to examine */
  int nusr;              /* number of us between reports */
  double thrs;           /* threshold in sigma units */
  double corr;           /* common component correlation amplitude */
  double csigma;         /* amplitude */
  int fftsps;            /* override samples/us heuristic */
  int slices;            /* number of 1us slices per fft */
  double tphs[MAX_NUM_TONES];
  double tfrq[MAX_NUM_TONES];
  double tamp[MAX_NUM_TONES];
  int ntones;
  char pathd[50];
  char pathr[50];
  /* station info */
  int nstn;              /* number of stations */
  VdifConf vcnf[MAX_NUM_STATIONS];
} SetupInfo;

#endif /* __ENOISE_H__ */

/*
 * eof
 */
