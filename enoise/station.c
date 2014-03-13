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
 *  station.c
 *
 *  Based on the anoise simulator developed by Geoff Crew
 *  enoise enhances anoise by simulating non-zero-baseline data
 *
 *  anoise:
 *  A deviant form of bnoise/vnoise which makes vdif formed data simulating
 *  VDIF data in a variety of flavors.  The goal here is to simulate ALMA
 *  data and correlate it against Mark5b-like or other flavors of data.
 *
 *  This file contains functions to create and output stations.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include "enoise.h"
#include "station.h"

/* station functions */

/*
 * Add a new station from the command-line specification
 */
int new_station(char *station, SetupInfo *setup)
{
  VdifConf *vdcp = parse_station(station, setup);
  if (!vdcp) return(1);
  assert(vdcp == setup->vcnf + setup->nstn);
  if (setup->nstn >= MAX_NUM_STATIONS) return(2);
  vdcp->ofp = fopen(vdcp->file, "w");
  starting_header((VdifExt*)vdcp->vhd, vdcp, setup->nstn);
  setup->nstn++;
  return(vdcp->ofp ? 0 : 1);
}

/*
 * Create work arrays for the station.  The complex allocations here
 * are twice as large as necessary for conservative interface safety.
 * fftw_malloc is identical to malloc (we are told) except for possible
 * special alignment restrictions.
 *
 * We could probably reuse some of these arrays, but the current
 * plan allows for eventual thread parallelism, perhaps.
 *
 * Returns null on malloc errors.
 */
int initialize_station(int ns, SetupInfo *setup)
{
  int ch, st;
  VdifConf *vdcp = &(setup->vcnf[ns]);
  vdcp->sigma = sqrt (1.0 + setup->csigma * setup->csigma);
  vdcp->thresh *= vdcp->sigma;
  vdcp->work = (double *)fftw_malloc(sizeof(double) * setup->sps);
  if (!vdcp->work) return(perror("malloc"),1);
  vdcp->workC = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * setup->sps);
  if (!vdcp->workC) return(perror("malloc"),1);
  vdcp->spec = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * setup->sps);
  if (!vdcp->spec) return(perror("malloc"),1);
  vdcp->chop = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * setup->sps);
  if (!vdcp->chop) return(perror("malloc"),1);
  vdcp->chanC = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * setup->sps);
  if (!vdcp->chanC) return(perror("malloc"),1);
  vdcp->chan = (double *)fftw_malloc(sizeof(double) * setup->sps);
  if (!vdcp->chan) return(perror("malloc"),1);
  /* zero statistics counters */
  for (ch = 0; ch < vdcp->chans; ch++) {
    vdcp->tm[ch] = 1.0;
    vdcp->sq[ch] = 0.0;
    vdcp->ab[ch] = 0.0;
    for (st = 0; st < 4; st++) vdcp->bs[ch][st] = 0;
    vdcp->bc[ch] = 0;
  }
  if (setup->verb>0) describe_station(vdcp, setup);
  vdcp->plan_ws = fftw_plan_dft_1d(setup->sps,
      vdcp->workC, vdcp->spec, FFTW_FORWARD, FFTW_ESTIMATE);      /* always FORWARD */
  vdcp->plan_cc = fftw_plan_dft_1d(setup->sps,
      vdcp->chop, vdcp->chanC, FFTW_BACKWARD, FFTW_ESTIMATE);     /* always BACKWARD */
  if (vdcp->name) vdcp->osp = fopen(vdcp->name, "w");
  return(0);
}

/*
 * End of life tear-down: close the file, free willie.
 */
void destroy_station(int ns, SetupInfo *setup)
{
  VdifConf *vdcp = &(setup->vcnf[ns]);
  fclose(vdcp->ofp);
  fftw_free(vdcp->work);
  fftw_free(vdcp->workC);
  fftw_free(vdcp->spec);
  fftw_free(vdcp->chop);
  fftw_free(vdcp->chan);
  fftw_free(vdcp->chanC);
  fftw_destroy_plan(vdcp->plan_ws);
  fftw_destroy_plan(vdcp->plan_cc);
}

/*
 * Fabricate station data per us
 */
//void fabricate_station(double ratediff, SetupInfo *setup, FILE *fn, StInfo *stinfo)
void fabricate_station(int cus, SetupInfo *setup, FILE *fn, StInfo *stinfo)
{
  //receive_data(setup->data, setup->sps, setup->vcnf[stinfo->ns].work, fn, setup->verb);
  receive_data(setup, fn, stinfo->ns);
  if (setup->verb > 2) {
    int ii = 0;
    for (ii = 0; ii < 5; ii++)
      fprintf(stdout,
        "Value in vcnf.work at postion %d is %f\n", ii, setup->vcnf[stinfo->ns].work[ii]);
  }

  //if (packetize_station(ratediff, setup, stinfo)) output_station(stinfo->ns, setup);
  if (packetize_station(cus, setup, stinfo)) output_station(stinfo->ns, setup);
  
  if (setup->vcnf[stinfo->ns].spec) {
    if (setup->vcnf[stinfo->ns].osp) fclose(setup->vcnf[stinfo->ns].osp);
    if (setup->vcnf[stinfo->ns].name) free(setup->vcnf[stinfo->ns].name);
    setup->vcnf[stinfo->ns].osp = 0;
    setup->vcnf[stinfo->ns].name = 0;
  }
}

/*
 * Read data from input file streams
 * check whether there are enough numbers left in the file
 * if not, reset the file position back to the beginning of the file
 * @param arr Start position of the array to read data in
 * @param numsamps Number of samples to read
 * @param fs File stream to read data from
 */
int read_data(float *arr, int numsamps, FILE *fs)
{
  size_t numread;
  // check whether there is enough random number left
  // if not, move the pointer back to the beginning of the data file
  numread = fread(arr, sizeof(float), numsamps, fs);
  if(numread != numsamps)
  {
    // check whether reached EOF
    // if so, the total number of elements successfully read will be zero
    if(fread(arr + numread, sizeof(float), numsamps - numread, fs) == 0)
    {
      // sets the file position to the beginning of the file
      fprintf(stdout, "Reset data file position!\n");
      rewind(fs);
      assert(fread(arr + numread, sizeof(float), numsamps - numread, fs)
        == (numsamps - numread));
      return 1;
    } else {
      fprintf(stdout, "Error while reading data!\n");
      return -1;
    }
  }
  return 0; 
}


/* support functions to create new stations */

/*
 * Returns log base 2 of the input
 */
static int mylog2(int nc)
{
  int rv = 0;
  while (nc > 1) { rv++; nc>>=1; }
  return(rv);
}

/*
 * The station definition should have set
 *  chans   -- the number of channels (power of 2)
 *  bytes   -- the number of bytes per logical packet
 *  srate   -- the sample rate for a single channel
 *
 * The packet rate is 31250 logical packets per second, but we need
 * to impose the ethernet limit; we divide by 2 until it is ok.
 *
 * Then we need to define the following
 *  smpus   -- total sample rate (leads to bits/second)
 *  opus    -- total octet rate (octets per second)
 *  opp     -- octets per packet
 *  load    -- load point in the packet
 *  shift   -- shift of bits in octet
 *  choff   -- fft-space channel offset (provided by user)
 *
 * SlicesFixme: unfortunately, this also has to come out
 * even on a per-packet basis, which is not guaranteed.
 */
static char *compute_derived(char *type, VdifConf *vdcp)
{
  int i, nb, gapcnt = 0;
  vdcp->pps = 31250;
  while (vdcp->bytes > 8192) { vdcp->bytes /= 2; vdcp->pps *= 2; }
  /* derived quantities */
  vdcp->smpus = vdcp->srate * vdcp->chans;
  vdcp->opus = (vdcp->smpus * 2.0) / 64.0;    /* bits/us / bits/octet */
  vdcp->opp = VDIF_HDR_OCTETS + vdcp->bytes/8;
  vdcp->load = &vdcp->vhd[VDIF_HDR_OCTETS];
  vdcp->shift = 0;
  if (*type == '-' || *type == '+') {
    if (0 < sscanf(type, "%lf%n", &vdcp->choff, &nb)) {
      type += nb;
    } else {
      fprintf(stderr, "Illegal channel shift %s\n", type);
      return(0);
    }
  }
  if (*type == ',') {
    while (0 < sscanf(type, ",%lf%n", &vdcp->chgap[gapcnt+1], &nb)) {
      type += nb;
      gapcnt++;
    }
    fprintf(stdout, "Number of gaps given is %d\n", gapcnt);
    /* if there are more than two channels and only one gap value is given
     * the gap is applied to all channels
     * gap between channel 0 and 1 is stored in vdcp->chgap[1] 
     */
    if (gapcnt == 1 && vdcp->chans > 2) {
      for (i = 2; i < vdcp->chans; i++)
        vdcp->chgap[i] = vdcp->chgap[1];
    } else {
      //fprintf(stderr, "Illegal channel gap %s\n", type);
      assert((gapcnt+1) <= vdcp->chans);
    }
  }
  return(type);
}

/*
 * Parse the command-line specification of each station data type.
 * Returns a pointer to the next part (time).
 */
static char *update_type(char *type, VdifConf *vdcp)
{
  if (!type) {
    fprintf(stderr, "Station type is required\n");
    return((char*)0);
  } else if (!strncmp(type,"vlbi512",7)) {
    type += 7;
    vdcp->chans = 1;
    vdcp->bytes = 8192;
    vdcp->srate = 1024;
  } else if (!strncmp(type,"alma500",7)) {
    type += 7;
    vdcp->chans = 1;
    vdcp->bytes = 8000;
    vdcp->srate = 1000;
  } else if (!strncmp(type,"trad8.0x",8)) {
    type += 8;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    if (vdcp->chans > 99) type++;
    vdcp->bytes = 128 * vdcp->chans;
    vdcp->srate = 16;
  } else if (!strncmp(type,"vlbi32.0x",9)) {
    type += 9;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    vdcp->bytes = 512 * vdcp->chans;
    vdcp->srate = 64;
  } else if (!strncmp(type,"vlbi64.0x",9)) {
    type += 9;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    vdcp->bytes = 1024 * vdcp->chans;
    vdcp->srate = 128;
  } else if (!strncmp(type,"alma62.5x",9)) {
    type += 9;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    vdcp->bytes = 1000 * vdcp->chans;
    vdcp->srate = 125;
  } else if (!strncmp(type,"sma32.0x",8)) {
    type += 8;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    vdcp->bytes = 512 * vdcp->chans;
    vdcp->srate = 64;
  } else if (!strncmp(type,"carma32.0x",10)) {
    type += 10;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    vdcp->bytes = 512 * vdcp->chans;
    vdcp->srate = 64;
  } else if (!strncmp(type,"smto32.0x",9)) {
    type += 9;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    vdcp->bytes = 512 * vdcp->chans;
    vdcp->srate = 64;
  } else if (!strncmp(type,"iram3032.0x",11)) {
    type += 11;
    vdcp->chans = atoi(type++);
    if (vdcp->chans > 9) type++;
    vdcp->bytes = 512 * vdcp->chans;
    vdcp->srate = 64;
  } else {
    fprintf(stderr, "Undefined station type %s\n", type);
    return((char*)0);
  }
  if (vdcp->chans != (1 << mylog2(vdcp->chans))) {
    fprintf(stderr, "%d is not a power of two\n", vdcp->chans);
    return((char*)0);
  }
  if (vdcp->chans > MAX_NUM_CHANNELS) {
    fprintf(stderr, "%d is too many channels\n", vdcp->chans);
    return((char*)0);
  }
  return(type = compute_derived(type, vdcp));
}

/*
 * Epoch of 2012 if the user can't be bothered to supply start time.
 */
static void update_start(char *start, VdifConf *vdcp)
{
  int ns = sscanf(start, ":%d@%d", &vdcp->epoch, &vdcp->esecs);
  time_t now = time(0);
  if (ns < 2) {
    /* 20120101_000000 == 24@0.0000 */
    vdcp->epoch = 24;
    vdcp->esecs  = now - 1325376000;
  }
}

/*
 * Pick apart the station specification, returning a pointer to the
 * completed header, or else a null pointer if there is a problem.
 */
static VdifConf *parse_station(char *station, SetupInfo *setup)
{
  char *nc;
  VdifConf *vdcp = &(setup->vcnf[setup->nstn]);
  memset(vdcp, 0, sizeof(VdifConf));
  vdcp->id = (station[0]) | (station[1]<<8);
  /* station[2] is a : */
  if (station[2] != ':') return(
    fprintf(stderr, "Illegal station specification %s\n", station),
      (VdifConf*)0);
  vdcp->file  = station+3;
  nc = strchr(vdcp->file, ':');
  if (nc) *nc++ = 0;
  nc = update_type(nc, vdcp);
  if (nc) update_start(nc, vdcp);
  else return((VdifConf*)0);
  vdcp->thresh = setup->thrs;        /* allows a per-station threshold */
  /* for fft monitoring */
  if (setup->verb>0) {
    vdcp->name = malloc(strlen(vdcp->file)+20);
    if (vdcp->name) {
      sprintf(vdcp->name, "%s.spec", vdcp->file);
      /* if the file has vdif in the name, replace "vdif" with spec */
      nc = strstr(vdcp->name, "vdif");
      if (nc) strncpy(nc, "spec\0", 5);
    } else {
      perror("malloc");
    }
  }
  return(vdcp);
}

/*
 * Load the header with initial values.  After this, it's just bumping
 * the counters as each packet is output.
 */
static void starting_header(VdifExt *vh, VdifConf *vdcp, int nstn)
{
  memset(vh, 0, sizeof(VdifExt));
  /* word0 */
  vh->legacy = vh->invalid = 0;
  vh->secs_inre = vdcp->esecs;
  /* word1 */
  vh->ref_epoch = vdcp->epoch;
  vh->df_num_insec = 0;
  vh->UA = 0;
  /* word2 */
  vh->df_len_octets = vdcp->opp;
  vh->num_chan_log2 = mylog2(vdcp->chans);
  vh->ver = 0;
  /* word3 */
  vh->stationID = vdcp->id;
  vh->threadID = nstn;
  vh->bpsm1 = 2 - 1;
  vh->dt = 0;
  /* extended header */
  vh->seql = vh->seqh = vh->edv = vh->status = 0;
  vh->magic = 0xbade1f;
}


/* support functions to initialize stations */

static void describe_station(VdifConf *vdcp, SetupInfo *setup)
{
  int i;
  int ns = (vdcp - setup->vcnf);
  assert(vdcp == &(setup->vcnf[ns]));
  fprintf(stdout,
    "Opened station %c%c[%d] to file %s (%s)\n"
    "  with %d sam/us = %d x %d sam/us, %d databytes/packet,\n"
    "  %d pkts/s and %g octets/us %d bytes/packet\n",
      (vdcp->id & 0x00FF), (vdcp->id & 0xFF00) >> 8,
      ns, vdcp->file, vdcp->ofp ? "ok" : "fail",
      vdcp->smpus, vdcp->chans, vdcp->srate, vdcp->bytes,
      vdcp->pps, vdcp->opus, vdcp->bytes + 32);
  fprintf(stdout,
    "  sigma %g, threshold %g, oversampling %d,\n"
    "  channel bw %.3f (%.2f MHz), shift %+.3f\n",
    vdcp->sigma, vdcp->thresh, vdcp->over,
    vdcp->chbw, 1.0, vdcp->choff         /* SlicesFixme */
  );

  fprintf(stdout, "  ");
  for(i = 0; i < vdcp->chans; i++)
    fprintf(stdout, "gap[%d] %.3f  ", i, vdcp->chgap[i]);
  fprintf(stdout, "\n");

  if (setup->verb>1) fprintf(stdout,
    "  Header bytes: %08X %08X %08X %08X (...)\n",
      (vdcp->vhd[0] & 0x00000000FFFFFFFF),
      (vdcp->vhd[0] & 0xFFFFFFFF00000000)>>32,
      (vdcp->vhd[1] & 0x00000000FFFFFFFF),
      (vdcp->vhd[1] & 0xFFFFFFFF00000000)>>32);
  if (setup->verb>1) fprintf(stdout,
    "  Work pointers %p %p %p %p\n",
    vdcp->work, vdcp->spec, vdcp->chop, vdcp->chan);
}
/* support functions to fabricate station data per us*/

/*
 * Receive the common signal with a noisy receiver
 */
//static void receive_data(double *data, int sps, double *work, FILE *fn, int verb)
static void receive_data(SetupInfo *setup, FILE *fn, int ns)
{
  int ss, rc;
  float noise[setup->sps];

  // read data from input file stream
  rc = read_data(noise, setup->sps, fn);
  assert(rc != -1);
  if(rc == 1)
    fprintf(stdout, "Station noise data file has been reset!\n");
 
  for(ss = 0; ss < setup->sps; ss++) {
    *(setup->vcnf[ns].work + ss) = *(setup->data + ss) + setup->stnoise * (double)noise[ss];

    if (setup->verb > 2)
      fprintf(stdout, "data %f plus noise %f is assigned to work with value%f.\n",
              *(setup->data + ss), noise[ss], *(setup->vcnf[ns].work + ss));
  }
  if (setup->verb > 1) fprintf(stdout, "Received data with station noise.\n");
}

/*
 * Channelize the data and populate the packet
 * Returns 0 if packet is incomplete.
 */
//static int packetize_station(double ratediff, SetupInfo *setup, StInfo *stinfo)
static int packetize_station(int cus, SetupInfo *setup, StInfo *stinfo)
{
  VdifConf *vc = setup->vcnf + stinfo->ns;
  int ch;
  double lofreq, phase, theta, diff;

  for (ch = 0; ch < vc->chans; ch++) {
    lofreq = setup->skyfreq;
    /* delay in (us) and lofreq in (MHz) */
    phase = lofreq * stinfo->delay;
    assert(fraction_of(phase) >= 0.0);
    theta = 2 * M_PI * fraction_of(phase);
    //diff = 2 * M_PI * fraction_of(lofreq * ratediff);

    channelize(ch, vc, setup->sps, cus, stinfo, lofreq, theta);
    //channelize(ch, vc, setup->sps, theta, diff);
    load_chan(ch, vc, setup->sps, setup->bsmx, setup->verb);
  }

  if (vc->load == &vc->vhd[vc->opp]) {      /* packet is full */
    vc->load = &vc->vhd[VDIF_HDR_OCTETS];   /* start of data */
    vc->shift = 0;
    return(1);
  }
  return(0);
}

/*
 * As packets are finished, output the data, and bump the header.
 */
static void output_station(int ns, SetupInfo *setup)
{
  VdifExt *hdr = (VdifExt*)setup->vcnf[ns].vhd;
  /* write the packet */
  if (setup->verb>2) fprintf(stdout, "Writing Station-%d packet...", ns);
  fwrite(hdr, sizeof(uint64_t), setup->vcnf[ns].opp, setup->vcnf[ns].ofp);
  if (setup->verb>2) fprintf(stdout,
    " sec %d fr %05d", hdr->secs_inre, hdr->df_num_insec);
  /* boost the counters */
  hdr->seql ++;
  if (hdr->seql == 0) hdr->seqh ++;
  /* boost the timestamp */
  hdr->df_num_insec++;
  if (hdr->df_num_insec == setup->vcnf[ns].pps) {
    hdr->df_num_insec = 0;
    hdr->secs_inre++;
  }
  if (setup->verb>2) fprintf(stdout,
    " ( %d fr %05d)\n", hdr->secs_inre, hdr->df_num_insec);
}

/*
 * Take the time series data in work, fft it to frequency space,
 * top-hat filter for the approprate channel, fft it back, and
 * put the result into the chan array.  The spectral array is
 * supposed to have sps/2 + 1 elements (with the missing elements
 * available as complex conjugates of existing ones since the input
 * is real:  i.e. spec[fq - sps] = conj(spec[fq]) for fq = 0..sps/2 ).
 * (I.e. we grab one slice, but we're implicitly grabbing the conjugate.)
 *
 * Because of the oversampling, our channels are limited to some
 * very low frequency components.  Note these are NZ 0 channels.
 *
 * We need to downshift prior to going back to the time domain
 * (equivalent to placing a LO at the lower edge of the band).
 *
 * We divide by sps because with fftw's normalization FORWARD+BACKWARD
 * is equivalent to multiplying by the number of samples (sps).
 */
//static void channelize(int ch, VdifConf *vc, int sps, double theta, double diff)
static void channelize(int ch, VdifConf *vc, int sps, int cus, StInfo *stinfo, double lofreq, double theta)
{
  int fq, fx, idx;
  int fq_min, fq_max;
  double fq_min_dbl;

  /* fft work array into spectral array */
  if (vc->osp) spec_out(0, ch, vc->work, 0, sps, vc->osp);
  real_to_complex(vc->work, vc->workC, sps);

  fftw_execute(vc->plan_ws);
  if (vc->osp) spec_out(1, ch, 0, vc->spec, sps, vc->osp);

  symmetry_check(vc->spec, sps);

  /* slice this channel to chop array */
  /* we need a copy because fft is destructive */
  memset(vc->chop, 0, sps * sizeof(double));
  //fq_min_dbl = (double)ch * (vc->chbw + vc->chgap[ch]) + vc->choff;
  fq_min_dbl = (double)ch * vc->chbw + vc->choff;
  for(idx = 0; idx <= ch; idx++)
    fq_min_dbl += vc->chgap[idx]; 
  //fprintf(stdout, "rint(%f) is %f\n", fq_min_dbl, rint(fq_min_dbl));
  //fprintf(stdout, "rint(%f) is %f\n", fq_min_dbl+vc->chbw, rint(fq_min_dbl+vc->chbw));
  fq_min = ((int)rint(fq_min_dbl)%2 == 0) ? rint(fq_min_dbl) + 1 : rint(fq_min_dbl);
  fq_max = rint(fq_min + vc->chbw);
  //fq_min = rint(fq_min_dbl) + 1;
  //fq_max = rint(fq_min_dbl + vc->chbw);
  if (fq_min < 1) fq_min = 1;               /* always punt DC */
  if (fq_max >= sps/2) fq_max = sps/2 - 1;  /* stay inbounds */

  //fprintf(stdout, "channel %d: fq_min_dbl is %f, fq_min is %d, fq_max is %d\n", ch, fq_min_dbl, fq_min, fq_max);

  for (fx = 1, fq = fq_min; fq <= fq_max; fq++, fx++) {
    vc->chop[fx] = vc->spec[fq] / (double)sps;
    vc->chop[sps-fx] = creal(vc->chop[fx]) - cimag(vc->chop[fx]) * I;
  }
  symmetry_check(vc->chop, sps);
  /* fft chopped array back to time array */
  if (vc->osp) spec_out(2, ch, 0, vc->chop, sps, vc->osp);
  fftw_execute(vc->plan_cc);
  //complex_to_real(vc->chanC, vc->chan, sps);

  hilbert_transform(vc->chanC, sps);
  applyfringerot(vc->chanC, vc->chan, sps, theta, cus, stinfo, lofreq);

  //for(idx = 0; idx < sps; idx++)
  //{
  //  change_phase(&vc->chanC[idx], idx, sps, theta, cus, delay_rate, lofreq);
  //  vc->chan[idx] = creal(vc->chanC[idx]);
  //}

  //real_to_complex(vc->chan, vc->chanC, sps);
  //applyinversefringerot(vc->chanC, vc->chan, sps, theta, cus, stinfo, lofreq);

  if (vc->osp) spec_out(3, ch, vc->chan, 0, sps, vc->osp);
}

/*
 * Do the bit banging to put the bits from the chan array into
 * the relevant bits in the packet.  The packet can be loaded
 * with garbage--we just preserve the bits we're not overwriting.
 *
 * If there are at most 32 channels, the 3 statements marked if/while
 * can be coded with "if".  With >32 channels, we extend into the
 * next octet (64 bits) and so we may need shift optr more than once.
 * (The original code had "if" and MAX_NUM_CHANNELS equal to 32.)
 *
 * Each slice has sps samples, but since we are oversampled, we
 * step by vc->over to collect only vc->smpus of them.
 */
static void load_chan(int ch, VdifConf *vc, int sps, uint64_t bsmx, int verb)
{
  int tick, shift = vc->shift + 2 * ch;
  double sample, thresh = vc->thresh * vc->tm[ch];
  uint64_t bits, mask, *optr = vc->load;

  while (shift > 63) { shift -= 64; optr++; }   /* if/while */
  for (tick = 0; tick < sps; tick += vc->over) {
    sample = vc->chan[tick];
    /* quantize it to 2 bits */
    if      (sample >  thresh) bits = 0x3;
    else if (sample < -thresh) bits = 0x0;
    else if (sample >     0.0) bits = 0x2;
    else if (sample <     0.0) bits = 0x1;
    else  /* sample == 0.0 */  bits = 0x1;    /* minor bias */
    /* update statistics up to limit */
    if (vc->bc[ch] < bsmx) {
      (vc->bs[ch][bits])++;
      (vc->bc[ch])++;
      vc->sq[ch] += sample*sample;
      vc->ab[ch] += fabs(sample);
    }
    /* install the octet bits at proper location */
    bits <<= shift;
    mask = 0x3;
    mask <<= shift;
    (*optr) &= ~mask;
    (*optr) |= bits;
    /* move to the next bit location for this channel */
    shift += 2 * vc->chans;
    while (shift > 63) { shift -= 64; optr++; } /* if/while */

    if (verb > 2) {
      if (tick < 5) {
        printf("sample is %f, thresh is %f\n"
               "vc->ab[%d] is %f\n", sample, thresh, ch, vc->ab[ch]);
      }
    }
  }

  if (vc->bc[ch] <= bsmx) update_thresh_mult(ch, vc);

  if (ch == (vc->chans-1)) {
    /* undo previous advancement, then forward 2 bits */
    shift -= 2 * (vc->chans - 1);
    while (shift < 0) { shift += 64; optr--; }  /* if/while */
    vc->shift = shift;          /* starting ch 0 shift */
    vc->load = optr;            /* starting ch 0 load */
  }
}

/*
 * If the channelized data departs significantly from gaussian,
 * the channel thresholds will need adjustment to achieve full
 * population of the 4 bit states.
 *
 * Here we compute and capture a running average rms value of the
 * sample data seen so far.  This value is used later to rescale
 * the samples so that the rms of the frequency channels is the
 * same as the rms of the original data.
 *
 * It should be possible to predict what the scaling factor is,
 * but it's simpler to just compute it.
 */
static void update_thresh_mult(int ch, VdifConf *vc)
{
  double rms = sqrt(vc->sq[ch] / (double)(vc->bc[ch]));
  vc->tm[ch] = rms;
}

/*
 * As a debugging aid, dump out the samples and spectra to the file
 * if it has been opened.  Once we've cycled through all the channels
 * on the first slice, it will be closed.
 */
static void spec_out(int wh, int ch,
    double *dd, fftw_complex *cd, int nn, FILE *sp)
{
  static char *ats[4] = {
      "input data", "spectral data", "channel spec", "output data" };
  static int index = 0;
  int ii;
  fprintf(sp, "# index %d\n", index++);
  fprintf(sp, "# %s, channel %d with %d per us\n", ats[wh], ch, nn);
  fprintf(sp, "# %d rows of %s\n", (wh == 0 || wh == 3) ? nn : nn/2,
    (wh == 0 || wh == 3)
        ? "index(time) data[index]"
        : "index(freq) cabs[index] carg[index]"
  );
  fprintf(sp, "#\n");
  if (wh == 0 || wh == 3) for (ii = 0; ii < nn; ii++)
    fprintf(sp, "%8d %e\n", ii, dd[ii]);
  else for (ii = 0; ii < nn/2; ii++)
    fprintf(sp, "%8d %e %e\n", ii, cabs(cd[ii]), carg(cd[ii])*180.0/M_PI);
  fprintf(sp, "\n\n");
}

static void applyfringerot(fftw_complex *chanC, double *chan, int sps, double theta, int cus, StInfo *stinfo, double lofreq)
{
  int idx;
  double diff, time, total;
  for(idx = 0; idx < sps; idx++)
  {
    total = stinfo->delay + stinfo->delay_rate * (cus + (double)idx / sps) / 1e6
                              + stinfo->accel * pow(cus + (double)idx / sps, 2) / 1e12;
    total /= 1e6;
    time = (cus + (double)idx / sps) / 1e6 + total;
    diff = 2 * M_PI * fraction_of(lofreq * (stinfo->delay_rate * time + stinfo->accel * pow(time, 2)));
    //diff = 2 * M_PI * fraction_of(lofreq * (stinfo->delay_rate * (cus + (double)idx / sps) / 1e6
    //                          + stinfo->accel * pow(cus + (double)idx / sps, 2) / 1e12));
    //diff = 2 * M_PI * lofreq * stinfo->delay_rate * cus / 1e6;
    //diff = 0;
    //theta = 0;
    //if(idx == 0)
    //  printf("theta+diff at cus %d is %.12f\n", cus, theta+diff);
    chan[idx] = cabs(chanC[idx]) * cos(carg(chanC[idx]) + theta + diff);
  }
}

static void applyinversefringerot(fftw_complex *chanC, double *chan, int sps, double theta, int cus, StInfo *stinfo, double lofreq)
{
  int idx;
  double diff;
  for(idx = 0; idx < sps; idx++)
  {
    diff = 2 * M_PI * fraction_of(lofreq * (stinfo->delay_rate * (cus + (double)idx / sps) / 1e6
                              + stinfo->accel * pow(cus + (double)idx / sps, 2) / 1e12));
    //diff = 2 * M_PI * lofreq * stinfo->delay_rate * cus / 1e6;
    chan[idx] = cabs(chanC[idx]) * cos(carg(chanC[idx]) - theta - diff);
  }
}
   

//static void change_phase(fftw_complex *cnum, double theta, double diff)
static void change_phase(fftw_complex *cnum, int idx, int sps, double theta, int cus, double delay_rate, double lofreq)
{
  double r, p, np, a, b, diff;
  double complex z;
  //theta = 0;

  //printf("val is %.9f, phasefrac is %.9f,"
  //"2*PI*phasefrac is %.9f\n", val, phasefrac, 2*M_PI*phasefrac);
  diff = 2 * M_PI * lofreq * delay_rate * (cus + (double)idx / sps) / 1e6;
  //if(cus < 4)
  //  printf("diff at cus %d idx %d is %.9f\n", cus, idx, diff);
  // In polar form
  r = cabs(*cnum);
  p = carg(*cnum);
  //printf("phase is %.9f\n", p);
  np = p + theta + diff;
  a = r * cos(np);
  b = r * sin(np);

  // In rectangular form
  // complex multiplication
  // the phase change is: cos(theta)-isin(theta)
  //a = creal(*cnum) * cos(theta) - cimag(*cnum) * (-sin(theta));
  //b = cimag(*cnum) * cos(theta) + creal(*cnum) * (-sin(theta)); 

  z = a + I * b;
  *cnum = (fftw_complex) z;
  //printf("phase after addition is %.9f\n", carg(z));
  //printf("phase after addition is %.9f\n", carg(*cnum));
}

static double fraction_of(double val)
{
  return val - rint(val - 0.5);
}

/*
 * Convert an array of real numbers to complex
 */
static void real_to_complex(double* work, fftw_complex* workC, int sps)
{
  int idx;
  for(idx = 0; idx < sps; idx++)
  {
    workC[idx] = (fftw_complex) work[idx] + 0.0 * I;
  }
}

/*
 * Convert an array of complex numbers to real
 * The imaginary part of the complex numbers has to be 0
 */
static void complex_to_real(fftw_complex* chanC, double* chan, int sps)
{
  int idx;
  for(idx = 0; idx < sps; idx++)
  {
    assert(cimag(chanC[idx]) < EPSILON);
    //printf("%f\n", cimag(chanC[idx]));
    chan[idx] = creal(chanC[idx]);
  }
}

/*
 * Check whether the freqency spectrum array is symmetric
 */
static void symmetry_check(fftw_complex* freq, int sps)
{
  int idx;
  for(idx = 1; idx < sps / 2; idx++)
  {
    assert((creal(freq[idx]) - creal(freq[sps-idx])) < EPSILON);
    assert((cimag(freq[idx]) + cimag(freq[sps-idx])) < EPSILON);
    //printf("%d %f == %d %f\n", idx, creal(freq[idx]), sps-idx, creal(freq[sps-idx]));
  }
}

static void hilbert_transform(fftw_complex* pSrc, int sps)
{
  fftw_complex* temp;
  int i;
  double r, p;
  complex z;
  temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sps);
  fftw_plan plan_fd = fftw_plan_dft_1d(sps, pSrc, temp, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan plan_bd = fftw_plan_dft_1d(sps, temp, pSrc, FFTW_BACKWARD, FFTW_ESTIMATE);

  // DFT
  fftw_execute(plan_fd);
  // Xa(f) =
  // 2X(f), for f>0
  // X(f),  for f=0,
  // 0,     for f<0
  for(i = 1; i < (sps/2); i++)
  {
    temp[i] *= 2.0;
    temp[sps - i] *= 0.0;
  }

  // inverse DFT
  // now pSrc is filled with the analytic signal
  fftw_execute(plan_bd);

  // normalization
  for(i = 0; i < sps; i++)
  {
    pSrc[i] /= sps;
  }

  fftw_free(temp);
  fftw_destroy_plan(plan_fd);
  fftw_destroy_plan(plan_bd);
}

/*
 * eof
 */
