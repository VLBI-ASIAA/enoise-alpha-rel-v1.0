enoise-alpha-rel-v1.0
=====================

<h2><em>enoise</em>: non-zero baseline (VLBI) data simulator</h2>

<h5>Author: Zheng Meyer-Zhao</h5>
<h5>Simulator design: Dr. C. Y. Kuo</h5>
<h5>Date:   2014/03/13</h5>

<p>Based on the <em>anoise</em> simulator developed by Geoff Crew,
enoise enhances anoise by simulating non-zero-baseline data.</p>

<p><em>anoise</em> is a deviant form of bnoise/vnoise which makes vdif formed data simulating
VDIF data in a variety of flavors.
The goal here is to simulate ALMA data and correlate it against Mark5b-like or other flavors of data.</p>

<h5>Prerequisite:</h5>
DiFX-2.3.0, GNU Scientific Library (gsl, gsl-devel) and FFTW3 (fftw, fftw-devel)

<h5>The eoisez software consists of the following components:</h5>

<em>random</em>

Generate random numbers used as common signals, delay signals, and station noise.

<em>input</em>

Example input files that will be used by the simulator.
Users can modify the input files for simulations.

<em>genv2dvex.sh</em>

Script to generate .v2d file, vex file, and .calc file to be used by DiFX.

<em>usemodel.py</em> and <em>difxmodel.py</em>

Scripts to calculate delay and delay_rate.

<em>enoise</em>

Generate non-zero baseline simulation data.
