/*--------------------------------------------------------------*/
/*	Simple code for testing grid search algorithms using 	*/
/*	simulate pulses. This code has the option for adding	*/
/*	multiple different errors to the pulses to account	*/
/*	for real life errors:					*/
/*		-> Gaussian noise defined by (mean,sd)		*/
/*		-> Time align errors define by random		*/
/*		   shift in start bin upto max NM		*/
/*		-> Ability to use superpulses to see effects	*/
/*		   of summing noisy pulses to improve signal	*/
/*		-> Ability to run multiple scans per position   */
/*		   and take average deviation			*/
/*		-> Options to use azimuthal segs in chi2	*/
/*		   or use all segs or stick to core + pc +	*/
/*		   hit_seg					*/
/*		-> Option to use drift time cut to reduce	*/
/*		   search space for chi2 calc			*/
/*		-> Drift time range can be set in +- TC ns	*/
/*								*/
/*	Originally written by Jonathan Wright Jan 2017		*/
/*	Editted							*/
/*		-> Feb 2017 - J.Wright				*/
/*		-> Apr 2017 - J.Wright				*/
/*--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "pdecomp.h"
#include "deconv.h"
int drift_time;
int main(int argc, char **argv) {
  clock_t begin, end;
  int i, j, k, aaa, az, test_seg, gs_res, p_index, seg, dt_range;
  int phi_val = 0;
  float chi2_test, chi2_core, chi2_p_cont, chi2_azi, chi2_dummy, chi2, chi2_min, rdiff, pdiff, zdiff, xdiff, ydiff, cart_diff, rz_diff, noise, tt, sig_dummy, azi, trace, core, p_cont;
  float rdiff_temp, pdiff_temp, zdiff_temp, xdiff_temp, ydiff_temp, drift_time;
  float dt_lower = 0.00;
  float dt_upper = 0.94;
  char line[256];
  FILE *file;
  FILE *fid;
  PDecomp pdd;

/* Define parameters in config file */
/* Control parameters --> 1: on
		      --> 0: off */
  int TA   = 1;									// Controls time align error
  int SP   = 1;									// Controld superpulse
  int NP   = 10;                                                                // Number of pulses to generate for superpulse
  float M  = 0.00;								// Mean noise added to pulses
  float SD = 0.02;								// Standard deviation of Gaussian noise - normalised
  int NM   = 2;									// Max number of bins shifted in time align error
  int PA   = 1;									// Multiple searches per position
  int NS   = 10;								// Number of searches per position
  int AZ   = 1;									// Control use of azimuthal segs for chi2
  int AS   = 1;									// Control use of all segs for chi2
  int ST   = 1;									// Control use of segment timing
  int DT   = 1;									// Control use of drift time cut
  int TC   = 50;								// Sets range for drift time cut = dt +- TC

//#define TESTING								// Turn on for TESTING mode - very noisy!

/* Read in configuration file */
  if (argc%2 != 1) {
    printf("Possible options:\n"
           "      -c config_file_name\n");
    return 1;
  }

  for (i=1; i<argc-1; i+=2) {
    if (strstr(argv[i], "-c")) {
      if (!(fid = fopen(argv[i+1], "r"))) {
        printf("\nERROR: config file %s does not exist?\n", argv[i+1]);
        return 1;
      }
      /* read config file */
      printf("\nReading values from config file %s\n", argv[i+1]);
      while (fgets(line, sizeof(line), fid)) {
        // printf("%s", line);
        if ((!strncmp(line, "TA", 2) && (1 != sscanf(line+2," %i", &TA ))) ||
	    (!strncmp(line, "SP", 2) && (1 != sscanf(line+2," %i", &SP ))) ||
            (!strncmp(line, "NP", 2) && (1 != sscanf(line+2," %i", &NP ))) ||
            (!strncmp(line,  "M", 2) && (1 != sscanf(line+2," %f", &M  ))) ||
            (!strncmp(line, "SD", 2) && (1 != sscanf(line+2," %f", &SD ))) ||
	    (!strncmp(line, "NM", 2) && (1 != sscanf(line+2," %i", &NM ))) ||
	    (!strncmp(line, "PA", 2) && (1 != sscanf(line+2," %i", &PA ))) ||
	    (!strncmp(line, "NS", 2) && (1 != sscanf(line+2," %i", &NS ))) ||
	    (!strncmp(line, "AZ", 2) && (1 != sscanf(line+2," %i", &AZ ))) ||
	    (!strncmp(line, "AS", 2) && (1 != sscanf(line+2," %i", &AS ))) ||
            (!strncmp(line, "ST", 2) && (1 != sscanf(line+2," %i", &ST ))) ||
            (!strncmp(line, "DT", 2) && (1 != sscanf(line+2," %i", &DT ))) ||
            (!strncmp(line, "TC", 2) && (1 != sscanf(line+2," %i", &TC )))) {
          printf("ERROR in config file %s\n"
                 "   ...line is: %s", argv[i+1], line);
          return 1;
        }
      }
      fclose(fid);
    } else {
      printf("Possible options:\n"
             "      -c config_file_name\n");
      return 1;
    }
  }

  #ifdef TESTING
    printf("\n\n"									// For use in testing config read
           "Configuration file read with paramaters"
  	   "\n TA = %i \n SP = %i \n NP = %i \n  M = %f"
	   "\n SD = %f \n NM = %i \n PA = %i \n NS = %i"
	   "\n AZ = %i \n AS = %i \n ST = %i \n\n",
	   TA,SP,NP,M,SD,NM,PA,NS,AZ,AS,ST);
  #endif

  if (!read_basis(BASIS_FILE, &pdd)) printf("Success!\n");

  file = fopen("pos_res_result.txt","w+");

  #ifdef TESTING
    fid  = fopen("pulses.txt","w+");
    for(p_index=0;p_index<5;p_index++) {
  #else
    for(p_index=0;p_index<GRID_PTS;p_index++) {
  #endif
    begin = clock();									// Start clock at beginning of an event - CHECK PLACEMENT
    rdiff = pdiff = zdiff = xdiff = ydiff = cart_diff = rz_diff = 0.0;

    if(!PA) NS = 1;
    if(AS) AZ = 0;									// If using all segs, turn off azi segs --> saves computational time

    #ifdef TESTING
      phi_val = 0;									// First 5 events in TESTING are phi=0 events
    #endif

    if(pdd.basis[p_index].ip != phi_val) continue;					// Prevent phi != 0 events from being written out

    for(aaa=0;aaa<NS;aaa++) {
      float hit_sigs[NUM_SIGS][TIME_STEPS_C];                                           // initialise test pulse array and noisy array
      memset(hit_sigs, 0, sizeof(float)*NUM_SIGS*TIME_STEPS_C);
      float noisy_sigs[NUM_SIGS][TIME_STEPS_C];
      memset(noisy_sigs, 0.0, sizeof(float)*TIME_STEPS_C*NUM_SIGS);
      float dt_sig[TIME_STEPS_C];                                           // initialise test pulse array and noisy array
      memset(dt_sig, 0, sizeof(float)*TIME_STEPS_C);

      for(i=0;i<NUM_SIGS;i++) {								// fill test pulse array with basis signals
        for(j=0;j<TIME_STEPS_C;j++) {
          hit_sigs[i][j] = pdd.basis[p_index].signal[i][j];
          test_seg       = pdd.basis[p_index].iseg;
        }
      }

      if(!SP) NP = 1;									// If superpulses are turned off, NP = 1 --> single pulse
      for(j=0; j<TIME_STEPS_C; j++) {							// Else use NP --> superpulse
        for(k=0;k<NP;k++) {
	  if(!AS) {
            if(test_seg!=19) {
              noise = r_norm(M, SD);
              noisy_sigs[test_seg][j] += hit_sigs[test_seg][j] + noise;			// Add noise to hit segment if test_seg!=core
	    }
            noise = r_norm(M, SD);
  	    noisy_sigs[0][j] += hit_sigs[0][j] + noise;					// Add noise pc
            noise = r_norm(M, SD);
            noisy_sigs[19][j] += hit_sigs[19][j] + noise;				// Add noise to core
	    if(AZ) {
	      for(az=1;az<9;az++) {
		if(az == test_seg) continue;						// Skip az if equal to test seg to prevent duplicating
                noise = r_norm(M, SD);
                noisy_sigs[az][j] += hit_sigs[az][j] + noise;				// Add noise to azimuthal segs
              }
            }
	  }
	  if(AS) {
	    for(seg=0;seg<NUM_SIGS;seg++) {
	      noise = r_norm(M, SD);
	      noisy_sigs[seg][j] += hit_sigs[seg][j] + noise;				// Add noise to all segments
	    }
 	  }
        }
	if(SP) for(seg=0;seg<NUM_SIGS;seg++) noisy_sigs[seg][j] /= NP;

        #ifdef TESTING									// Print traces to fid for TESTING
          fprintf(fid,"\n %i %f %f %f",(p_index*TIME_STEPS_C)+j, noisy_sigs[0][j], noisy_sigs[test_seg][j], noisy_sigs[19][j]);
        #endif
      }

      //-------------Time Align Error-----------//

      if(TA) {
	float shift_sigs[NUM_SIGS][TIME_STEPS_C+10];					// Add 10 bins to allow for shift between pulses
	memset(shift_sigs, 0.0, sizeof(float)*NUM_SIGS*(TIME_STEPS_C+10));

	int t_shift = rand_interval(0, NM);						// Generated random number between 0 and NM

	for(k=0;k<t_shift;k++) {							// Set pre t_shift baseline to prevent strange results in chi2 calc
	  for(seg=0;seg<NUM_SIGS;seg++) {
  	    shift_sigs[seg][k] = noisy_sigs[seg][k];
	  }
	}
	for(k=t_shift;k<TIME_STEPS_C+t_shift;k++) {
	  if(!AS) {
            shift_sigs[0][k] = noisy_sigs[0][k-t_shift];
	    if(test_seg !=19) shift_sigs[test_seg][k] = noisy_sigs[test_seg][k-t_shift];
	    shift_sigs[19][k] = noisy_sigs[19][k-t_shift];
            if(AZ) {
  	      for(az=1;az<9;az++) shift_sigs[az][k] = noisy_sigs[az][k-t_shift];
            }
          }
	  if(AS) {
	    for(seg=0;seg<NUM_SIGS;seg++) shift_sigs[seg][k] = noisy_sigs[seg][k-t_shift];
	  }
        }

        for(k=0;k<TIME_STEPS_C;k++) {							// Set noisy_sigs->shift_sigs since chi2 calc still uses noisy_sigs array
	  if(!AS) {
            noisy_sigs[0][k] = shift_sigs[0][k];
	    if(test_seg!=19) noisy_sigs[test_seg][k] = shift_sigs[test_seg][k];
	    noisy_sigs[19][k] = shift_sigs[19][k];
            if(AZ) {
  	      for(az=1;az<9;az++) noisy_sigs[az][k] = shift_sigs[az][k];
	    }
	  }
	  if(AS) {
	    for(seg=0;seg<NUM_SIGS;seg++) noisy_sigs[seg][k] = shift_sigs[seg][k];
	  }
        }
      }

      //-----------Drift Time Calculation-------//
      if(DT) {										// If DT, calc drift time from sample pulse
	for(k=0;k<TIME_STEPS_C;k++) {							// else, set drift time as actual value from basis.dat
	  dt_sig[k] = noisy_sigs[0][k];
        }
	drift_time = calc_dt(dt_sig, dt_lower, dt_upper);
//      printf("\nActual drift time = %f \n New drift time = %f",pdd.basis[p_index].t_drift, drift_time);
      } else drift_time = pdd.basis[p_index].t_drift;

      //-------------Noisy Grid Search----------//
      // Compare pulses from hit segment, pc    //
      // and core to get chi2		        //
      //----------------------------------------//
      chi2 = 0.0;
      chi2_min = 10000.0;
      for(i=0; i<GRID_PTS; i++) {
        chi2 = 0.0;
	if(DT) if(pdd.basis[i].t_drift <= drift_time-TC || pdd.basis[i].t_drift >= drift_time+TC) continue;
        for(j=0;j<NUM_SIGS;j++) {
	  if(!AS) {
	    if(!AZ) {
	      if(j!=0 && j!=19 && j!=test_seg) continue;
	    }
	    if(AZ) {
	      if(j>8 && j!=19 && j!=test_seg) continue;
	    }
	  }
          for(k=0; k<TIME_STEPS_C; k++) {
            sig_dummy = pdd.basis[i].signal[j][k];
	    chi2_dummy = sig_dummy - noisy_sigs[j][k];
	    chi2 += chi2_dummy*chi2_dummy;
          }
	}
	chi2 /= TIME_STEPS_C;

        if(chi2 < chi2_min) {
	  chi2_min = chi2;
	  gs_res = i;
	}
      }

      rdiff_temp = pdd.basis[p_index].ir - pdd.basis[gs_res].ir;				// Needs temp var to prevent cart_diff and rz_diff from adding squared summation
      pdiff_temp = pdd.basis[p_index].ip - pdd.basis[gs_res].ip;				// and only dividing by linear value
      zdiff_temp = pdd.basis[p_index].iz - pdd.basis[gs_res].iz;
      xdiff_temp = pdd.basis[p_index].x - pdd.basis[gs_res].x;
      ydiff_temp = pdd.basis[p_index].y - pdd.basis[gs_res].y;

      rdiff += fabs(rdiff_temp);
      pdiff += fabs(pdiff_temp);
      zdiff += fabs(zdiff_temp);
      xdiff += fabs(xdiff_temp);
      ydiff += fabs(ydiff_temp);

      cart_diff += sqrt(xdiff_temp*xdiff_temp + ydiff_temp*ydiff_temp + zdiff_temp*zdiff_temp);
      rz_diff += sqrt(rdiff_temp*rdiff_temp + zdiff_temp*zdiff_temp);

      #ifdef TESTING
        printf("\nGrid Search Algorithm using Noisy Pulses\n");
        printf("Test pulse id = %d, Grid search result id = %d\n", p_index, gs_res);
        if(rdiff != 0.0 || zdiff != 0.0 || pdiff != 0.0) {
          printf("Rdiff = %f, Pdiff = %f, Zdiff = %f \n", rdiff, pdiff, zdiff);
          printf("RZdiff = %f, XYZdiff = %f \n", rz_diff, cart_diff);
	}
      #endif


    } //aaa PA loop

    end = clock();
    tt = (float)(end-begin)/CLOCKS_PER_SEC;                                     	// Times for a single event to run

    if(p_index == 2) printf("\n\nThe first event took %f s to run\n",tt);

    rdiff /= NS;
    pdiff /= NS;
    zdiff /= NS;
    xdiff /= NS;
    ydiff /= NS;
    cart_diff /= NS;
    rz_diff /= NS;

    #ifdef TESTING
      printf("\n\n\n %hd %hd %hd %f %f %f %f %f\n\n\n",pdd.basis[p_index].ir, pdd.basis[p_index].ip, pdd.basis[p_index].iz, rdiff, pdiff, zdiff, cart_diff, rz_diff);
    #endif

    fprintf(file,"\n %hd %hd %hd %f %f %f %f %f",pdd.basis[p_index].ir, pdd.basis[p_index].ip, pdd.basis[p_index].iz, rdiff, pdiff, zdiff, cart_diff, rz_diff);
  } //p_index loop
  fclose(file);
  #ifdef TESTING
    fclose(fid);
  #endif
} //main loop

//----------Normal Noise-------------------//
double r_norm (double mu, double sigma) {
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;
  return (mu + sigma * (double) X1);
}

//-----------Random Number Generator within limits-----------//
int rand_interval(int min, int max)
{
    int r;
    const unsigned int range = 1 + max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;

    /* Create equal size buckets all in a row, then fire randomly towards
     * the buckets until you land in one of them. All buckets are equally
     * likely. If you land off the end of the line of buckets, try again. */
    do
    {
        r = rand();
    } while (r >= limit);

    return min + (r / buckets);
}
//------------Calcs drift time between lower and upper---------//
static int calc_dt(float *signal, float lower, float upper) {
  int c1, c2, t1, t2, i;
  c1 = c2 = t1 = t2 = 0;

  for(i=0;i<TIME_STEPS_C;i++) {
    float test_sig = fabs(signal[i]);
    if(test_sig >= lower && c1 == 0) {
      t1 = i;
      c1++;
    }
    if(test_sig >= upper && c2 == 0) {
      t2 = i;
      c2++;
    }
    if(c1 != 0 && c2 != 0) drift_time = 10*(t2 - t1);					// Factor 10 accounts for bin size when using 100 MHz digitiser
  }
return drift_time;
}

