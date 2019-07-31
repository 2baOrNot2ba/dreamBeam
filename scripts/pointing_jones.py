#!/usr/bin/env python
"""Model of a LOFAR station. Gets Jones matrix towards a given direction
   and frequency.


   Example:
   $ pointing_jones.py print LOFAR LBA SE607 Hamaker 2012-04-01T01:02:03 60 1 6.11 1.02 60E6

   This prints out the Jones matrices for the LOFAR LBA antenna at 60.e6 Hz for
   the SE607 station tracking a source at RA-DEC 6.11 1.02 (radians) for 60s
   starting at 2012-04-01T01:02:03 using the Hamaker model.
"""
import sys
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from antpat.dualpolelem import plot_polcomp_dynspec
from dreambeam.rime.scenarios import on_pointing_axis_tracking, compute_paral
from dreambeam.telescopes.rt import TelescopesWiz


SCRIPTNAME = sys.argv[0].split('/')[-1]
USAGE = """Usage:\n  {} print|plot telescope band stnID beammodel beginUTC \
duration timeStep pointingRA pointingDEC [frequency]""".format(SCRIPTNAME)

# Startup a telescope wizard
TW = TelescopesWiz()


def printJonesFreq(timespy, Jnf):
    # Select one frequency
    print "Time, Freq, J11, J12, J21, J22"  # header for CSV
    for ti in range(len(timespy)):
        # Creates and prints a comma separated string
        jones_1f_outstring = ",".join(map(
                            str, [
                             timespy[ti].isoformat(), freq,
                             Jnf[ti, 0, 0], Jnf[ti, 0, 1],
                             Jnf[ti, 1, 0], Jnf[ti, 1, 1]]))
        print jones_1f_outstring
        # Print out data for BST-mode comparison (ie powers of p & q channels):
        #print("{0} {1} {2}".format( (timespy[ti]-timespy[0]).total_seconds(), np.abs(Jnf[ti,0,0])**2+np.abs(Jnf[ti,0,1])**2, np.abs(Jnf[ti,1,0])**2+np.abs(Jnf[ti,1,1])**2) )


def plotJonesFreq(timespy, Jnf):
    p_ch = np.abs(Jnf[:, 0, 0].squeeze())**2+np.abs(Jnf[:, 0, 1].squeeze())**2
    q_ch = np.abs(Jnf[:, 1, 1].squeeze())**2+np.abs(Jnf[:, 1, 0].squeeze())**2
    # p_ch = -np.real(Jnf[:,0,0].squeeze())   #+ 0,1 => 0,0
    # q_ch = -np.real(Jnf[:,1,1].squeeze())  #- 1,0 0> 1,1
    # In dB:
    # p_ch = 10*np.log10(p_ch)
    # q_ch = 10*np.log10(q_ch)
    plt.figure()
    plt.subplot(211)
    plt.plot(timespy, p_ch)
    plt.title('p-channel')
    plt.subplot(212)
    plt.plot(timespy, q_ch)
    plt.title('q-channel')
    plt.xlabel('Time')
    plt.show()


def printAllJones(timespy, freqs, Jn):
    """Print all the Jones matrices over time & frequency."""
    print "Time, Freq, J11, J12, J21, J22"  # header for CSV
    # duration.seconds/ObsTimeStp.seconds
    for ti in range(0, len(timespy)):
        for fi, freq in enumerate(freqs):
            # Create and print a comma-separated string
            jones_nf_outstring = ",".join(map(
                                 str, [
                                    timespy[ti].isoformat(), freq,
                                    Jn[fi, ti, 0, 0], Jn[fi, ti, 0, 1],
                                    Jn[fi, ti, 1, 0], Jn[fi, ti, 1, 1]]))
            print jones_nf_outstring


def plotAllJones(timespy, freqs, Jn):
    """Plot all the Jones matrices over time & frequency."""
    plot_polcomp_dynspec(timespy, freqs, Jn)


def main(telescopeName, band, antmodel, stnID, bTime, duration, stepTime,
         CelDir, freq=None):
    """An python entry_point for the pointing_jones command."""
    # Get the telescopeband instance:
    telescope = TW.getTelescopeBand(telescopeName, band, antmodel)
    # Compute the Jones matrices
    timespy, freqs, Jn, srcfld, res, pjonesOfSrc = \
        on_pointing_axis_tracking(telescope, stnID, bTime, duration, stepTime,
                                  CelDir, xtra_results=True)
    if (freq < freqs[0] or freq > freqs[-1]) and freq is not None:
        raise ValueError("Requested frequency {} Hz outside of band {}"
                         .format(freq, band))
    # Do something with resulting Jones according to cmdline args
    if action == "plot":
        stnrot = telescope['Station'][stnID].stnRot
        compute_paral(srcfld, stnrot, res, pjonesOfSrc, bTime)
    if freq is None:
        if action == "plot":
            plotAllJones(timespy, freqs, Jn)
        else:
            printAllJones(timespy, freqs, Jn)
    else:
        frqIdx = np.where(np.isclose(freqs, freq, atol=190e3))[0][0]
        Jnf = Jn[frqIdx, :, :, :].squeeze()
        if action == "plot":
            plotJonesFreq(timespy, Jnf)
        else:
            printJonesFreq(timespy, Jnf)


if __name__ == "__main__":
    # Process cmd line arguments
    args = sys.argv[1:]
    try:
        try:
            action = args.pop(0)
        except IndexError:
            raise RuntimeError("Specify output-type:\n  'print' or 'plot'")
        try:
            telescope = args.pop(0)
        except IndexError:
            raise RuntimeError("Specify telescope:\n  " \
                               + ', '.join(TW.get_telescopes()))
        try:
            band = args.pop(0)
        except IndexError:
            raise RuntimeError("Specify band/feed:\n  " \
                               + ', '.join(TW.get_bands(telescope)))
        try:
            stnID = args.pop(0)
        except IndexError:
            raise RuntimeError("Specify station-ID:\n  " \
                               + ', '.join(TW.get_stations(telescope, band)))
        try:
            antmodel = args.pop(0)
        except IndexError:
            raise RuntimeError("Specify beam-model:\n  " \
                               + ', '.join(TW.get_beammodels(telescope, band)))

        try:
            btime = datetime.strptime(args[0], "%Y-%m-%dT%H:%M:%S")
        except IndexError:
            raise RuntimeError("Specify start-time (UTC in ISO format: yyyy-mm-ddTHH:MM:SS )")
        except ValueError:
            raise RuntimeError("Wrong start-time format (yyyy-mm-ddTHH:MM:SS).")
        try:
            duration = timedelta(0, float(args[1]))
        except (IndexError,ValueError):
            raise RuntimeError("Specify duration (in seconds).")
        try:
            steptime = timedelta(0, float(args[2]))
        except (IndexError, ValueError):
            raise RuntimeError("Specify step-time (in seconds).")
        try:
            celdir = (float(args[3]), float(args[4]), 'J2000')
        except (IndexError, ValueError):
            raise RuntimeError("Specify pointing direction (in radians): RA DEC")
        if len(args) > 5:
            try:
                freq = float(args[5])
            except ValueError:
                raise RuntimeError("Specify frequency (in Hz).")
        else:
            freq = None
    except RuntimeError as mess:
        print(mess)
        print(USAGE)
        sys.exit(2)

    main(telescope, band, antmodel, stnID, btime, duration, steptime, celdir,
         freq)
