"""Script to generate LOFAR antenna response data."""
import sys
#sys.path.append('../..')
import numpy
import re
import pickle
from antpat.dualpolelem import DualPolElem
from antpat.reps.hamaker import HamakerPolarimeter
from native.parseAntennaField import getArrayBandParams, list_stations
from telescopes.rt import TelescopeStnBnd

TELESCOPE_NAME = 'LOFAR'
nr_pols = 2
sampfreq = 100e6
nr_channels = 512
bands = ['LBA', 'HBA']
antmodels = ['Hamaker']
oosr2 = 1./numpy.sqrt(2)
# + -
# + +
polcrdrot = numpy.array([[+oosr2, -oosr2,  0.],
                         [+oosr2, +oosr2,  0.],
                         [    0.,     0.,  1.]])
LOFAR_HAdata_dir = './share/' #Directory for native telescope project data.
TELEDATADIR = 'data/'         #Directory for telescope data for RIME level work.
HA_LBAfile_default = TELEDATADIR+'HA_LOFAR_elresp_LBA.p'
HA_HBAfile_default = TELEDATADIR+'HA_LOFAR_elresp_HBA.p'
DP_LBAfile_default = TELEDATADIR+'DP_model_LBA.p'
DP_HBAfile_default = TELEDATADIR+'DP_model_HBA.p'


def read_LOFAR_HAcc(coefsccfilename):
    """Read Hamaker-Arts coefficients from c++ header files used in the
    "lofar_element_response" code developed at ASTRON for LOFAR. It contains
    LOFAR specific constructs such as reference to "lba" and "hba", so it is not
    suitable for other projects.
    """
    re_fcenter = r'[lh]ba_freq_center\s*=\s*(?P<centerstr>.*);'
    re_frange  = r'[lh]ba_freq_range\s*=\s*(?P<rangestr>.*);'
    re_shape   = r'default_[lh]ba_coeff_shape\[3\]\s*=\s*\{(?P<lstshp>[^\}]*)\}'
    re_hl_ba_coeffs_lst = r'(?P<version>\w+)(?P<band>[hl]ba)_coeff\s*\[\s*(?P<nrelem>\d+)\s*\]\s*=\s*\{(?P<cmplstr>[^\}]*)\}'
    re_cc_cmpl_coef = r'std::complex<double>\((.*?)\)'
    with open(coefsccfilename, 'r') as coefsccfile:
        coefsfile_content = coefsccfile.read()
    searchres = re.search(re_fcenter, coefsfile_content)
    freq_center = float(searchres.group('centerstr'))
    searchres = re.search(re_frange, coefsfile_content)
    freq_range = float(searchres.group('rangestr'))
    searchres = re.search(re_shape, coefsfile_content)
    lstshp = [int(lstshpel) for lstshpel in searchres.group('lstshp').split(',')]
    lstshp.append(nr_pols)
    searchres = re.search(re_hl_ba_coeffs_lst, coefsfile_content, re.M)
    HAcoefversion = searchres.group('version')
    HAcoefband = searchres.group('band')
    HAcoefnrelem = searchres.group('nrelem')
    lstofCmpl = re.findall(re_cc_cmpl_coef, searchres.group('cmplstr'))
    cmplx_lst = []
    for reimstr in lstofCmpl:
        reimstrs = reimstr.split(',')
        cmplx_lst.append( complex(float(reimstrs[0]), float(reimstrs[1])) )
    coefs = numpy.reshape(numpy.array(cmplx_lst), lstshp)
    #The coefficients are order now as follows:
    #coefs[k,theta,freq,spherical-component].shape == (2,5,5,2)
    artsdata={'coefs': coefs, 'HAcoefversion': HAcoefversion,
              'HAcoefband': HAcoefband, 'HAcoefnrelem': HAcoefnrelem,
              'freq_center': freq_center, 'freq_range': freq_range}
    return artsdata


def convLOFARcc2HA(inpfile, outfile, channels):
    """Convert a .cc file of the Hamaker-Arts model to a file with a pickled
    dict of a Hamaker-Arts instance."""
    artsdata = read_LOFAR_HAcc(inpfile)
    artsdata['channels'] = channels
    pickle.dump(artsdata, open(outfile, 'wb'))


def convHA2DPE(inp_HA_file, out_DP_file):
    """Convert a file with a pickled dict of a Hamaker-Arts instance to a
    file with a pickled DualPolElem object."""
    artsdata = pickle.load(open(inp_HA_file, 'rb'))
    HLBA = HamakerPolarimeter(artsdata)
    stnDPolel = DualPolElem(HLBA)
    pickle.dump(stnDPolel, open(out_DP_file, 'wb'))


def gen_antmodelfiles(inpfileL=LOFAR_HAdata_dir+'DefaultCoeffLBA.cc',
                       inpfileH=LOFAR_HAdata_dir+'DefaultCoeffHBA.cc',
                       outfileL=HA_LBAfile_default,
                       outfileH=HA_HBAfile_default
                       ):
    """A convenience function to produce the pickled 'artsdata' for the default
    LOFAR model data stored in the 'lofar_elem_resp' packages c++ header files.
    Also adds nominal LOFAR frequency channels."""
    
    #Adding nominal frequency channels. The HA model for LOFAR has two bands
    #while data recording S/W has 3 intervals based on sampling frequency,
    #namely (0,100), (100,200), (200,300), each with 512 channels.
    #Here I concatenate the two latter intervals.

    channels = numpy.linspace(0., sampfreq, nr_channels, endpoint=False)
    convLOFARcc2HA(inpfileL, outfileL, channels)
    inpfileL = outfileL
    convHA2DPE(inpfileL, DP_LBAfile_default)
    channels = numpy.linspace(sampfreq, 3*sampfreq, 2*nr_channels, endpoint=False)
    convLOFARcc2HA(inpfileH, outfileH, channels)
    inpfileH = outfileH
    convHA2DPE(inpfileH, DP_HBAfile_default)


def create_stnbnd_telescope(stnID, ArrBand, antmodel='Hamaker'):
    """Create a telescope-station-band object for station stnId and band
    ArrBand using antenna model antmodel."""
    #    *Setup station Jones*
    ##Get the metadata of the LOFAR station. stnRot is the transformation matrix
    ##  ITRF_crds = stnRot*LOFAR_crds
    stnPos, stnRot, stnRelPos = getArrayBandParams(stnID, ArrBand)
    ##Create station's antenna model
    if antmodel == 'NEC_TC':
        #This is an example of dual-pol antenna object constructed
        #from two single-pol radiation patterns
        basefile = 'Lofar-dipole-FREE.out'
        fileloc = NECdir+basefile
        atvfd = NECread.readNECout_tvecfuns(fileloc)
        lofarAntX = RadFarField(atvfd)
        lofarAntY = lofarAntX.rotate90()
        lofarpolel = DualPolElem(lofarAntX, lofarAntY)
    elif antmodel == "idl_dip":
        #ideal dipoles:
        tstantX = theoreticalantennas.ideal_dipole_grid_X(np.linspace(10.e6,90e6,11))
        tstantY = tstantX.rotate90()
        stnDPolel = DualPolElem(tstantX, tstantY)
    elif antmodel == 'Hamaker':
        #This is an example of dual-pol element built from a monolithic
        #Jones representation.
        if ArrBand == 'LBA':
            DP_LHBAfile = DP_LBAfile_default
        else:
            DP_LHBAfile = DP_HBAfile_default
        stnDPolel = pickle.load(open(DP_LHBAfile, 'rb'))
    else:
        print("Error no such LOFAR model")
        exit(1)
    #Rotate 45 degrees since LOFAR elements are 45 degrees to meridian:
    stnDPolel.rotateframe(polcrdrot)
    
    #Create a StationBand object for this
    stnbd = TelescopeStnBnd(stnDPolel, stnPos, stnRot)
    return stnbd


def savetelescope(stnlst, antmodel):
    """Save all the data relevant to the telescope beam modelling into
    one file."""
    telescope = {'Name': TELESCOPE_NAME, 'PatternModel': antmodel}
    print("Generating telescope beam model data for:")
    telescope['Station'] = {}
    for stnId in stnlst:
        telescope['Station'][stnId] = {}
        for band in bands:
            print(stnId, band, antmodel)
            stnbndtel = create_stnbnd_telescope(stnId, band, antmodel)
            telescope['Station'][stnId][band] = stnbndtel
    saveName="teldat_"+TELESCOPE_NAME+"_"+antmodel+".p"
    pickle.dump(telescope, open(TELEDATADIR+saveName, 'wb'))


if __name__ == "__main__":
    """Use this to produce telescope data files for use in Myra. Run this as:
    $ python telwizhelper_LOFAR.py
    """
    gen_antmodelfiles()
    stnlst = list_stations()
    antmodel = 'Hamaker'
    savetelescope(stnlst, antmodel)
