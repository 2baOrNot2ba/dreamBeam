#!/usr/bin/env python
"""A module to read LOFAR antenna field files.
   The ROTATION_MATRIX field is the rotation matrix of the station;
   it maps the the local station coordinates to the ITRF coordinates
   such that
     r_ITRF = ROTATION_MATRIX * r_stn where r are 3-element column vectors.
"""
#This file was taken from the mscorpol.py modules. (It has been modified.)
import sys
import numpy as np
import os
from os.path import dirname
import argparse

#ANTENNAFIELDDIR='src/NDPPP/LOFAR/MAC/Deployment/data/StaticMetaData/AntennaFields/'
STATICMETADATA = dirname(__file__)+'/../share/StaticMetaData/'
ANTENNAFIELDDIR = STATICMETADATA
IHBADELTASDIR=dirname(__file__)+'/../share/iHBADeltas/'
IHBADELTASDIR=STATICMETADATA
COMMENT_CHAR = '#'
AFfileNameType=2


def _getAntennaFieldFile(stationName ):
    antenna_field_dir = ANTENNAFIELDDIR
    if AFfileNameType==2 or AFfileNameType==3:
       basename=stationName+'-'+'AntennaField'+'.conf'
    else:
       basename='AntennaField'+stationName+'.conf'
    filepath=antenna_field_dir+'/'+basename
    print filepath
    return filepath


def _getiHBADeltafile(stationName):
    """Get file path to iHBADelta for given stationName. """
    basename=stationName+'-'+'iHBADeltas'+'.conf'
    filepath=IHBADELTASDIR+'/'+basename
    return filepath


def parseAntennaField(stationName):
    filepath=_getAntennaFieldFile(stationName)
    return parseAntennaFieldFile(filepath)


def parseAntennaFieldFile(filename, AFfileNameType=2):
    """Parse LOFAR AntennaField file by name and return data as a dict.
    Note that this reads in both LBA and HBA parameters.
    """
    AntFldData={'LBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA0': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA1': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]}
               }
    try:
        f = open(filename)
    except IOError:
        print "Error: "+filename+" does not exist."
        raise
    line=f.readline()
    while line:
        if COMMENT_CHAR in line:
           line, comment = line.split(COMMENT_CHAR, 1)
        if "HBA0" in line:
              AntBand="HBA0"
        elif "HBA1" in line:
              AntBand="HBA1"
        elif "HBA" in line:
              AntBand="HBA"
        elif "LBA" in line:
              AntBand="LBA"
        else:
              line=f.readline()
              continue
        where, rest = line.split(AntBand, 1)
        where=where.strip()
        if where == '':
            #Read absolute position of station origin
            line = f.readline()
            elementposshape, elementposLine = line.split(' ',1)
            elementposLine = elementposLine.strip('[] \n').split()
            position = [float(v) for v in elementposLine]
            AntFldData[AntBand]['POSITION']=position
            if AntBand!="HBA0"  and AntBand!="HBA1":
            #Read relative position of each element
                line = f.readline()
                dimstr,rest = line.split('[',1)
                shp = dimstr.split('x')
                if shp[0][0] == '(':
                    (idxbeg,idxend) = shp[0].strip('() ').split(',')
                    nrrows = int(idxend)+1
                else:
                    nrrows = int(shp[0])
                for elementNr in range(0,nrrows):
                    line = f.readline()
                    vals = line.split()
                    posxpol = [float(v) for v in vals[0:3]]
                    posypol = [float(v) for v in vals[3:6]]
                    AntFldData[AntBand]['REL_POS'].append(posxpol) #Note: Skip ypol as it is identical to xpol
                #Read ending ']' line
                line=f.readline()
        elif where=='NORMAL_VECTOR':
            line = f.readline()
            elementposshape, elementposLine = line.split(' ',1)
            elementposLine = elementposLine.strip('[] \n').split()
            nrmv = [float(v) for v in elementposLine]
            AntFldData[AntBand][where]=nrmv
        elif where=='ROTATION_MATRIX':
              line=f.readline()
              elementposshape, elementposLine = line.split(' ',1)
              elementposLine = elementposLine.strip('[] \n').split()
              elementposLine=line.split()
              dimstr,rest=line.split('[',1)
              shp=dimstr.split('x')
              for xyz in range(3):
                  line=f.readline()
                  rowstr=line.split()
                  row=[float(v) for v in rowstr]
                  AntFldData[AntBand][where].append(row)
              #Read ending ']' line
              line=f.readline()
        line=f.readline()
    return AntFldData


def parseiHBADeltasfile(stationName):
    """Parse iHBADelta file."""
    iHBADeltasdata=[]
    filepath = _getiHBADeltafile(stationName)
    f = open(filepath)
    line=f.readline()
    while line:
        if COMMENT_CHAR in line:
            line, comment = line.split(COMMENT_CHAR, 1)
        if "HBADeltas" in line:
            HBADeltaLine = True
            break
        else:
            line=f.readline()
            continue
    if HBADeltaLine:
        line=f.readline()
        
        nrelems, nrdims = line.strip().rstrip('[').split('x')
        if nrelems[0] == '(':
            nrelems = int(nrelems.strip('() ').split(',')[1])+1
            nrdims = int(nrdims.strip('() ').split(',')[1])+1
        else:
            nrelems, nrdims = int(nrelems), int(nrdims)
        for elemnr in range(nrelems):
            line=f.readline()
            xpos, ypos, zpos = line.lstrip().split()
            iHBADeltasdata.append([float(xpos), float(ypos), float(zpos)])
            #elempos[elemnr,0], elempos[elemnr,1], elempos[elemnr,2] =\
            #                 float(xpos), float(ypos), float(zpos)
    else:
        print "Error: iHBADeltas file is corrupt."
        raise
    return iHBADeltasdata


def getHBAsepton(stnName, HBAsingleelems):
    """Get configuration parameters for an HBA station in SEPTON mode. SEPTON stands for
    single element per tile ON. In this mode the HBA can be seen as an instantenously
    omni-directional interferometer capable of producing allsky snapshots."""
    stnPos, stnRot, stnRelPos, stnIntilePos = getArrayBandParams(stnName, 'HBA')
    nrtiles = len(HBAsingleelems)
    relpospertile = np.zeros((nrtiles,3))
    for tilenr in range(nrtiles):
        relpospertile[tilenr] = stnIntilePos[HBAsingleelems[tilenr]]
    stnRelPos += relpospertile
    return stnPos, stnRot, stnRelPos


def getArrayBandParams(stnName, ArrBand):
    """Get configuration parameters for an array band of the station. Array band
    can be HBA or LBA."""
    AntFld=parseAntennaField(stnName)
    stnLoc=stnName[0:2]
    if ArrBand == 'LBA':
        AntBand = 'LBA'
        hbadeltas = [0.,0.,0.]
    elif ArrBand == 'HBA':
        if stnLoc=='CS' or stnLoc=='RS':
            AntBand = 'HBA0'
        else:
            AntBand = 'HBA'
        hbadeltas = parseiHBADeltasfile(stnName)
    else:
        print("Error Array band not known. Only 'LBA' or 'HBA' are valid.")
        exit(1)
    stnPos=np.matrix(AntFld[AntBand]['POSITION']).T
    stnRot=np.matrix(AntFld[AntBand]['ROTATION_MATRIX'])
    stnRelPos=np.matrix(AntFld[AntBand]['REL_POS'])
    stnIntilePos = np.matrix(hbadeltas)

    return stnPos, stnRot, stnRelPos, stnIntilePos


def list_stations(antenna_field_dir=ANTENNAFIELDDIR):
    """List all the available LOFAR station-ids."""
    dirlist = os.listdir(antenna_field_dir)
    stnId_list = []
    for f in dirlist:
        splitres = f.split("-AntennaField.conf")
        if len(splitres) == 2 and splitres[1] == '':
            stnId = splitres[0]
            stnId_list.append(stnId)
    return stnId_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("antenna_field_file")
    args = parser.parse_args()
    AFD=parseAntennaFieldFile(args.antenna_field_file)
    print(AFD)


if __name__ == '__main__':
    main()

