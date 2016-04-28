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
ANTENNAFIELDDIR=dirname(__file__)+'/../share/AntennaFields/'
COMMENT_CHAR = '#'

def _getAntennaFieldFile(stationName, antenna_field_dir=ANTENNAFIELDDIR,
                         AFfileNameType=2):
    if AFfileNameType==2:
       basename=stationName+'-'+'AntennaField'+'.conf'
    else:
       basename='AntennaField'+stationName+'.conf'
    filepath=antenna_field_dir+'/'+basename
    return filepath


def parseAntennaField(stationName):
    filepath=_getAntennaFieldFile(stationName)
    return parseAntennaFieldFile(filepath)
    

def parseAntennaFieldFile(filename):
    AntFldData={'LBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA0': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA1': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]}
               }
    f = open(filename)
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
            line=f.readline()
            elementposLine=line.split()[2:5]
            position=[float(v) for v in elementposLine]
            AntFldData[AntBand]['POSITION']=position
            if AntBand!="HBA0"  and AntBand!="HBA1":
            #Read relative position of each element
                line=f.readline()
                dimstr,rest=line.split('[',1)
                shp=dimstr.split('x')
                for elementNr in range(0,int(shp[0])):
                  line=f.readline()
                  vals=line.split()
                  posxpol=[float(v) for v in vals[0:3]]
                  posypol=[float(v) for v in vals[3:6]]
                  AntFldData[AntBand]['REL_POS'].append(posxpol) #Note: Skip ypol as it is identical to xpol
                #Read ending ']' line
                line=f.readline()
        elif where=='NORMAL_VECTOR':
              line=f.readline()
              elementposLine=line.split()[2:5]
              nrmv=[float(v) for v in elementposLine]
              AntFldData[AntBand][where]=nrmv
        elif where=='ROTATION_MATRIX':
              line=f.readline()
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


def getArrayBandParams(stnName,ArrBand):
   AntFld=parseAntennaField(stnName)
   stnLoc=stnName[0:2]
   if ArrBand=='LBA':
      AntBand='LBA'
   elif ArrBand=='HBA':
      if stnLoc=='CS' or stnLoc=='RS':
         AntBand='HBA0'
      else:
         AntBand='HBA'
   else:
     print("Error Array band not known. Only 'LBA' or 'HBA' are valid.")
     exit(1)
   stnPos=np.matrix(AntFld[AntBand]['POSITION']).T
   stnRot=np.matrix(AntFld[AntBand]['ROTATION_MATRIX'])
   stnRelPos=np.matrix(AntFld[AntBand]['REL_POS'])
   return stnPos, stnRot, stnRelPos


def list_stations(antenna_field_dir=ANTENNAFIELDDIR):
    """List all the available LOFAR station-ids."""
    dirlist = os.listdir(antenna_field_dir)
    stnId_list = []
    for f in dirlist:
        (stnId, blank) = f.split("-AntennaField.conf",1)
        if blank is "":
            stnId_list.append(stnId)
    return stnId_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("antenna_field_file")
    args = parser.parse_args()
    AFD=parseAntennaFieldFile(args.antenna_field_file)
    print(AFD)
