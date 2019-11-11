#!/usr/bin/env python

### Author: Qi Feng (based on Matt Buchovecky's sourceList.py)
# script that takes an optional argument for the date and target collection and calculates angular separation and elevation of each target from the moon.

import sys, operator, argparse
import ephem, subprocess
from math import ceil
from astropy.coordinates import SkyCoord
from astropy import units as u




# setting up ephem observer object for veritas
veritas = ephem.Observer()
veritas.lat = '31:40.51'
veritas.lon = '-110:57.132'
veritas.elevation = 1268


# argument parser
parser = argparse.ArgumentParser(
    description="Takes optional arguments to specify date and source collection, and min / max moon distances. If no arguments are specified, will choose from all psf sources to make an ordered list appropriate for taking a PSF measurement, as well as suggested exposure times for each source. For a more general look at elevation and moon distance, see moonDist.py..")

parser.add_argument('--date', default=veritas.date,
                    help="Specify DATE (in UT) in the format \"YYYY/MM/DD HH:MM:SS\"   don\'t forget the quotation marks. The default value is today's date.")

parser.add_argument('--minMoonDist', default=0., type=float,
                    help="The minimum distance in degrees that a source should be from the moon to include it in the list. The default value is 30 degrees.")

parser.add_argument('--maxMoonDist', default=180., type=float,
                    help="The maximum distance in degrees that a source should be from the moon, to prevent backlighting and arm shadows. The default value is 90 degrees.")

parser.add_argument('--minElevation', default=20., type=float,
                    help="The minimum elevation in degrees you would like to look at. The default is 20 degrees.")

parser.add_argument('--maxElevation', default=90., type=float,
                    help="The minimum elevation in degrees you would like to look at. The default is 90 degrees.")

parser.add_argument('-s', '--star', default="Vega",
                    help="The name of a star. Default is Vega. ")

args = parser.parse_args()

# setting date/time to user-spefied value (or default to current date/time)
veritas.date = args.date
# letting user know the date and target collection used.
print
print "Date and time used (in UT): %s" % veritas.date

moonlightSources = {}

try:
    aRADec = SkyCoord.from_name(args.star, frame='icrs')
    realRADec = aRADec.to_string('hmsdms')
    sourceRA = realRADec.split()[0]
    sourceDEC = realRADec.split()[1]
    sourceRA_deg = aRADec.ra.degree
    sourceDec_deg = aRADec.dec.degree
    sourceRA_rad = aRADec.ra.radian
    sourceDec_rad = aRADec.dec.radian
except:
    print("Cannot find coordinates for name {}".format(args.star))
    exit(1)

sourceEpoch = 2000

# makes sure same epoch is used
veritas.epoch = float(sourceEpoch)

# Define ephem moon object and calculate position (ra, dec) and phase
TheMoon = ephem.Moon(veritas)
TheMoon.compute(veritas)
illum = TheMoon.moon_phase * 100.

moonReflection = ephem.FixedBody(TheMoon.ra, TheMoon.dec)

# Get angular separation of moon and target
degFromMoon = 180. / ephem.pi * ephem.separation((TheMoon.ra, TheMoon.dec), (sourceRA_rad, sourceDec_rad))

# Define ehpem object for source, to get elevation
sourceObj = ephem.FixedBody()
sourceObj._ra = sourceRA_rad
sourceObj._dec = sourceDec_rad
sourceObj.compute(veritas)


sourceEl = sourceObj.alt * 180. / ephem.pi  # elevation of source
sourceAz = sourceObj.az * 180. / ephem.pi  # azimuth of source

moonlightSources[args.star]=(sourceEl, sourceAz, sourceRA, sourceDEC, degFromMoon)


# for first column
columnTitle = "Source"
print("")
print(columnTitle),
print('\t'),
col_string = "Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist"

print(col_string)
print("-------------------------------------------------------------------------------------------------")
name = args.star
el = sourceEl  # el
az = sourceAz # azimuth
ra = sourceRA  # ra
dec = sourceDEC  # dec
dist = degFromMoon

print(name),
print("\t"),
row_string = "%0.2f\t\t%0.2f\t\t%s\t%s\t%0.2f" %(el, az, ra, dec, dist)
print(row_string)

print("-------------------------------------------------------------------------------------------------")
print("The Moon is %0.2f%% illuminated" % illum)
print(TheMoon.dec)


sys.exit(0)  # great job

