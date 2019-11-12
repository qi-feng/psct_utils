#!/usr/bin/env python

### Author: Qi Feng (based on Matt Buchovecky's sourceList.py)
# script that takes an optional argument for the date and target collection and calculates angular separation and elevation of each target from the moon.

import sys, operator, argparse
import ephem, subprocess
from math import ceil
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np



# setting up ephem observer object for veritas
veritas = ephem.Observer()
veritas.lat = '31:40.51'
veritas.lon = '-110:57.132'
veritas.elevation = 1268


def calc_source(date, star, wobble_offset=0, wob_az_el=False):
    # setting date/time to user-spefied value (or default to current date/time)
    veritas.date = date
    # letting user know the date and target collection used.
    print("")
    print("Date and time used (in UT): %s" % date)


    try:
        aRADec = SkyCoord.from_name(star, frame='icrs')
        realRADec = aRADec.to_string('hmsdms')
        sourceRA = realRADec.split()[0]
        sourceDEC = realRADec.split()[1]
        sourceRA_deg = aRADec.ra.degree
        sourceDec_deg = aRADec.dec.degree
        sourceRA_rad = aRADec.ra.radian
        sourceDec_rad = aRADec.dec.radian
    except:
        print("Cannot find coordinates for name {}".format(star))
        exit(1)

    el, az = calc_from_ra_dec(veritas, sourceRA_rad, sourceDec_rad)
    print_info(veritas, star, el, az, sourceRA, sourceDEC,
               col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
               columnTitle="Source")

    if wobble_offset > 0 and wob_az_el: # this is only useful for optical alignment stuff
        print("Calculate wobble with offset {} deg...".format(wobble_offset))
        columnTitle = "Source"
        print("-------------------------------------------------------------------------------------------------")
        columnLength = ceil(float(len("{}+EL{}deg".format(star, wobble_offset))) / 8.) * 8.
        columnTabs = int(ceil((columnLength - len(columnTitle)) / 8.))
        print(
            "-------------------------------------------------------------------------------------------------")
        # sys.stdout.write( columnTitle )
        print(columnTitle),
        for x in range(0, columnTabs):
            print('\t'),
        print("Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist")
        print("-------------------------------------------------------------------------------------------------")
        #EL up
        el_wob = el + wobble_offset
        az_wob = az
        ra, dec = veritas.radec_of(az_wob/180.*ephem.pi, el_wob/180.*ephem.pi)
        #c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
        coord = SkyCoord(ra, dec, frame='icrs', unit=u.radian)
        el_, az_ = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

        #print(el_, az_, el_wob, az_wob)

        realRADec = coord.to_string('hmsdms')
        sourceRA = realRADec.split()[0]
        sourceDEC = realRADec.split()[1]

        print_info(veritas,"{}+EL{}deg".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
               col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
               columnTitle="Source", concise=True)

        #EL down
        el_wob = el - wobble_offset
        az_wob = az
        ra, dec = veritas.radec_of(az_wob / 180. * ephem.pi, el_wob / 180. * ephem.pi)
        # c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
        coord = SkyCoord(ra, dec, frame='icrs', unit=u.radian)
        el_, az_ = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

        #print(el_, az_, el_wob, az_wob)

        realRADec = coord.to_string('hmsdms')
        sourceRA = realRADec.split()[0]
        sourceDEC = realRADec.split()[1]

        print_info(veritas, "{}-EL{}deg".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
                   col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
                   columnTitle="Source", concise=True)

        #AZ right
        el_wob = el
        az_wob = az + wobble_offset/np.cos(el_wob/ 180. * ephem.pi)

        ra, dec = veritas.radec_of(az_wob / 180. * ephem.pi, el_wob / 180. * ephem.pi)
        # c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
        coord = SkyCoord(ra, dec, frame='icrs', unit=u.radian)
        el_, az_ = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

        #print(el_, az_, el_wob, az_wob)

        realRADec = coord.to_string('hmsdms')
        sourceRA = realRADec.split()[0]
        sourceDEC = realRADec.split()[1]

        print_info(veritas, "{}+AZ{}deg".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
                   col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
                   columnTitle="Source", concise=True)


        #AZ left
        el_wob = el
        az_wob = az - wobble_offset/np.cos(el_wob/ 180. * ephem.pi)

        ra, dec = veritas.radec_of(az_wob / 180. * ephem.pi, el_wob / 180. * ephem.pi)
        # c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
        coord = SkyCoord(ra, dec, frame='icrs', unit=u.radian)
        el_, az_ = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

        #print(el_, az_, el_wob, az_wob)

        realRADec = coord.to_string('hmsdms')
        sourceRA = realRADec.split()[0]
        sourceDEC = realRADec.split()[1]

        print_info(veritas, "{}-AZ{}deg".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
                   col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
                   columnTitle="Source", concise=True)
        print("-------------------------------------------------------------------------------------------------")

    if wobble_offset > 0 and not wob_az_el:  # this is the real wobble in RA Dec
            print("Calculate wobble with offset IN RA and Dec by {} deg...".format(wobble_offset))
            print("-------------------------------------------------------------------------------------------------")
            columnTitle = "Source"
            columnLength = ceil(float(len("{}-{}N".format(star, wobble_offset))) / 8.) * 8.
            columnTabs = int(ceil((columnLength - len(columnTitle)) / 8.))
            print(
                "-------------------------------------------------------------------------------------------------")
            # sys.stdout.write( columnTitle )
            print(columnTitle),
            for x in range(0, columnTabs):
                print('\t'),
            print("Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist")
            print("-------------------------------------------------------------------------------------------------")

            # N
            coord = SkyCoord(sourceRA_rad, sourceDec_rad+wobble_offset/180*ephem.pi, frame='icrs', unit=u.radian)
            el_wob, az_wob = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

            realRADec = coord.to_string('hmsdms')
            sourceRA = realRADec.split()[0]
            sourceDEC = realRADec.split()[1]

            print_info(veritas, "{}-{}N".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
                       col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
                       columnTitle="Source", concise=True)

            # S
            coord = SkyCoord(sourceRA_rad, sourceDec_rad-wobble_offset/180*ephem.pi, frame='icrs', unit=u.radian)
            el_wob, az_wob = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

            realRADec = coord.to_string('hmsdms')
            sourceRA = realRADec.split()[0]
            sourceDEC = realRADec.split()[1]

            print_info(veritas, "{}-{}S".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
                       col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
                       columnTitle="Source", concise=True)

            # E
            coord = SkyCoord(sourceRA_rad +  wobble_offset / 180 * ephem.pi / np.cos(sourceDec_rad), sourceDec_rad, frame='icrs', unit=u.radian)
            el_wob, az_wob = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

            realRADec = coord.to_string('hmsdms')
            sourceRA = realRADec.split()[0]
            sourceDEC = realRADec.split()[1]

            print_info(veritas, "{}-{}E".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
                       col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
                       columnTitle="Source", concise=True)

            # W
            coord = SkyCoord(sourceRA_rad - wobble_offset / 180 * ephem.pi / np.cos(sourceDec_rad), sourceDec_rad,
                             frame='icrs', unit=u.radian)
            el_wob, az_wob = calc_from_ra_dec(veritas, coord.ra.radian, coord.dec.radian)

            realRADec = coord.to_string('hmsdms')
            sourceRA = realRADec.split()[0]
            sourceDEC = realRADec.split()[1]

            print_info(veritas, "{}-{}W".format(star, wobble_offset), el_wob, az_wob, sourceRA, sourceDEC,
                       col_string="Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
                       columnTitle="Source", concise=True)

            print("-------------------------------------------------------------------------------------------------")
    print("")

def calc_from_ra_dec(veritas, sourceRA_rad, sourceDec_rad, sourceEpoch = 2000):

    # makes sure same epoch is used
    veritas.epoch = float(sourceEpoch)

    # Define ehpem object for source, to get elevation
    sourceObj = ephem.FixedBody()
    sourceObj._ra = sourceRA_rad
    sourceObj._dec = sourceDec_rad
    sourceObj.compute(veritas)


    sourceEl = sourceObj.alt * 180. / ephem.pi  # elevation of source
    sourceAz = sourceObj.az * 180. / ephem.pi  # azimuth of source


    return sourceEl, sourceAz


def print_info(veritas, star, sourceEl, sourceAz, sourceRA, sourceDEC,
               col_string = "Elevation\tAzimuth\t\tRA\t\tDec\t\tMoonDist",
               columnTitle = "Source", concise=False):

    # for first column

    columnLength = ceil(float(len(star))/8.)*8.
    columnTabs = int( ceil( (columnLength-len(columnTitle))/8.) )
    if not concise:
        print("-------------------------------------------------------------------------------------------------")
        #sys.stdout.write( columnTitle )
        print(columnTitle),
        for x in range(0, columnTabs):
            print('\t'),

    # Define ephem moon object and calculate position (ra, dec) and phase
    TheMoon = ephem.Moon(veritas)
    TheMoon.compute(veritas)

    illum = TheMoon.moon_phase * 100.

    moonReflection = ephem.FixedBody(TheMoon.ra, TheMoon.dec)

    # Get angular separation of moon and target
    source_coord = SkyCoord(sourceRA, sourceDEC, frame='icrs')
    sourceRA_rad = source_coord.ra.radian
    sourceDec_rad = source_coord.dec.radian
    degFromMoon = 180. / ephem.pi * ephem.separation((TheMoon.ra, TheMoon.dec), (sourceRA_rad, sourceDec_rad))


    name = star
    el = sourceEl  # el
    az = sourceAz # azimuth
    ra = sourceRA  # ra
    dec = sourceDEC  # dec
    dist = degFromMoon

    length = len(name)
    numTabs = int( ceil( (columnLength-length-1)/8. ) )
    if not concise:
        print(col_string)
        print("-------------------------------------------------------------------------------------------------")

    print(name),
    for i in range (0, numTabs):
        print("\t"),

    row_string = "%0.2f\t\t%0.2f\t\t%s\t%s\t%0.2f" %(el, az, ra, dec, dist)
    print(row_string)

    if not concise:
        print("-------------------------------------------------------------------------------------------------")
        print("The Moon is %0.2f%% illuminated" % illum)
        print(TheMoon.dec)



if __name__ == "__main__":
    # argument parser
    parser = argparse.ArgumentParser(
        description="Takes optional arguments to specify date and source collection, and min / max moon distances. If no arguments are specified, will choose from all psf sources to make an ordered list appropriate for taking a PSF measurement, as well as suggested exposure times for each source. For a more general look at elevation and moon distance, see moonDist.py..")

    parser.add_argument('-d', '--date', default=veritas.date,
                        help="Specify DATE (in UT) in the format \"YYYY/MM/DD HH:MM:SS\"   don\'t forget the quotation marks. The default value is today's date.")

    parser.add_argument('-s', '--star', default="Vega",
                        help="The name of a star. Default is Vega. ")

    parser.add_argument('-w', '--wobble', default=0.5, type=float,
                        help="Wobble offset in deg. Default is 0.5 deg. ")

    parser.add_argument('-t', '--wobble_az_el', action='store_true',
                        help="Wobble in AZ EL instead of RA Dec, non-conventional. Default is False. ")

    args = parser.parse_args()

    calc_source(args.date, args.star, wobble_offset=args.wobble, wob_az_el=args.wobble_az_el)

    sys.exit(0)  # great job

