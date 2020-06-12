# -*- coding: utf-8 -*-

# APPROX TIME without doing the actual koala reduction
# BOTH FILES ~ 1 s
# No FLUX: 167s, 149s, 155s
# No Tell: 165s, 163s,
# neither file: 172s, 184s
# just doing reduce: 163s,  176s
# From quick testing koala_reduce is 80s per object, telluric is like 7s and flux is like 4s

from timeit import default_timer as timer

import os.path as pth
import os

os.chdir("..")
os.chdir("..")

import matplotlib.pyplot as plt

plt.ion()
DO_PLOTTING = False

from koala.utils.utils import SciFrame
from koala.loader import load_science_images, load_calibaration_files
from koala.__main__functions import (
    obtain_throughput_correction_blake,
    reduce_science_images,
    create_reduction_files,
    get_sky,obtain_flux_calibration_blake, obtain_telluric_correction_blake
)

from koala.__main__ import (
    DATA_PATH,
    version,
)

DATE = "20180310"
GRATING = "385R"
PIXEL_SIZE = 0.6  # Just 0.1 precision
KERNEL_SIZE = 1.25
OBJECT = "POX4"
DESCRIPTION = "POX4 CUBE"

data_sci = load_science_images("config_science_images.yml")
data_cal = load_calibaration_files("config_calibration_stars.yml")

FILE_SKY_FLAT_RED = data_sci.sci_calib_data.skyflat
# FILE NOT DIVIDED BY THE FLAT
THROUGHPUT_FILE_RED = data_sci.sci_calib_data.throughput_correction
FLUX_CAL_FILE = data_sci.sci_calib_data.flux_calibration
TELLURIC_CORRECTION_FILE = data_sci.sci_calib_data.telluric_correction
"""
so for telluric correction is appears to be like, xxx_xxx_xxx_telluric_correction.dat is the correction FROM a single star
and then the telluric_correction_xx_xxx_xxx.dat is the combined telluric correction object thing. it is just a 
nanmedian of all the telluric corrections given. 
"""
SCIENCE_RED_1_FILENAME = data_sci.sci_rss[SciFrame.frame1]
SCIENCE_RED_2_FILENAME = data_sci.sci_rss[SciFrame.frame2]
SCIENCE_RED_3_FILENAME = data_sci.sci_rss[SciFrame.frame3]
SCIENCE_RED_FILENAMES = [
    SCIENCE_RED_1_FILENAME,
    SCIENCE_RED_2_FILENAME,
    SCIENCE_RED_3_FILENAME,
]

start = timer()


if __name__ == "__main__":

    print("\n> Running Pykoala. Version:", version)

    pk = "_{}p{}_{}k{}".format(
        int(PIXEL_SIZE),
        int((abs(PIXEL_SIZE) - abs(int(PIXEL_SIZE))) * 10),
        int(KERNEL_SIZE),
        int((abs(KERNEL_SIZE) - abs(int(KERNEL_SIZE))) * 100),
    )  # new name? or at least what does pk stand for?

    # ----------------------------------------------------------------------------------------------------------------------
    # THROUGHPUT CORRECTION USING SKYFLAT
    # ----------------------------------------------------------------------------------------------------------------------
    # The very first thing that we need is to get the throughput correction.
    # IMPORTANT: We use a skyflat that has not been divided by a flatfield in 2dFdr !!!!!!
    # If this has been done before, we can read the file containing the throughput correction
    # throughput_red = read_table(THROUGHPUT_FILE_RED, ["f"] )
    # TODO: throughput_red isn't used anywhere maybe at some point you could feed that in instead of the skyflat
    # TODO: There is no code which allows you to use this file that you read. Maybe at some point it was just used for
    # TODO: plotting purposes. but currently and in the monolithic file there is no mention of how to use it.

    # Now we read the RSS file, we ONLY correct for ccd defects and high cosmics
    skyflat_red = obtain_throughput_correction_blake(FILE_SKY_FLAT_RED, THROUGHPUT_FILE_RED)

    # ############################## ############################## ############################## ############################## #

    if not data_sci.sci_calib_data.reduction_files_exist:  # The files don't exist and we need to create them
        flux_calibration, telluric_correction = create_reduction_files(skyflat_red, pk, data_cal, data_sci)
    else:
        w_star, flux_calibration = obtain_flux_calibration_blake(FLUX_CAL_FILE)
        w_star, telluric_correction = obtain_telluric_correction_blake(TELLURIC_CORRECTION_FILE)

    # ############################## ############################## ##############################  ############################## #
    # ---------------------------------------------------------------------------
    #  OBTAIN SKY SPECTRA IF NEEDED
    # ---------------------------------------------------------------------------
    # Using the same files than objects but choosing fibres without object emission
    sky_r = []
    sky = []
    for file in SCIENCE_RED_FILENAMES:
        sky_from_in, sky_r_from_in = get_sky(file, skyflat_red)
        sky.append(sky_from_in)
        sky_r.append(sky_r_from_in)


# ------------------------------------------------------------------------------------------------------------------
# TIME FOR THE OBJECT !!
# ------------------------------------------------------------------------------------------------------------------
rss_list = [
    SCIENCE_RED_1_FILENAME,
    SCIENCE_RED_2_FILENAME,
    SCIENCE_RED_3_FILENAME,
]

sky_list = sky  

fits_file_red = pth.join(DATA_PATH, GRATING, "POX4_A_red_combined_cube_2_TEST_GitHub.fits")
reduced_hikids = reduce_science_images(
    rss_list, fits_file_red, skyflat_red, sky_list, telluric_correction, flux_calibration, OBJECT, DESCRIPTION
)

end = timer()
print("\n> Elapsing time = ", end - start, "s")
