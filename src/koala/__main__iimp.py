# -*- coding: utf-8 -*-


def obtain_throughput_correction_iimp(data_sci):
    # To obtain the throughput_correction we need to reduce the skyflat
    if not data_sci.skyflat.reduced():
        reduce_science_images_iimp(skyflat=True)
    generate_throughput_correctoin_iimp()


def obtain_flux_calibration_iimp(data_cal):
    # First we check the input files to see if they are already reduced.
    for calib_star in data_cal.calib_star:
        if not calib_star.is_reduced_data:
            # Generate reduced data
            reduce_science_images_iimp(skyflat=False)
        # Now check if the reduced file exists.
        if not calib_star.abs_flux_cal:
            # Generated the absolute flux calibration.
            generate_flux_calibration_iimp()
    # We can now obtain the flux calibration file.
    generate_combined_flux_cal_iimp()  # flux_calibration_file is now created and placed in the .yml file


def obtain_telluric_correction_iimp(data_cal):
    # First we check the input files to see if they are already reduced.
    for calib_star in data_cal.calib_star:
        if not calib_star.is_reduced_data:
            # Generate reduced data
            reduce_science_images_iimp(skyflat=False)
        # Now check if the reduced file exists.
        if not calib_star.telluric_cal:
            # Generated the telluric calibration.
            generate_tel_calibration_iimp()
    # We can now obtain the flux calibration file.
    generate_combined_tel_cal_iimp()  # flux_calibration_file is now created and placed in the .yml file


def reduce_science_images_iimp(skyflat):
    # What I want this function to do is run KOALA_reduce, but we don't want it to return anything
    # What we want it to do is update the internal file.
    pass


def generate_flux_calibration_iimp():
    # THis function generates the absolute flux cal. can save it as a file. But also
    # Updates the config to have absolute_flux_cal for this star.
    pass


def generate_tel_calibration_iimp():
    # THis function generates the telluric correction file for each calibration star.
    # updates the .yml file
    pass


def generate_combined_flux_cal_iimp():
    # This takes in all the absolute flux calibration data and generates the file needed
    # for reduction.
    pass


def generate_combined_tel_cal_iimp():
    # THis takes in all the telluric calibration files and geneates the main file needed. updates the .yml
    pass


def generate_throughput_correctoin_iimp():
    # THis takes in the REDUCED skyflat and then generates the throughput correction, have this update
    # the .yml and everything.
    pass


def obtain_sky_iimp():
    # This is the function that takes in the science files. and reduces obtains the sky from them.
    # TODO: might need to check if we have to reduce the data.
    pass
