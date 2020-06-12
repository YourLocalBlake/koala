from koala.utils.io import read_table, array_to_text_file, spectrum_to_text_file
from koala.__main__ import DATA_PATH, KOALA_RSS, KOALA_reduce
from koala.utils.flux import obtain_flux_calibration
from koala.utils.plots import plot_response
from koala.utils.spectrum_tools import obtain_telluric_correction

DO_PLOTTING = False

DATE = "20180310"
GRATING = "385R"
PIXEL_SIZE = 0.6  # Just 0.1 precision
KERNEL_SIZE = 1.25
OBJECT = "POX4"
DESCRIPTION = "POX4 CUBE"


def obtain_flux_calibration_blake(file_loc):
    """ Checks if flux_calibration file exists, else signal the file to be created

    Check to see if the flux_cal
    ibration file location specified in config_science_images.yml exists. If it
    does, returns the file. Else the function will attempt to created a combined flux_calibration file for
    data reduction.

    Attributes
    -------
    file_loc: str
        location of the star file to calculate the flux_calibration of

    Returns
    -------
    str
        location of the combined flux calibration file.
    """
    if file_loc is None:  # Then the flux_calibration file does not exist.
        # I think this is hard to create at the time being due to how KOALA_reduce works.
        # As in I would have to rerun several bits of code a lot.
        pass
    # What this should do it create it, then add it to the config_science_images file.

    # AS A TEMP SOLUTION TILL TOP IS IMPLEMENTED TO REMOVE CODE FROM __main__ I am doing
    wavelength, flux_calibration = read_table(file_loc, ["f", "f"])
    return wavelength, flux_calibration


def obtain_throughput_correction_blake(
    skyflat_file_loc, throughput_correction_file_loc
):
    """ Obtain the relative_throughput of the skyflat

        Calculates the relative_thoughput of the skyflat file. If the .yml config file is missing the
        throughput_correction .dat file adds it to the file (TODO).
        Figure out if we want to actually read in the throughput_correction at some point...
        ToDO: figure out best way to save to the .yml file.

        Attributes
        -------
        skyflat_file_loc: str
            location of the skyflat file.
        throughput_correction_file_loc: str
            location of the already calculated throughput_correction.dat file - Can be None

        Returns
        -------
        RSSObject
            The same skyflat that was passed in, except with throughput_correction calculated
        """
    skyflat_red = KOALA_RSS(
        skyflat_file_loc,
        flat="",
        apply_throughput=False,
        sky_method="none",
        # skyflat = skyflat_red,  # TODO: what is this for? this should be here cause it is the flat...
        do_extinction=False,
        correct_ccd_defects=False,
        correct_high_cosmics=False,
        clip_high=100,
        step_ccd=50,
        plot=DO_PLOTTING,
    )

    # If the data have been normalized by the FLATFIELD, we only need a SCALE between fibres, We consider the median
    # value in the range [wave_min_scale, wave_max_scale] for all fibres and scale # TODO: figure out what this mean
    skyflat_red.find_relative_throughput(
        ymin=0, ymax=800000, wave_min_scale=6300, wave_max_scale=6500, plot=DO_PLOTTING,
    )
    # The relative throughput is an array stored in skyflat_red.relative_throughput

    # TODO: need to do the whole if it's None then we need a good default file name.
    if throughput_correction_file_loc is None:
        filename_default = DATA_PATH + "{}_{}_throughput_correction.dat".format(
            DATE, GRATING
        )
        array_to_text_file(skyflat_red.relative_throughput, filename=filename_default)

    return skyflat_red


def obtain_telluric_correction_blake(file_loc):
    """ Obtain the telluric correction for a given star file.

        Check to see if the telluric_correction file location specified in config_science_images.yml exists. If it
        does, returns the file. Else the function will attempt to created a combined telluric_correction file for
        data reduction.

        Attributes
        -------
        file_loc: str
            location of the star file to calculate the telluric_correction of

        Returns
        -------
        str
            location of the combined telluric correction file.
        """
    if file_loc is None:  # Then the telluric_correction file does not exist.
        # I think this is hard to create at the time being due to how KOALA_reduce works.
        # As in I would have to rerun several bits of code a lot.
        pass
    # What this should do it create it, then add it to the config_science_images file.

    wavelength, telluric_correction = read_table(file_loc, ["f", "f"])
    return wavelength, telluric_correction


def reduce_koala_blake(star, skyflat_red, pk):
    """ temp holder for reducing koala, star is a calib_star_gen object... I think 10/6/2020

    Parameters
    ----------
    star: gencalib star object, probs.

    Returns
    -------

    """
    # todo: fits_file_Red here should prob become like the name of something.
    fits_file_red = DATA_PATH + GRATING + star.star_name + "_used_in_cal"  # THIS IS THE FILENAME
    positions = star.star_files  # Star position 1 through to i.e. 3.
    star_red = KOALA_reduce(  # TODO: blake, red for reduced cause I'm laz.y
        positions,
        fits_file=fits_file_red + ".fits",
        obj_name=star.star_name,
        description=star.star_name,
        apply_throughput=True,
        skyflat=skyflat_red,
        correct_ccd_defects=True,
        correct_high_cosmics=True,
        clip_high=50,
        step_ccd=50,
        sky_method="self",
        n_sky=400,
        pixel_size_arcsec=PIXEL_SIZE,
        kernel_size_arcsec=KERNEL_SIZE,
        ADR=False,
        valid_wave_min=6085,
        valid_wave_max=9305,
        plot=DO_PLOTTING,
        warnings=False,
    )
    star_red.combined_cube.half_light_spectrum(r_max=5, plot=DO_PLOTTING)
    intergrated_star_flux_file_name = "{}/{}/CALIBRATION_STARS/{}_{}{}_integrated_star_flux.dat".format(
        DATA_PATH, GRATING, star.star_name, GRATING, pk
    )
    spectrum_to_text_file(
        star_red.combined_cube.wavelength,
        star_red.combined_cube.integrated_star_flux,
        filename=intergrated_star_flux_file_name,
    )
    return star_red


def preform_telluric_correction(star_red, star_name, pk):
    telluric_correction_star1 = star_red.get_telluric_correction(
        n_fibres=15,
        correct_from=6830.0,
        correct_to=8400.0,
        exclude_wlm=[
            [6000, 6350],
            [6460, 6720],
            [6830, 7450],
            [7550, 7750],
            [8050, 8400],
        ],
        apply_tc=True,
        combined_cube=True,
        weight_fit_median=1.0,
        step=15,
        wave_min=6085,
        wave_max=9305,
    )
    telluric_file = "{}/{}/CALIBRATION_STARS/{}_{}{}_telluric_correction.dat".format(
        DATA_PATH, GRATING, star_name, GRATING, pk
    )
    # We can save this calibration as a text file
    spectrum_to_text_file(
        star_red.combined_cube.wavelength,
        telluric_correction_star1,
        filename=telluric_file,
    )
    return telluric_correction_star1


def preform_flux_calibration(star_red, star_name, pk):
    # ----------------------------------------------------------------------------------------------------------------------
    # -------------------   FLUX CALIBRATION (COMMON) -----------
    # Now we read the absolute flux calibration data of the calibration star and get the response curve
    # (Response curve: correspondence between counts and physical values)
    # Include exp_time of the calibration star, as the results are given per second
    # For this BE CAREFUL WITH ABSORPTIONS (Halpha) and check behaviour in the edges of the CCD
    # Change fit_degree (3,5,7), step, min_wave, max_wave to get better fits !!!

    # BLAKE: so we are doing response curves now.
    star_red.combined_cube.do_response_curve(
        DATA_PATH
        + "/FLUX_CAL/fhilt600_edited.dat",  # THIS IS HARD CODED. this will need to be another config
        # file, to get e.g. what stars we want to use for response_curves.
        plot=DO_PLOTTING,
        min_wave=6110.0,
        max_wave=9305.0,
        step=20,
        exp_time=120.0,
        fit_degree=7,
    )
    response_file_red = "{}/{}/CALIBRATION_STARS/{}_{}{}_response.dat".format(
        DATA_PATH, GRATING, star_name, GRATING, pk
    )
    # Now we can save this calibration as a text file
    spectrum_to_text_file(
        star_red.combined_cube.response_wavelength,
        star_red.combined_cube.response_curve,
        filename=response_file_red,
    )
    return star_red


def obtained_combined_flux_cal(reduced_stars_list, pk):
    # Define in "stars" the 2 cubes we are using, and plotting their responses to check
    stars = [star.combined_cube for star in reduced_stars_list]
    plot_response(stars, scale=[1, 1.14, 1.48, 1])
    # We obtain the flux calibration applying:
    flux_calibration = obtain_flux_calibration(stars)

    # And we save this absolute flux calibration as a text file
    flux_calibration_file = (
            DATA_PATH + "/flux_calibration_" + DATE + "_" + GRATING + pk + ".dat"
    )
    spectrum_to_text_file(
        reduced_stars_list[0].combined_cube.wavelength,   # THis is wavelength
        flux_calibration,  # This is the data
        filename=flux_calibration_file,
    )
    return flux_calibration, flux_calibration_file


def obtain_combined_telluric_correction(reduced_tell_list, reduce_stars, pk):
    # TODO: blake, I put in the [x] here just kind of correspond the the stars in it. I think it was specified ones
    #  Angel did each time he manually wrote a star reduction.

    # TODO: we used to have reduce_star_list_b[1].combined_cube.wavelength as the first argument to
    #  obtain_telluric_correction and specturm_to-text_file. I changed it to the new list we have in, They are just
    #  the wavelengths I believe, hence they should be interchangable. we might need the flux calibration wavelengths
    #  tho? to be determined.

    telluric_correction = obtain_telluric_correction(
        reduce_stars[1].combined_cube.wavelength, reduced_tell_list
    )
    # Save this telluric correction to a file
    telluric_correction_file = (
            DATA_PATH + "/telluric_correction_" + DATE + "_" + GRATING + pk + ".dat"
    )
    spectrum_to_text_file(
        reduce_stars[1].combined_cube.wavelength,
        telluric_correction,
        filename=telluric_correction_file,
    )
    return telluric_correction, telluric_correction_file

def reduce_science_images(rss_list, fits_file_red, skyflat_red, sky_list, telluric_correction, flux_calibration,
                          OBJECT, DESCRIPTION):
    hikids_red = KOALA_reduce(
        rss_list,
        obj_name=OBJECT,
        description=DESCRIPTION,
        # rss_clean=True,
        fits_file=fits_file_red,
        # save_rss_to_fits_file_list=save_rss_list,
        apply_throughput=True,
        skyflat=skyflat_red,
        plot_skyflat=False,
        correct_ccd_defects=True,
        correct_high_cosmics=False,
        clip_high=100,
        step_ccd=50,
        # fix_wavelengths=True,
        # sol=[0.119694453613, -0.000707644207572, 2.03806478671e-07],
        # sky_method="1Dfit",
        sky_method="1D",
        sky_list=sky_list,
        scale_sky_1D=1.0,
        auto_scale_sky=True,
        brightest_line="Ha",
        brightest_line_wavelength=6641.0,
        id_el=False,
        high_fibres=10,
        cut=1.5,
        plot_id_el=True,
        broad=1.8,
        id_list=[
            6300.30,
            6312.1,
            6363.78,
            6548.03,
            6562.82,
            6583.41,
            6678.15,
            6716.47,
            6730.85,
            7065.28,
            7135.78,
            7318.39,
            7329.66,
            8750.47,
            8862.79,
            9014.91,
            9069.0,
        ],
        # clean_sky_residuals=False,
        # dclip=3.0,
        # extra_w=1.3,
        # step_csr = 25,
        telluric_correction=telluric_correction,
        do_extinction=True,
        correct_negative_sky=False,
        pixel_size_arcsec=PIXEL_SIZE,
        kernel_size_arcsec=KERNEL_SIZE,
        # offsets=[-0.54, -0.87, 1.58, -1.26] # EAST-/WEST+  NORTH-/SOUTH+
        ADR=False,
        flux_calibration=flux_calibration,
        # size_arcsec=[60,60],
        valid_wave_min=6085,
        valid_wave_max=9305,
        plot=DO_PLOTTING,
        warnings=False,
    )
    return hikids_red


def create_reduction_files(skyflat_red, pk, data_cal, data_sci):
    # TODO: so we have data_cal.calib_stars. THIS IS A LIST OF THE CAL STARS.

    # So first we reduce the stars.
    reduced_stars = [reduce_koala_blake(star, skyflat_red, pk) for star in data_cal.calib_stars]
    # TODO: in the future I want to be able to jsut have the data_cal.calib_stars objects be changed. Like the
    #  internal rss stores the changes.

    reduce_star_list_b = []
    reduce_tel_cur_b = []

    for index, star_red in enumerate(reduced_stars):

        star_name = data_cal.calib_stars[index].star_name
        if not data_sci.sci_calib_data.telluric_correction_exists:
            reduce_tel_cur_b.append(preform_telluric_correction(star_red, star_name, pk))  # TODO: once again we
            # don't want another list.

        if not data_sci.sci_calib_data.flux_calibration_exists:
            reduce_star_list_b.append(preform_flux_calibration(star_red, star_name, pk))
        # #  CHECK AND GET THE FLUX CALIBRATION FOR THE NIGHT  RED
        # #  First we take another look to the RSS data ploting the integrated fibre values in a map
        # TODO: blake, make this some function or something.
        #    star1r.RSS_map(star1r.integrated_fibre, norm=colors.PowerNorm(gamma=1./4.))  # Dead fibre!!!
        #    star2r.RSS_map(star2r.integrated_fibre, norm=colors.PowerNorm(gamma=1./4.))
        #    star3r.RSS_map(star3r.integrated_fibre, norm=colors.PowerNorm(gamma=1./4.))

        # # We check again that star1 is on a dead fibre, we don't use this star for absolute flux calibration

    if not data_sci.sci_calib_data.flux_calibration_exists:
        flux_calibration, flux_cal_file = obtained_combined_flux_cal(reduce_star_list_b, pk)
        FLUX_CAL_FILE = flux_cal_file
    else:
        FLUX_CAL_FILE = data_sci.sci_calib_data.flux_calibration

    # Similarly, provide a list with the telluric corrections and apply:
    if not data_sci.sci_calib_data.telluric_correction_exists:
        telluric_correction, telluric_correction_file = obtain_combined_telluric_correction(
            reduce_tel_cur_b, reduced_stars, pk
        )
        TELLURIC_CORRECTION_FILE = telluric_correction_file
    else:
        TELLURIC_CORRECTION_FILE = data_sci.sci_calib_data.telluric_correction


    w_star, flux_calibration = obtain_flux_calibration_blake(FLUX_CAL_FILE)
    w_star, telluric_correction = obtain_telluric_correction_blake(TELLURIC_CORRECTION_FILE)

    # TODO: do we ever use wavelength of the star?

    # TODO: with this section we need to write the telluric correction to to the .yml file, prob just to make life
    #  easier?

    return flux_calibration, telluric_correction


def get_sky(file, skyflat_red):
    sky_r_in = KOALA_RSS(
        file,
        apply_throughput=True,
        skyflat=skyflat_red,
        do_extinction=False,
        correct_ccd_defects=True,
        correct_high_cosmics=False,
        clip_high=100,
        step_ccd=50,
        sky_method="none",
        is_sky=True,
        win_sky=151,
        plot=DO_PLOTTING,
    )
    sky = sky_r_in.plot_combined_spectrum(
        list_spectra=[
            870,
            871,
            872,
            873,
            875,
            876,
            877,
            878,
            879,
            880,
            881,
            882,
            883,
            884,
            885,
            886,
            887,
            888,
            889,
            900,
        ],
        median=True,
    )
    return sky, sky_r_in
