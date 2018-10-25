"""
This module provides a series of tools that will be used for the study of the Scanning Strategy of 
the LSPE/STRIP experiment. 
"""
import healpy as hp
import numpy as np
from astropy import units as u
import time as timing
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time


### INSTRUMENT CHARACTERISTICS ###
sampling_rate = 50 #Hz
###

### LOCATION INFORMATION ###
LAT = np.array([28, 16, 24]) #deg
LONG = np.array([-16, 38, 32]) #deg
Height = 2400 #m
###

### TIME INFORMATION ###
LCT_start = (0, 0, 0) #h, m, s
LCD_start = (1, 1, 2015) #d, m, y
UTC, DST = (0 , 0) #h
obs_time = (1, 0, 0, 0, 0) #y, d, h, m, s
###


def period2sec(years=0, days=0, hours=0, min=0, sec=0, sidereal=False):
    """
    It converts a period of time in seconds. If sidereal=True it returns the mean 
    sidereal time (sec); otherwise it returns the mean solar time (sec). 
   
    Assumptions (from: "An introduction to modern astrophysics", (Carroll, Ostlie), 2nd Ed.):

    1 min = 60 sec
    1 hour = 60 min 
    1 (mean) solar day = 24 hours
    1 (mean) sidereal day = 23.93447 (23h 56m 4.1s) (with precession of the equinoxes)
    1 solar year = 365 solar days  *Different from (Carroll, Ostlie)
    1 (mean) sidereal year = 365.24219 solar days (without precession of the equinoxes)
    1 (mean) sidereal year = 365.25631 solar days (with precession of the equinoxes)

    N.B. (Carroll, Ostlie) assumes 1 (mean) solar year = 365.24219 solar days.

    This function takes into account the precession of the equinoxes.
    """
    if sidereal:
        return np.around(
            ((((years * 365.25631 * 24) + (days * 23.93447)) + hours) * 60 + min) * 60 + sec)
    return np.around((((years * 365 + days) * 24 + hours) * 60 + min) * 60 + sec)


def sex2dec(angles, radians=False):    
    """
    It converts angles (or time) from sexagesimal to decimal and then radians (if radians = True).
    
    Example
    -------
    >>> 16°30'35"
    >>> sex2dec(np.array([16, 30, 35])) = 16.509722° = 0.28815 rad
    """
    
    if len(angles.shape) == 1:
        sign = -1 if np.sum(angles<0) == 1 else 1
        deg, min, sec = angles[0], angles[1], angles[2]
    else:
        sign = np.ones(len(angles))
        sign[np.sum(angles<0, axis=1) == 1] = -1 
        deg = angles[..., 0]
        min = angles[..., 1]
        sec = angles[..., 2]
    if radians:
        return np.radians((((sec / 60) + min) / 60) + np.abs(deg)) * sign
    return ((((sec / 60) + min) / 60) + np.abs(deg)) * sign


def dec2sex(t):
    """
    It converts times (or angles) from decimal to sexagesimal.

    Example
    -------
    >>> 17.2563888888888 hours
    >>> decimal2sexagesimal(17.2563888888888) = np.array([17, 15, 23]) [hours, min, sec] if times
    >>> decimal2sexagesimal(17.2563888888888) = np.array([17, 15, 23]) [deg, min, sec] if angles
    """
    sign = t / np.abs(t)
    h = np.floor(np.abs(t))
    m = np.floor((np.abs(t) - h) * 60)
    s = (((np.abs(t) - h) * 60) - m) * 60
    if np.array([t]).shape[-1] == 1:
        return np.column_stack((np.floor(h) * sign, m, s))[0]
    return np.column_stack((np.floor(h) * sign, m, s))


def degrees2hours(angles, decimal=False):
    """
    It converts angles from degrees to hours.

    Example
    -------
    >>> 23° 27' 43.56"
    >>> degrees2hours(np.array([23, 27, 43.56]), decimal=True) = 1.56414 hours
    >>> degrees2hours(np.array([23, 27, 43.56])) = np.array([1, 33, 50.904]) [hours, min, sec]
    """
    decimal_degrees = sex2dec(angles)
    h = 24 * decimal_degrees / 360
    if decimal:
        return h
    return dec2sex(h)


def hours2degrees(angles, decimal=False):
    """
    It converts angles from hours to degrees.

    Example
    -------
    >>> 1h 33min 50.904sec
    >>> hours2degrees(np.array([1, 33, 50.904]), decimal=True) = 1.56414 hours 
    >>> hours2degrees(np.array([1, 33, 50.904])) = np.array([1, 33, 50.904]) [hours, min, sec] 
    """
    decimal_hours = sex2dec(angles)
    d = 360 * decimal_hours / 24
    if decimal:
        return d
    return dec2sex(d)


def GreenwichCalendarDate2JulianDate(GD, GM, GY):
    """
    It converts the Greenwich Calendar Date (GCD) to the Julian Date. 

    Parameters
    ----------
    GD : float
         It is the day in the GCD system.
    GM : integer
         It is the month in the GCD system.
    GY : integer 
         It is the year in the GCD system.

    Returns
    -------
    out : float
    """
    if GM < 3 :
        y, m = (GY - 1, GM + 12)
    else:
        y, m = (GY, GM)
    B = 2 - np.int(y / 100) + np.int(np.int(y / 100) / 4) 
    if GY < 1582:
        if GM < 10:
            if GY < 15:
                B = 0
    if y < 0 :
        C = np.int((365.25 * y) - 0.75)
    else:
        C = np.int(365.25 * y)
    D = np.int(30.6001 * (m + 1))
    return B + C + D + GD + 1720994.5


def LocalCivilTime2GreenwichCalendarDay(LCT, LCD, UTC=0, DST=0):
    """
    It converts the Local Civil Time (LCT) to the Greenwich Calendar Date.

    Parameters
    ----------
    LCT : tuple of length 3
          It is the Local Civil Time in the format (hours, minutes, seconds). 
    LCD : tuple of length 3
          It is the Local Civil Day in the format (day, month, year). 
    UTC : integer in the range (-12, 12)
          It is the number which describes the Earth time zone with respect to the Greenwich time 
          zone (UTC=0).
    DST : integer in the range (0, 1)
          It is the number which says if the Daylight Saving Time is on or off.

    Returns
    -------
    out : float, integer, integer [day, month, year]
    """
    (h, m, s), (D, M, Y)  = LCT, LCD
    UniversalTime = (h - DST - UTC) + m / 60 + s / 3600 
    return UniversalTime / 24 + D, M, Y


def LocalCivilTime2JulianDay(LCT, LCD, UTC=0, DST=0):
    """
    It converts the Local Civil Time (LCT) to the Julian Day.

    Parameters
    ----------
    LCT : tuple of length 3
          It is the Local Civil Time in the format (hours, minutes, seconds). 
    LCD : tuple of length 3
          It is the Local Civil Day in the format (day, month, year). 
    UTC : integer in the range (-12, 12)
          It is the number which describes the Earth time zone with respect to the Greenwich time 
          zone (UTC=0).
    DST : integer in the range (0, 1)
          It is the number which says if the Daylight Saving Time is on or off.

    Returns
    -------
    out : float
    """
    GreenwichCalendarDay, M, Y = LocalCivilTime2GreenwichCalendarDay(LCT, LCD, UTC=UTC, DST=DST)
    return GreenwichCalendarDate2JulianDate(GreenwichCalendarDay, M, Y)


def LocalCivilTime2LocalSiderealTime(LCT, LCD, LONG, UTC=0, DST=0):
    """
    It converts the Local Civil Time (LCT) to the Local Sidereal Time.

    Parameters
    ----------
    LCT : tuple of length 3
          It is the Local Civil Time in the format (hours, minutes, seconds). 
    LCD : tuple of length 3
          It is the Local Civil Day in the format (day, month, year). 
    LONG   : tuple of length 3
             It is the observation site Longitude in the format (deg, min, sec). 
    UTC : integer in the range (-12, 12)
          It is the number which describes the Earth time zone with respect to 
          the Greenwich time zone (UTC=0).
    DST : integer in the range (0, 1)
          It is the number which says if the Daylight Saving Time is on or off.

    Returns
    -------
    out : float
    """
    jd = LocalCivilTime2JulianDay(LCT, LCD, UTC=0, DST=0)
    d = jd - 2451543.5
    w = 282.9404 + 4.70935e-5 * d
    M = 356.0470 + 0.9856002585 * d
    L = w + np.mod(M, 360)
    GMST0 = np.mod(L, 360)/15 + 12
    rest = jd - np.floor(jd) + 0.5
    UT = (rest - np.floor(rest)) * 24
    LST = np.mod(GMST0, 24) + UT + sex2dec(LONG)/15
    return dec2sex(np.mod(LST, 24))


def get_nside_eff(fwhm_beam):
    """
    It returns the nside corresponding to the given angular resolution of the beam.
    N.B. The results are slightly differents from hp.nside2resol. A cross-check is strongly 
    raccommended.
    """
    fwhm = sex2dec(fwhm_beam, radians=True)
    area_pix_eff = np.pi * fwhm**2 / 4
    N_pix_eff = 4 * np.pi / area_pix_eff
    nside_eff = np.sqrt(N_pix_eff / 12)
    return 2**(np.int(np.log2(nside_eff)) + 1) 


def get_full_fp(fp_theta_path, fp_phi_path):
    """
    It reads 2 files.txt which are the focal plane coordinates in spherical coordinates of the first
    4 modules of STRIP (I, Y, O, R). It returns the full focal plane coordinates in cartesian 
    coordinates, since the last 3 modules (V, B, G) are symmetrical with respect the central one. 
    IT returns the complete focal plane coordinates and the total number of horns in the focal plane
    with the following conventions: 

    x_fp[0:7]   = module named I
    x_fp[7:14]  = module named Y
    x_fp[14:21] = module named O
    x_fp[21:28] = module named R
    x_fp[28:35] = module named v
    x_fp[35:42] = module named B
    x_fp[42:49] = module named G

    Returns
    -------
    out : numpy array of shape (49, 3), float
    """
    theta = np.loadtxt(fp_theta_path)
    phi = np.loadtxt(fp_phi_path)
    full_theta = np.append(theta, [theta[-7:], theta[-14:-7], theta[-21:-14]])
    full_theta[29:34] = full_theta[29:34][::-1]
    full_theta[36:41] = full_theta[36:41][::-1]
    full_theta[43:48] = full_theta[43:48][::-1]
    full_phi = np.append(phi, [-phi[-7:], -phi[-14:-7], -phi[-21:-14]])
    full_phi[29:34] = full_phi[29:34][::-1]
    full_phi[36:41] = full_phi[36:41][::-1]
    full_phi[43:48] = full_phi[43:48][::-1]
    x_fp = hp.ang2vec(full_theta, full_phi)
    return x_fp, len(x_fp)


def get_full_fp_polarization_angles(fp_psi_path):
    """
    It reads 1 files.txt which cointains the focal plane polarization angles for the first 4 modules
    of STRIP (I, Y, O, R). It returns the full focal plane polarization angles, since the last 3 
    modules (V, B, G) are symmetrical with respect the central one. It uses the same conventions of 
    get_full_fp function. And, furthermore, it uses the following conventions to construct the 
    polarization versor:

    x = np.cos(psi)
    y = np.sin(psi)
    z = 0

    Returns
    -------
    out : numpy array of shape (49,), numpy array of shape (49, 3)
    """
    psi = np.loadtxt(fp_psi_path)
    full_psi = np.append(psi, [-psi[-7:], -psi[-14:-7], -psi[-21:-14]])
    full_psi[29:34] = full_psi[29:34][::-1]
    full_psi[36:41] = full_psi[36:41][::-1]
    full_psi[43:48] = full_psi[43:48][::-1]
    polarization_versors = np.column_stack((np.cos(full_psi), np.sin(full_psi),
                                           np.full_like(full_psi, 0)))
    return full_psi, polarization_versors


def get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=0, DST=0, day=None):
    """
    It returns the total observation time and the time samples expressed in seconds, and the Julian 
    Dates starting from a given Local Civil Time and Day.

    Parameters
    ----------
    LCT_start     : tuple of length 3
                    It is the Local Civil Time in the format (hours, minutes, seconds). 
    LCD_start     : tuple of length 3
                    It is the Local Civil Day in the format (day, month, year). 
    sampling_rate : float, Hertz
                    It is the instrument sampling rate.
    obs_time      : tuple of lenght 5
                    It is the observation time in the format (years, days, hours, minutes, seconds)
    UTC           : integer in the range (-12, 12)
                    It is the number which describes the Earth time zone with respect to the 
                    Greenwich time zone (UTC=0).
    DST           : integer in the range (0, 1)
                    It is the number which says if the Daylight Saving Time is on or off.
    day           : float, default = None
                    In the case of observation time greater than 1 day it is the observation day at 
                    which the time is computed. Otherwise, it must be setted to None.
                    If None it will return the time computed over all the given observation time. 
   
    Returns
    -------
    out : float, numpy array of shape (obs_time * sampling_rate,) or (86400 * sampling_rate,) if day
          is not None, numpy array of shape (obs_time * sampling_rate,) or (86400 * sampling_rate,) 
          if day is not None 
    """
    obs_t = period2sec(years=obs_time[0], days=obs_time[1], hours=obs_time[2], min=obs_time[3],
                       sec=obs_time[4], sidereal=False)
    JD_start = LocalCivilTime2JulianDay(LCT_start, LCD_start, UTC=UTC, DST=DST)
    LCT_step = (0, 0, 1 / sampling_rate)
    JD_step = (
        LocalCivilTime2JulianDay(np.array(LCT_start) + np.array(LCT_step), LCD_start, UTC=UTC,
        DST=DST) - LocalCivilTime2JulianDay(LCT_start, LCD_start, UTC=UTC, DST=DST))                
    if day is None:
        time = np.linspace(0, obs_t, obs_t * sampling_rate + 1) #sec
        JD = np.arange(JD_start, JD_start + JD_step * len(time), JD_step)
    else:
        if obs_t < 86400:
            raise ValueError("If the obs_time is lower than 1 day the 'day' parameter must be None")
        t = period2sec(days=1)
        time = np.linspace((day - 1) * t, day * t, t * sampling_rate + 1)
        JD = np.arange(JD_start + JD_step * (len(time) * (day - 1) - 1),
                       JD_start + JD_step * len(time) * day, JD_step)[:-1]      
    return obs_t, time[:-1], JD[:-1]


def get_location(LAT, LONG, Height):
    """
    It returns an astropy object which gets information about the location of observation site.

    Parameters
    ----------
    LAT    : tuple of length 3
             It is the observation site Latitude in the format (deg, min, sec). 
    LONG   : tuple of length 3
             It is the observation site Longitude in the format (deg, min, sec). 
    Height : integer, meters
             It is the observation site Height expressed in meters.       
    """
    Lat = sex2dec(LAT) #deg
    Long = sex2dec(LONG) #deg
    return EarthLocation(lon=Long, lat=Lat, height=Height)


def spin_generator(time, rpm):
    """
    It returns the phi angle in the case of spinning of the telescope.

    Parameters
    ----------
    time          : numpy array of shape (obs_time * sampling rate), sec
                    It is the time sample expressed in seconds.
    rpm           : integer
                    It is the number of revolutions per minute of the telescope.
    
    Returns
    -------
    out : numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day is not 
          None
    """
    spin_freq = rpm / 60  #Hz
    spin_T = 1 / spin_freq #sec
    return 2 * np.pi * (time / spin_T - (time / spin_T).astype(int))


def euler_rotation_matrix(phi, theta, psi=0):
    """
    It returns the rotation matrix correspondig to the rotation of the angles phi, theta, phi about 
    respectively the axes (z, x, z). So that, the following rotation is performed:
    
    R_z(phi) * R_x(theta) * R_z(psi)
    
    Parameters
    ----------
    phi    : float or numpy array, radians 
             Rotation about the x axis.
    theta  : float or numpy array, radians
             Rotation about the z axis.
    psi    : float or numpy array, radians, default = 0
             Rotation about the x axis.

    Returns
    -------
    out : numpy array of shape (obs_time * sampling_rate) 
    """
    c1, c2, c3 = (np.cos(phi), np.cos(theta), np.cos(psi))
    s1, s2, s3 = (np.sin(phi), np.sin(theta), np.sin(psi))
    matrix =  np.array([[c1*c3 - c2*s1*s3, -c1*s3 - c2*c3*s1, s1*s2],
                        [c3*s1 + c1*c2*s3, c1*c2*c3 - s1*s3, -c1*s2],
                        [s2*s3, c3*s2, c2]])
    return np.rollaxis(matrix, 2, 0)


def get_engine_rotations(time, rpm, zenith_distance, polarization_angle):
    """
    It returns the temporal sequence of the engine rotation angles. It supposes that the zenith 
    distance and the polarization angle are constants while the azimth angle spins around the z 
    axis.

    Parameters
    ----------
    time               : numpy array of shape (obs_time * sampling rate), sec
                         It is the time sample expressed in seconds.
    rpm                : integer
                         It is the number of revolutions per minute of the telescope.
    zenith_distance    : integer, degrees
                         It is the zenith distance in the telescope frame of reference.
    polarization_angle : integer, degrees
                         It is the polarization angle in the telescope frame of reference.

    Returns
    -------
    out : numpy array of shape (obs_time * sampling_rate), 
          numpy array of shape (obs_time * sampling_rate), 
          numpy array of shape (obs_time * sampling_rate) 
    """
    theta = np.full_like(time, np.radians(zenith_distance)) #rad
    phi = spin_generator(time, rpm) #rad
    psi = np.full_like(time, np.radians(polarization_angle)) #rad
    return theta, phi, psi


def get_fp_rotations(phi, theta, psi, x_fp, n_horns, time, n=None, cartesian=False):
    """
    It returns the temporal sequence of pointings for each horn. 

    Parameters
    ----------
    phi       : numpy array of shape (obs_time * sampling_rate), radians 
                The Azimith angle.
    theta     : numpy array of shape (obs_time * sampling_rate), radians 
                The Zenith distance.
    psi       : numpy array of shape (obs_time * sampling_rate), radians 
                The polarization angle.
    x_fp      : numpy array of shape (n_horns, 3)
                They are the positions (versors) of the horns in the frame of reference of the 
                focal plane (the one in which the central horn points towards the local Zenith with 
                versor (0, 0, 1)).
    n_honrs   : integer
                It is the total number of horns.
    time      : numpy array of shape (obs_time * sampling rate), sec
                It is the time sample expressed in seconds.
    n         : integer in the range (0, n_horns - 1), default = None
                It is the horn for which the pointings are computed. If None, will be returned 
                the pointings for all the horns.
    cartesian : boolean, dafault = False
                If True will be returned the pointings in cartesian coordinates. Otherwise, the 
                spherical coordinates of the pointings will be returned.

    Returns
    -------
    out : numpy array of shape (n_horns, obs_time * sampling_rate, 2), 
          numpy array of shape (n_horns, obs_time * sampling_rate, 3) if cartesian, 
          numpy array of shape (obs_time * sampling_rate, 2) if n,
          numpy array of shape (obs_time * sampling_rate, 3) if n and cartesian
    """
    tc, tw = (timing.clock(), timing.time())
    fp_rotations = euler_rotation_matrix(phi, theta, psi)
    if n is None:
        fp_pointings = np.sum(fp_rotations[None, ...] * x_fp[:, None, None, :], axis=-1)
        fp_pointings_spherical = np.column_stack(hp.vec2ang(fp_pointings)).reshape(
            n_horns, len(time), 2)
    else:
        fp_pointings = np.sum(fp_rotations * x_fp[n, None, :], axis=-1)
        fp_pointings_spherical = np.column_stack(hp.vec2ang(fp_pointings))
    clock_time, wall_time = (timing.clock() - tc, timing.time() - tw)
    print ('fp conversion clock time [sec]:', clock_time, '\n', 'fp conversion wall time [sec]:',
           wall_time)
    if cartesian:
        return fp_pointings
    return fp_pointings_spherical 


def get_horizon_coordinates(fp_pointings_spherical):
    """
    It converts from spherical to Horizon coordinates, with the conventions:
    
    Altitute = np.pi / 2 - zenith angle (theta)
    Azimuth  = 2 * np.pi - phi

    Parameters
    ----------
    fp_pointings_spherical : numpy array of shape (..., 2), radians
                             They are the spherical coordinates (theta, phi) that will be converted.

    Returns
    -------
    out : numpy array of shape (..., ), numpy array of shape (..., )  
    """
    Alt = np.pi/2  - fp_pointings_spherical[..., 0] #rad
    Az = 2 * np.pi - fp_pointings_spherical[..., 1] #rad
    return Alt, Az


def get_practical_icrs_coordinates(JD, loc, Alt, Az, hours=False):
    """
    It converts from Horizon to Equatorial (icrs) coordinates, according to the astropy conventions.
    
    Parameters
    ----------
    JD    : numpy array of shape (obs_time * sampling_rate) or (n_horns, obs_time * sampling_rate)
            The Julian Dates corresponding to the temporal sequence of pointings.
    loc   : astropy EarthLocation object
            An astropy object which gets information about the location of observation site.
    Alt   : numpy array of shape (obs_time * sampling_rate) or (n_horns, obs_time * sampling_rate)
            The temporal sequence of Altitude.
    Az    : numpy array of shape (obs_time * sampling_rate) or (n_horns, obs_time * sampling_rate)
            The temporal sequence of Azimuth.
    hours : boolean, default=False
            If True it returns the Right Ascension in hours; otherwise decimal.

    Returns
    -------
    out : numpy array of shape (..., ), numpy array of shape (..., )  
    """
    tc, tw = (timing.clock(), timing.time())
    LST = Time(JD, format='jd', location=loc).sidereal_time('mean').value
    Lat = loc.lat.rad
    Long = loc.lon.rad
    Dec = np.arcsin(np.sin(Alt) * np.sin(Lat) + np.cos(Alt) * np.cos(Lat) * np.cos(Az))
    HourAngle = np.arccos((np.sin(Alt) - np.sin(Dec) * np.sin(Lat)) / (np.cos(Dec) * np.cos(Lat)))
    index = np.sin(Az) < 0
    h = (360 - np.degrees(HourAngle)) / 15
    h[index] = np.degrees(HourAngle[index]) / 15
    Ra = LST - h
    Ra[Ra < 0] += 24
    clock_time, wall_time = (timing.clock() - tc, timing.time() - tw)
    print ('practical icrs conversion clock time [sec]:', clock_time, '\n',
           'practical icrs conversion wall time [sec]:', wall_time)
    if hours:
        return Dec, Ra
    return Dec, np.radians(360 * Ra / 24)


def get_icrs_coordinates(JD, loc, Alt, Az):
    """
    It converts from Horizon to Equatorial (icrs) coordinates, according to the astropy conventions.
    
    Parameters
    ----------
    JD  : numpy array of shape (obs_time * sampling_rate) or (n_horns, obs_time * sampling_rate)
          The Julian Dates corresponding to the temporal sequence of pointings.
    loc : astropy EarthLocation object
          An astropy object which gets information about the location of observation site.
    Alt : numpy array of shape (obs_time * sampling_rate) or (n_horns, obs_time * sampling_rate)
          The temporal sequence of Altitude.
    Az  : numpy array of shape (obs_time * sampling_rate) or (n_horns, obs_time * sampling_rate)
          The temporal sequence of Azimuth.
    
    Returns
    -------
    out : numpy array of shape (..., ), numpy array of shape (..., )  
    """
    tc, tw = (timing.clock(), timing.time())
    times = Time(JD, format='jd', location=loc)
    pointings = SkyCoord(alt=Alt*u.rad, az=Az*u.rad, obstime=times, frame='altaz', location=loc)
    Ra = pointings.transform_to('icrs').ra.rad #rad
    Dec = pointings.transform_to('icrs').dec.rad #rad
    clock_time, wall_time = (timing.clock() - tc, timing.time() - tw)
    print ('icrs conversion clock time [sec]:', clock_time, '\n',
           'icrs conversion wall time [sec]:', wall_time)
    return Dec, Ra


def get_polarization_angles(phi, theta, psi, x_fp_pol_versors, n_horns, time, n=None):
    """
    It returns the polarization angles projected in the sky.

    Parameters
    ----------
    phi              : float or numpy array, radians 
                       Rotation about the x axis.
    theta            : float or numpy array, radians
                       Rotation about the z axis.
    psi              : float or numpy array, radians, default = 0
                       Rotation about the x axis.
    x_fp_pol_versors : numpy array of shape (n_horns, 3)
                       They are the polarization versors for each horn in the frame of reference of 
                       the focal plane (the one in which the central horn points towards the local 
                       Zenith with versor (0, 0, 1)).
    n_honrs          : integer
                       It is the total number of horns.
    time             : numpy array of shape (obs_time * sampling rate), sec
                       It is the time sample expressed in seconds.
    n                : integer in the range (0, n_horns - 1), default = None
                       It is the horn for which the pointings are computed. If None, will be 
                       returned the pointings for all the horns.

    Returns
    -------
    out : numpy array of shape (n_horns, obs_time * sampling_rate) or 
          numpy array of shape (obs_time * sampling_rate, ) if n  
    """
    fp_pol_pointings = get_fp_rotations(phi, theta, psi, x_fp_pol_versors, n_horns, time, n=n,
                                        cartesian=True) #rad
    return np.arctan2(fp_pol_pointings[..., 1], fp_pol_pointings[..., 0])


def get_scanning_strategy(obs_time, sampling_rate, zenith_distance, boresight_angle, rpm, n=None,
                          day=None, LCT_start=(0, 0, 0), LCD_start=(1, 1, 2018), UTC=0, DST=0,
                          LAT=np.array([28, 16, 24]), LONG=np.array([-16, 38, 32]), Height=2400,
                          fp_theta_path='./fp_theta.txt', fp_phi_path='./fp_phi.txt',
                          fp_psi_path='./fp_psi.txt'):
    """
    It returns all the parameters of the STRIP Scanning Strategy, in the following order:
    
    - The location of the horns in Cartesian coordinates        : x_fp;
    - The total number of horns                                 : n_honrs;
    - The time sample                                           : time;
    - The Julian Dates sample                                   : JD;
    - The sample of the engine rotation angles                  : theta, phi, psi;
    - The focal plane pointings sample in spherical coordinates : fp_pointings_spherical;
    - The Horizon coordinates                                   : Alt, AZ;
    - The Equatorial (icrs) coordinates                         : Dec, Ra.

    In particular, if n is specified will be returned only the values for the horn n. Otherwise, 
    will be returned the values for all the horns. In the same way, if day is specified will be 
    returned the values for that obsarvation day. Otherwise, will be returned the values for the 
    whole observation period. 

    Parameters
    ----------
    obs_time        : tuple of lenght 5
                      It is the observation time in the format (years, days, hours, minutes, 
                      seconds)
    sampling_rate   : float, Hertz
                      It is the instrument sampling rate.
    zenith_distance : integer, degrees
                      It is the zenith distance in the telescope frame of reference.
    boresight_angle : integer, degrees
                      It is the polarization angle in the telescope frame of reference.
    rpm             : integer
                      It is the number of revolutions per minute of the telescope.
    n               : integer in the range (0, n_horns - 1), default = None
                      It is the horn for which the pointings are computed. If None, will be 
                      returned the pointings for all the horns.
    day             : float, default = None
                      In the case of observation time greater than 1 day it is the observation 
                      day at which the time is computed. Otherwise, it must be setted to None. If
                      None it will return the time computed over all the given observation time. 
    LCT_start       : tuple of length 3
                      It is the Local Civil Time in the format (hours, minutes, seconds). 
    LCD_start       : tuple of length 3
                      It is the Local Civil Day in the format (day, month, year). 
    UTC             : integer in the range (-12, 12)
                      It is the number which describes the Earth time zone with respect to the 
                      Greenwich time zone (UTC=0).
    DST             : integer in the range (0, 1)
                      It is the number which says if the Daylight Saving Time is on or off.
    LAT             : tuple of length 3
                      It is the observation site Latitude in the format (deg, min, sec). 
    LONG            : tuple of length 3
                      It is the observation site Longitude in the format (deg, min, sec). 
    Height          : integer, meters
                      It is the observation site Height expressed in meters.       
    fp_theta_path   : string
                      It is the path to the file.txt where are stored the theta positions of the 
                      focal plane horns.  
    fp_phi_path     : string
                      It is the path to the file.txt where are stored the phi positions of the 
                      focal plane horns.
    fp_psi_path     : string
                      It is the path to the file.txt where are stored the polarization angles of
                      the focal plane horns.

    Returns
    -------
    out : numpy array of shape (n_horns, 3),
          numpy array of shape (n_horns, ), 
          float, 
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,  
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate, 2) or (86400 * sampling_rate, 2) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day,
          numpy array of shape (obs_time * sampling_rate) or (86400 * sampling_rate) if day
    """
    
    ### FOCAL PLANE LOAD ###
    x_fp, n_horns = get_full_fp(fp_theta_path, fp_phi_path)
    x_fp_pol_angles, x_fp_pol_versors = get_full_fp_polarization_angles(fp_psi_path)
    ### GET LOCATION ###
    loc = get_location(LAT=LAT, LONG=LONG, Height=Height)
    ### GET TIME INFORMATION ###
    obs_t, time, JD = get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC, DST=DST,
                                 day=day)
    ### GET ENGINE ROTATIONS ###
    theta, phi, psi = get_engine_rotations(time, rpm, zenith_distance, boresight_angle) #rad
    ### GET FOCAL PLANE POINTINGS ###
    fp_pointings_spherical = get_fp_rotations(phi, theta, psi, x_fp, n_horns, time, n=n) #rad
    ### GET HORIZON COORDINATES ###
    Alt, Az = get_horizon_coordinates(fp_pointings_spherical) #rad
    ### EQUATORIAL (ICRS) COORDINATES CONVERSION ###
    Dec, Ra =  get_practical_icrs_coordinates(JD, loc, Alt, Az) #rad
    ### GET POLARIZATION ANGLES ###
    polarization_angles = get_polarization_angles(phi, theta, psi, x_fp_pol_versors, n_horns, time,
                                                  n=n)
    
    return (x_fp, x_fp_pol_angles, n_horns, time, JD, theta, phi, psi, fp_pointings_spherical, Alt,
            Az, Dec, Ra, polarization_angles)


