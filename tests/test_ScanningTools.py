""""
STRIP Scanning Strategy Tools test module.
"""

import unittest
import healpy as hp
import numpy as np
from ScanningTools import ScanningTools as st
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz 
from ScanningTools.Quaternions import Quaternion as q


angles = np.array([[-10, 45, 59],
                 [30, 35, 15],
                 [-180, 25, 20],
                 [3, 4, 5]])

hours = np.array([[23, 59, 16],
                  [7, 56, 59]])

t = np.array([1.546585, -0.56, 0.3333333333333333333, -1.001])

### INSTRUMENT CHARACTERISTICS ###
pointing_accuracy = np.array([0, 0, 25]) #deg (arcsec)
###

### LOCATION INFORMATION ###
LAT = np.array([28, 16, 24]) #deg
LONG = np.array([-16, 38, 32]) #deg
Height = 2400 #m

loc = st.get_location(LAT, LONG, Height)
###

### TIME INFORMATION ###
LCT_start = (0, 0, 0) #h, m, s
LCD_start = (1, 1, 2015) #g, m, y
UTC, DST = (0 , 0) #h
###

### ENGINE ROTATION INFORMATION ###
zenith_distance = 31 #deg
polarization_angle = 60 #deg
###


class TestScanningTools(unittest.TestCase):

    
    # def test_period2sec(self):
        
    #     one_sidereal_year = st.period2sec(years=1, days=0, hours=0, min=0, sec=0, sidereal=True)
    #     one_solar_year = st.period2sec(years=1, days=0, hours=0, min=0, sec=0)
    #     one_sidereal_day = st.period2sec(years=0, days=1, hours=0, min=0, sec=0, sidereal=True)
    #     one_solar_day = st.period2sec(years=0, days=1, hours=0, min=0, sec=0)
    #     period_0 = st.period2sec(years=1, days=1, hours=0, min=0, sec=0, sidereal=True)
    #     period_1 = st.period2sec(years=5, days=30, hours=0, min=0, sec=0, sidereal=True)
    #     period_2 = st.period2sec(years=2, days=17, hours=0, min=0, sec=0, sidereal=True)
    #     period_3 = st.period2sec(years=10, days=21, hours=15, min=3, sec=25, sidereal=True)
    #     self.assertEqual(one_sidereal_year, 31558145)
    #     self.assertEqual(one_solar_year, 31536000)
    #     self.assertEqual(one_sidereal_day, 86164)
    #     self.assertEqual(one_solar_day, 86400)
    #     self.assertEqual(period_0, 31644309)
    #     self.assertEqual(period_1, 160375649)
    #     self.assertEqual(period_2, 64581080)
    #     self.assertEqual(period_3, 317445103) 


    # def test_sex2dec(self):
        
    #     ang0 = st.sex2dec(angles)
    #     ang1 = st.sex2dec(angles[0], radians=True)
    #     self.assertTrue(np.allclose(ang0, np.array([-10.76638889, 30.587500, -180.422222,
    #                                                 3.06805556]))) 
    #     self.assertEqual(ang1, np.radians(ang0[0])) 


    # def test_dec2sex(self):
        
    #     t0 = st.dec2sex(t)
    #     t00 = st.dec2sex(t[0])
    #     self.assertTrue(np.allclose(t0, np.array([[1, 32, 47.706], [-0, 33, 36], [0, 20, 0],
    #                                               [-1, 0, 3.6]])))
    #     self.assertTrue(np.allclose(t00, np.array([1, 32, 47.706])))

        
    # def test_degrees2hours(self):
        
    #     ang0 = st.degrees2hours(angles)
    #     ang1 = st.degrees2hours(angles[2], decimal=True)
    #     self.assertTrue(np.allclose(ang0, st.dec2sex(st.sex2dec(angles) / 15)))
    #     self.assertTrue(np.allclose(ang1, st.sex2dec(angles)[2] / 15))

        
    # def test_hours2degrees(self):
        
    #     ang0 = st.hours2degrees(hours[1])
    #     ang1 = st.hours2degrees(hours, decimal=True)
    #     self.assertTrue(np.allclose(ang0, st.dec2sex(st.sex2dec(hours[1]) * 15)))
    #     self.assertTrue(np.allclose(ang1, st.sex2dec(hours) * 15))

        
    # def test_LocalCivilTime2JulianDay(self):
        
    #     "Integrated Test: it includes also the LCT2GCD and GCD2JD function conversion"
    #     Jul_1_2013 = st.LocalCivilTime2JulianDay((3, 37, 0), (1, 7, 2013), UTC=4, DST=1)
    #     Jun_19_2009 = st.LocalCivilTime2JulianDay((18, 0, 0), (19, 6, 2009), UTC=0, DST=0)
    #     self.assertTrue(np.allclose(Jul_1_2013, 2456474.442))
    #     self.assertTrue(np.allclose(Jun_19_2009, 2455002.25))
    #     t = Time(['2015-1-1 00:00:10', '2018-1-3 5:15:24.3', '1980-4-22 19:30:2']).jd
    #     T = np.array([st.LocalCivilTime2JulianDay((0, 0, 10), (1, 1, 2015), UTC=0, DST=0),
    #                   st.LocalCivilTime2JulianDay((5, 15, 24.3), (3, 1, 2018), UTC=0, DST=0),
    #                   st.LocalCivilTime2JulianDay((19, 30, 2), (22, 4, 1980), UTC=0, DST=0)])
    #     self.assertTrue(np.allclose(t, T))

        
    # def test_get_nside_eff(self):
        
    #     fwhm_beam0 = np.array([0, 5, 0]) #deg (arcmin)
    #     fwhm_beam1 = np.array([0, 21, 0]) #deg (arcmin)
    #     fwhm_beam2 = np.array([0, 32, 0]) #deg (arcmin)
    #     self.assertEqual(st.get_nside_eff(fwhm_beam0), 1024)
    #     self.assertEqual(st.get_nside_eff(fwhm_beam1), 256)
    #     self.assertEqual(st.get_nside_eff(fwhm_beam2), 128)

        
    # def test_get_full_fp(self):
        
    #     def general_test(x_fp, i, j):
    #         self.assertTrue(np.allclose(x_fp[i, 0], x_fp[j, 0]))
    #         self.assertTrue(np.allclose(x_fp[i, 1], -x_fp[j, 1]))
    #         self.assertTrue(np.allclose(x_fp[i, 2], x_fp[j, 2]))
            
    #     x_fp, n_horns = st.get_full_fp('./ScanningTools/fp_data/fp_theta.txt',
    #                                    './ScanningTools/fp_data/fp_phi.txt')
    #     self.assertTrue(np.allclose(np.sum(x_fp**2, axis=1), 1))
    #     self.assertEqual(n_horns, 49)
    #     general_test(x_fp, 7, 42)
    #     general_test(x_fp, 8, 47)
    #     general_test(x_fp, 9, 46)        
    #     general_test(x_fp, 10, 45)
    #     general_test(x_fp, 11, 44)
    #     general_test(x_fp, 12, 43)
    #     general_test(x_fp, 13, 48)        
    #     general_test(x_fp, 14, 35)
    #     general_test(x_fp, 15, 40)
    #     general_test(x_fp, 16, 39)
    #     general_test(x_fp, 17, 38)
    #     general_test(x_fp, 18, 37)
    #     general_test(x_fp, 19, 36)
    #     general_test(x_fp, 20, 41)
    #     general_test(x_fp, 21, 28)
    #     general_test(x_fp, 22, 33)
    #     general_test(x_fp, 23, 32)
    #     general_test(x_fp, 24, 31)
    #     general_test(x_fp, 25, 30)
    #     general_test(x_fp, 26, 29)
    #     general_test(x_fp, 27, 34)

        
    # def get_full_fp_polarization_angles(self):
        
    #     def general_test(x_fp, i, j):
    #         self.assertTrue(np.allclose(x_fp[i, 0], x_fp[j, 0]))
    #         self.assertTrue(np.allclose(x_fp[i, 1], -x_fp[j, 1]))
    #         self.assertTrue(np.allclose(x_fp[i, 2], x_fp[j, 2]))
            
    #     full_psi, polarization_versor = st.get_full_fp_polarization_angles(
    #         './ScanningTools/fp_data/fp_psi.txt')
    #     self.assertTrue(np.allclose(np.sum(polarization_versor**2, axis=1), 1))
    #     self.assertEqual(len(full_psi), 49)
    #     self.assertEqual(len(polarization_versor), 49)        
    #     general_test(polarization_versor, 7, 42)
    #     general_test(polarization_versor, 8, 47)
    #     general_test(polarization_versor, 9, 46)        
    #     general_test(polarization_versor, 10, 45)
    #     general_test(polarization_versor, 11, 44)
    #     general_test(polarization_versor, 12, 43)
    #     general_test(polarization_versor, 13, 48)        
    #     general_test(polarization_versor, 14, 35)
    #     general_test(polarization_versor, 15, 40)
    #     general_test(polarization_versor, 16, 39)
    #     general_test(polarization_versor, 17, 38)
    #     general_test(polarization_versor, 18, 37)
    #     general_test(polarization_versor, 19, 36)
    #     general_test(polarization_versor, 20, 41)
    #     general_test(polarization_versor, 21, 28)
    #     general_test(polarization_versor, 22, 33)
    #     general_test(polarization_versor, 23, 32)
    #     general_test(polarization_versor, 24, 31)
    #     general_test(polarization_versor, 25, 30)
    #     general_test(polarization_versor, 26, 29)
    #     general_test(polarization_versor, 27, 34)
        
        
    # def test_get_timeJD(self):
        
    #     def general_tests(time, sampling_rate, JD, JD_step, t0, t1):
    #         self.assertTrue(np.allclose(time[1:] - time[0:-1], 1 / sampling_rate))
    #         self.assertEqual(len(JD), len(time))
    #         self.assertEqual(np.sum(np.diff(JD_step)), 0)
    #         self.assertTrue(np.allclose((t1-t0).sec, 1 / sampling_rate, rtol=1e-3))
            
    #     def tests_1h(obs_t, time, sampling_rate, JD, JD_step, t0, t1):
    #         self.assertEqual(obs_t, 3600)
    #         self.assertEqual(len(time), obs_t * sampling_rate)       
    #         general_tests(time, sampling_rate, JD, JD_step, t0, t1)
            
    #     def tests_1d(LCT_start, LCD_start, obs_t, time, sampling_rate, JD, JD_step, t0, t1, UTC=UTC,
    #                  DST=DST):
    #         general_tests(time, sampling_rate, JD, JD_step, t0, t1)
    #         obs_t0, time0, JD0 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time,
    #                                            UTC=UTC, DST=DST, day=None)
    #         self.assertEqual(obs_t, 86400)
    #         self.assertEqual(obs_t, obs_t0)
    #         self.assertEqual(len(time), obs_t * sampling_rate)       
    #         self.assertTrue(len(time), len(time0))
    #         self.assertTrue(len(JD), len(JD0))

    #     def tests_1y(LCT_start, LCD_start, obs_t, time, sampling_rate, JD, JD_step, t0, t1, UTC=UTC,
    #                  DST=DST, day=None):
    #         general_tests(time, sampling_rate, JD, JD_step, t0, t1)
    #         self.assertEqual(obs_t, 86400 * 365)
    #         if day:
    #             self.assertEqual(len(time), 86400 * sampling_rate)
    #             if day > 1:
    #                 self.assertTrue(time[0] != 0)
    #             else:
    #                 self.assertTrue(time[0] == 0)
    #         else:
    #             self.assertEqual(len(time), obs_t * sampling_rate)       

    #     sampling_rate = 50 #Hz
    #     obs_time = (0, 0, 1, 0, 0) #y, d, h, m, s
    #     day = None
    #     obs_t, time, JD = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                     DST=DST, day=day)
    #     JD_step = JD[1:] - JD[:-1]
    #     t0 = Time(JD[0], format='jd', location=loc)
    #     t1 = Time(JD[0] + JD_step[0], format='jd', location=loc)
    #     tests_1h(obs_t, time, sampling_rate, JD, JD_step, t0, t1)
        
    #     sampling_rate = 5 #Hz
    #     obs_t, time, JD = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                     DST=DST, day=day)
    #     JD_step = JD[1:] - JD[:-1]
    #     t0 = Time(JD[0], format='jd', location=loc)
    #     t1 = Time(JD[0] + JD_step[0], format='jd', location=loc)
    #     tests_1h(obs_t, time, sampling_rate, JD, JD_step, t0, t1)
        
    #     sampling_rate = 3 #Hz
    #     obs_time = (0, 1, 0, 0, 0) #y, d, h, m, s
    #     day = 1
    #     obs_t, time, JD = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                     DST=DST, day=day)
    #     JD_step = JD[1:] - JD[:-1]
    #     t0 = Time(JD[0], format='jd', location=loc)
    #     t1 = Time(JD[0] + JD_step[0], format='jd', location=loc)
    #     tests_1d(LCT_start, LCD_start, obs_t, time, sampling_rate, JD, JD_step, t0, t1, UTC=UTC,
    #              DST=DST)

    #     sampling_rate = 1 #Hz
    #     obs_time = (1, 0, 0, 0, 0) #y, d, h, m, s
    #     day0, day1, day2, day3, day4 = (1, 5, 364, None, None) 
    #     obs_t0, time0, JD0 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day0)
    #     JD_step0 = JD0[1:] - JD0[:-1]
    #     t00 = Time(JD0[0], format='jd', location=loc)
    #     t10 = Time(JD0[0] + JD_step0[0], format='jd', location=loc)
    #     tests_1y(LCT_start, LCD_start, obs_t0, time0, sampling_rate, JD0, JD_step0, t00, t10,
    #              UTC=UTC, DST=DST, day=day0)
    #     obs_t1, time1, JD1 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                     DST=DST, day=day1)
    #     JD_step1 = JD1[1:] - JD1[:-1]
    #     t01 = Time(JD1[0], format='jd', location=loc)
    #     t11 = Time(JD1[0] + JD_step1[0], format='jd', location=loc)
    #     tests_1y(LCT_start, LCD_start, obs_t1, time1, sampling_rate, JD1, JD_step1, t01, t11,
    #              UTC=UTC, DST=DST, day=day1)
    #     obs_t2, time2, JD2 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day2)
    #     JD_step2 = JD2[1:] - JD2[:-1]
    #     t02 = Time(JD2[0], format='jd', location=loc)
    #     t12 = Time(JD2[0] + JD_step2[0], format='jd', location=loc)
    #     tests_1y(LCT_start, LCD_start, obs_t2, time2, sampling_rate, JD2, JD_step2, t02, t12,
    #              UTC=UTC, DST=DST, day=day2)
    #     obs_t3, time3, JD3 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day3)
    #     JD_step3 = JD3[1:] - JD3[:-1]
    #     t03 = Time(JD3[0], format='jd', location=loc)
    #     t13 = Time(JD3[0] + JD_step3[0], format='jd', location=loc)
    #     tests_1y(LCT_start, LCD_start, obs_t3, time3, sampling_rate, JD3, JD_step3, t03, t13,
    #              UTC=UTC, DST=DST, day=day3)
    #     LCT_start4 = (12, 0, 0)
    #     obs_t4, time4, JD4 = st.get_timeJD(LCT_start4, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day4)
    #     JD_step4 = JD4[1:] - JD4[:-1]
    #     t04 = Time(JD4[0], format='jd', location=loc)
    #     t14 = Time(JD4[0] + JD_step4[0], format='jd', location=loc)
    #     tests_1y(LCT_start4, LCD_start, obs_t4, time4, sampling_rate, JD4, JD_step4, t04, t14,
    #              UTC=UTC, DST=DST, day=day4)
        
        
    # def test_spin_generator(self):
        
    #     def general_spin_tests(phi, obs_time, time, sampling_rate, rpm, day=None):
    #         if day:
    #             self.assertEqual(len(phi), 86400 * sampling_rate)
    #         else:
    #             self.assertEqual(len(phi), obs_time * sampling_rate)
    #         self.assertEqual(
    #             np.sum(np.r_[True, phi[1:] > phi[:-1]] & np.r_[phi[:-1] > phi[1:], True]),
    #             rpm * len(phi) / sampling_rate / 60)
    #         self.assertEqual(phi.min(), 0)
    #         self.assertTrue(phi.max() < 2  * np.pi)
            
    #     obs_time1, obs_time2, obs_time3 = ((0, 30, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0))
    #     sampling_rate1, sampling_rate2, sampling_rate3 = (1, 3, 50)
    #     rpm1, rpm2, rpm3 = (13, 1, 5)
    #     day1, day2, day3 = (2, None, None)  
    #     obs_t1, time1, JD1 = st.get_timeJD(LCT_start, LCD_start, sampling_rate1, obs_time1, UTC=UTC,
    #                                        DST=DST, day=day1)
    #     phi1 = st.spin_generator(time1, rpm1)
    #     general_spin_tests(phi1, obs_t1, time1, sampling_rate1, rpm1, day=day1)
    #     obs_t2, time2, JD2 = st.get_timeJD(LCT_start, LCD_start, sampling_rate2, obs_time2, UTC=UTC,
    #                                        DST=DST, day=day2)
    #     phi2 = st.spin_generator(time2, rpm2)
    #     general_spin_tests(phi2, obs_t2, time2, sampling_rate2, rpm2, day=day2)
    #     obs_t3, time3, JD3 = st.get_timeJD(LCT_start, LCD_start, sampling_rate3, obs_time3, UTC=UTC,
    #                                        DST=DST, day=day3)
    #     phi3 = st.spin_generator(time3, rpm3)
    #     general_spin_tests(phi3, obs_t3, time3, sampling_rate3, rpm3, day=day3)


    # def test_euler_rotation_matrix(self):
        
    #     phi1, theta1, psi1 = np.radians(([10, 10, 10], [30, 30, 30], [0, 0, 0]))
    #     m1 = st.euler_rotation_matrix(phi1, theta1, psi1)
    #     M1 = np.array([[0.98480775301220802, -0.1503837331804353, 0.086824088833465152],
    #                    [0.17364817766693033, 0.85286853195244328, -0.49240387650610395],
    #                    [0, 0.49999999999999994, 0.86602540378443871]])
    #     phi2, theta2, psi2 = np.radians(([10, 10, 10], [30, 30, 30], [45, 45, 45]))
    #     m2 = st.euler_rotation_matrix(phi2, theta2, psi2)
    #     M2 = np.array([[0.59002688280798476, -0.80270159783205308, 0.086824088833465152],
    #                    [0.72585692637316113, 0.4802813184352156, -0.49240387650610395],
    #                    [0.35355339059327368, 0.35355339059327373, 0.86602540378443871]])
    #     phi3, theta3, psi3 = np.radians(([10, 10, 10], [30, 30, 30], [45, 0, 45]))
    #     m3 = st.euler_rotation_matrix(phi3, theta3, psi3)
    #     M3 = np.array([[[0.59002688280798476, -0.80270159783205308, 0.086824088833465152],
    #                    [0.72585692637316113, 0.4802813184352156, -0.49240387650610395],
    #                     [0.35355339059327368, 0.35355339059327373, 0.86602540378443871]],
    #                   [[0.98480775301220802, -0.1503837331804353, 0.086824088833465152],
    #                    [0.17364817766693033, 0.85286853195244328, -0.49240387650610395],
    #                    [0, 0.49999999999999994, 0.86602540378443871]],
    #                    [[0.59002688280798476, -0.80270159783205308, 0.086824088833465152],
    #                    [0.72585692637316113, 0.4802813184352156, -0.49240387650610395],
    #                     [0.35355339059327368, 0.35355339059327373, 0.86602540378443871]]])
    #     self.assertTrue(np.allclose(m1, np.repeat(M1[None, ...], 3, axis=0)))
    #     self.assertTrue(np.allclose(m2, np.repeat(M2[None, ...], 3, axis=0)))
    #     self.assertTrue(np.allclose(m3, M3))

        
    # def test_engine_rotations(self):
        
    #     obs_time = (0, 0, 0, 30, 0)
    #     obs_t, time, JD = st.get_timeJD(LCT_start, LCD_start, 50, obs_time, UTC=UTC, DST=DST,
    #                                     day=None)
    #     rpm = 7
    #     theta, phi, psi = st.get_engine_rotations(time, rpm, zenith_distance, polarization_angle)
    #     self.assertEqual(len(theta), len(time))
    #     self.assertEqual(len(phi), len(time))
    #     self.assertEqual(len(psi), len(time))
    #     self.assertTrue(np.allclose(theta[1:] - theta[:-1], 0))
    #     self.assertTrue(np.allclose(psi[1:] - psi[:-1], 0))

        
    # def test_fp_rotations(self):

    #     def general_tests(fp_pointings, fp_pointings_c):
    #         self.assertTrue(np.allclose(np.diff(fp_pointings_c[..., 2], axis=-1), 0))
    #         self.assertTrue(np.degrees(fp_pointings[..., 0]).max() <= 365)
    #         self.assertTrue(np.degrees(fp_pointings[..., 1]).max() <= 365)
    #         self.assertTrue(np.allclose(np.sum(fp_pointings_c**2, axis=-1), 1))

    #     x_fp, n_horns = st.get_full_fp('./ScanningTools/fp_data/fp_theta.txt',
    #                                    './ScanningTools/fp_data/fp_phi.txt')
    #     obs_time1, obs_time2 = ((0, 0, 0, 30, 0), (0, 90, 0, 0, 0))
    #     rpm1, rpm2 = (7, 2)
    #     day1, day2 = (None, 10)
    #     sampling_rate1, sampling_rate2 = (50, 1)
    #     obs_t1, time1, JD1 = st.get_timeJD(LCT_start, LCD_start, sampling_rate1, obs_time1, UTC=UTC,
    #                                        DST=DST, day=day1)
    #     theta1, phi1, psi1 = st.get_engine_rotations(time1, rpm1, zenith_distance,
    #                                                  polarization_angle)
    #     obs_t2, time2, JD2 = st.get_timeJD(LCT_start, LCD_start, sampling_rate2, obs_time2, UTC=UTC,
    #                                        DST=DST, day=day2)
    #     theta2, phi2, psi2 = st.get_engine_rotations(time2, rpm2, zenith_distance,
    #                                                  polarization_angle)
    #     n1 = 30
    #     n2 = None
    #     fp_rot1 = st.euler_rotation_matrix(phi1, theta1, psi1)
    #     fp_rot2 = st.euler_rotation_matrix(phi2, theta2, psi2)
    #     fp_pointings1 = st.get_fp_rotations(phi1, theta1, psi1, x_fp, n_horns, time1, n=n1,
    #                                         cartesian=False)
    #     fp_pointings1_c = st.get_fp_rotations(phi1, theta1, psi1, x_fp, n_horns, time1, n=n1,
    #                                           cartesian=True)
    #     fp_pointings2 = st.get_fp_rotations(phi2, theta2, psi2, x_fp, n_horns, time2, n=n2,
    #                                         cartesian=False)
    #     fp_pointings2_c = st.get_fp_rotations(phi2, theta2, psi2, x_fp, n_horns, time2, n=n2,
    #                                         cartesian=True)
    #     i = np.random.randint(0, len(time1))
    #     self.assertTrue(np.allclose(np.dot(fp_rot1[i], x_fp[n1]), fp_pointings1_c[i]))
    #     rot1 = q.get_quaternion_from_euler(phi1[i], theta1[i], psi1[i])
    #     self.assertTrue(np.allclose(rot1.rotate_vector_by_quaternion(x_fp[n1]).get_versor(),
    #                                 fp_pointings1_c[i, :]))
    #     general_tests(fp_pointings1, fp_pointings1_c)
    #     j = np.random.randint(0, len(time2))
    #     p = np.random.randint(0, n_horns)
    #     self.assertTrue(np.allclose(np.dot(fp_rot2[j], x_fp[p]), fp_pointings2_c[p][j]))
    #     rot2 = q.get_quaternion_from_euler(phi2[j], theta2[j], psi2[j])
    #     self.assertTrue(np.allclose(rot2.rotate_vector_by_quaternion(x_fp[p]).get_versor(),
    #                                 fp_pointings2_c[p, j, :]))
    #     general_tests(fp_pointings2, fp_pointings2_c)
        

    # def test_get_horizon_coordinates(self):

    #     def general_tests(Alt, Az):
    #         self.assertTrue(np.degrees(Alt.max()) <= 90)
    #         self.assertTrue(np.degrees(Alt.min()) >= 0)
    #         self.assertTrue(np.degrees(Az.max()) <= 360)
    #         self.assertTrue(np.degrees(Az.min()) >= 0)

    #     x_fp, n_horns = st.get_full_fp('./ScanningTools/fp_data/fp_theta.txt',
    #                                    './ScanningTools/fp_data/fp_phi.txt')
    #     obs_time = (0, 2, 0, 0, 0)
    #     sampling_rate = 1
    #     rpm = 4
    #     day = 1
    #     n1, n2 = (0, 15)
    #     obs_t, time, JD = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                     DST=DST, day=day)
    #     theta, phi, psi = st.get_engine_rotations(time, rpm, zenith_distance, polarization_angle)
    #     fpp = st.get_fp_rotations(phi, theta, psi, x_fp, n_horns, time, n=None, cartesian=False)
    #     fpp1 = st.get_fp_rotations(phi, theta, psi, x_fp, n_horns, time, n=n1, cartesian=False)
    #     fpp2 = st.get_fp_rotations(phi, theta, psi, x_fp, n_horns, time, n=n2, cartesian=False)
    #     Alt, Az = st.get_horizon_coordinates(fpp)
    #     Alt1, Az1 = st.get_horizon_coordinates(fpp1)
    #     Alt2, Az2 = st.get_horizon_coordinates(fpp2)


    # def test_get_icrs_coordinates(self):
        
    #     x_fp, n_horns = st.get_full_fp('./ScanningTools/fp_data/fp_theta.txt',
    #                                    './ScanningTools/fp_data/fp_phi.txt')
    #     obs_time = (0, 2, 0, 0, 0)
    #     sampling_rate = 1
    #     rpm = 3
    #     day1, day2 = (2, None)
    #     n1, n2 = (48, None)
    #     obs_t1, time1, JD1 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day1)
    #     obs_t2, time2, JD2 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day2)
    #     theta1, phi1, psi1 = st.get_engine_rotations(time1, rpm, zenith_distance,
    #                                                  polarization_angle)
    #     theta2, phi2, psi2 = st.get_engine_rotations(time2, rpm, zenith_distance,
    #                                                  polarization_angle)  
    #     fpp1 = st.get_fp_rotations(phi1, theta1, psi1, x_fp, n_horns, time1, n=n1, cartesian=False)
    #     fpp2 = st.get_fp_rotations(phi2, theta2, psi2, x_fp, n_horns, time2, n=n1, cartesian=False)
    #     fpp3 = st.get_fp_rotations(phi1, theta1, psi1, x_fp, n_horns, time1, n=n2, cartesian=False)
    #     fpp4 = st.get_fp_rotations(phi2, theta2, psi2, x_fp, n_horns, time2, n=n2, cartesian=False)
    #     Alt1, Az1 = st.get_horizon_coordinates(fpp1)
    #     Alt2, Az2 = st.get_horizon_coordinates(fpp2)
    #     Alt3, Az3 = st.get_horizon_coordinates(fpp3)
    #     Alt4, Az4 = st.get_horizon_coordinates(fpp4)
    #     Dec1, Ra1 = st.get_icrs_coordinates(JD1, loc, Alt1, Az1)      #1 day 2; n = 48 
    #     Dec2, Ra2 = st.get_icrs_coordinates(JD2, loc, Alt2, Az2)      #2 day all; n = 48
    #     Dec3, Ra3 = st.get_icrs_coordinates(JD1, loc, Alt3, Az3)      #3 day 2; n = all
    #     Dec4, Ra4 = st.get_icrs_coordinates(JD2, loc, Alt4, Az4)      #4 day all; n = all 
    #     Dec5, Ra5 = st.get_icrs_coordinates(JD1[0], loc, Alt1, Az1)   #5 day 2 [t=0]; n = 48
    #     self.assertTrue(np.allclose(Dec1, Dec3[n1]))
    #     self.assertTrue(np.allclose(Ra1, Ra3[n1]))
    #     self.assertTrue(np.allclose(Dec1, Dec2[len(Dec1):]))
    #     self.assertTrue(np.allclose(Ra1, Ra2[len(Ra1):]))
    #     self.assertTrue(np.allclose(Dec1, Dec4[n1, len(Dec1):]))
    #     self.assertTrue(np.allclose(Ra1, Ra4[n1, len(Dec1):]))
    #     self.assertTrue(np.allclose(Dec1[0], Dec5[0]))
    #     self.assertTrue(np.allclose(Ra1[0], Ra5[0]))
    #     self.assertFalse(np.allclose(Dec1[1:], Dec5[1:]))
    #     self.assertFalse(np.allclose(Ra1[1:], Ra5[1:]))

    #     def pointing_test(name, JD, loc):
    #         object = SkyCoord.from_name(name)
    #         object_AltAz = object.transform_to(AltAz(obstime=Time(JD, format='jd'), location=loc))
    #         object_Alt = object_AltAz.alt.rad
    #         object_Az = object_AltAz.az.rad
    #         object_Dec, object_Ra = st.get_icrs_coordinates(JD1, loc, object_Alt, object_Az)
    #         self.assertTrue(np.allclose(object.dec.rad, object_Dec))
    #         self.assertTrue(np.allclose(object.ra.rad, object_Ra))        

    #     name1, name2, name3 = ('M33', 'crab', 'NCG67')
    #     pointing_test(name1, JD1, loc)
    #     pointing_test(name2, JD1, loc)
    #     pointing_test(name3, JD1, loc)

        
    # def test_get_practical_icrs_coordinates(self):

    #     def general_tests(Dec, Ra, PDec, PRa, accuracy):
    #         self.assertTrue((np.abs(Dec - PDec) <= accuracy).all())
    #         diff_Ra = np.abs(Ra - PRa).ravel()
    #         diff_Ra[diff_Ra > 6] = 2 * np.pi - diff_Ra[diff_Ra > 6]
    #         self.assertTrue((diff_Ra <= accuracy).all())
            
    #     x_fp, n_horns = st.get_full_fp('./ScanningTools/fp_data/fp_theta.txt', './ScanningTools/fp_data/fp_phi.txt')
    #     obs_time = (0, 2, 0, 0, 0)
    #     sampling_rate = 1
    #     rpm = 3
    #     day1, day2 = (2, None)
    #     n1, n2 = (48, None)
    #     obs_t1, time1, JD1 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day1)
    #     obs_t2, time2, JD2 = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day2)
    #     theta1, phi1, psi1 = st.get_engine_rotations(time1, rpm, zenith_distance,
    #                                                  polarization_angle)
    #     theta2, phi2, psi2 = st.get_engine_rotations(time2, rpm, zenith_distance,
    #                                                  polarization_angle)  
    #     fpp1 = st.get_fp_rotations(phi1, theta1, psi1, x_fp, n_horns, time1, n=n1, cartesian=False)
    #     fpp3 = st.get_fp_rotations(phi1, theta1, psi1, x_fp, n_horns, time1, n=n2, cartesian=False)
    #     Alt1, Az1 = st.get_horizon_coordinates(fpp1)
    #     Alt3, Az3 = st.get_horizon_coordinates(fpp3)
    #     Dec1, Ra1 = st.get_icrs_coordinates(JD1, loc, Alt1, Az1)      #1 day 2; n = 48 
    #     Dec3, Ra3 = st.get_icrs_coordinates(JD1, loc, Alt3, Az3)      #3 day 2; n = all
    #     PDec1, PRa1 = st.get_practical_icrs_coordinates(JD1, loc, Alt1, Az1)    #1 day 2; n = 48 
    #     PDec3, PRa3 = st.get_practical_icrs_coordinates(JD1, loc, Alt3, Az3)    #3 day 2; n = all
    #     accuracy = st.sex2dec(pointing_accuracy)
    #     general_tests(Dec1, Ra1, PDec1, PRa1, accuracy)
    #     general_tests(Dec3, Ra3, PDec3, PRa3, accuracy)

        
    # def test_get_polarization_angles(self):

    #     def general_tests(x_fp_pol_versors, pol_ang_proj, fp_pol_pointings):
    #         self.assertTrue((np.max(pol_ang_proj, axis=-1) <= np.pi).all())
    #         self.assertTrue((np.min(pol_ang_proj, axis=-1) >= -np.pi).all())
            
    #     zenith_distance, zenith_distance1 = (0, 10)
    #     boresight_angle = 0
    #     obs_time, obs_time1 = ((0, 0, 0, 1, 0), (0, 1, 0, 0, 0))        
    #     sampling_rate, sampling_rate1 = (50, 1)
    #     rpm, rpm1 = (1, 5)
    #     day, day1 = (None, 1)
    #     n, n1 = (0, None)
    #     obs_t, time, JD = st.get_timeJD(LCT_start, LCD_start, sampling_rate, obs_time, UTC=UTC,
    #                                        DST=DST, day=day)
    #     obs_t1, time1, JD1 = st.get_timeJD(LCT_start, LCD_start, sampling_rate1, obs_time1, UTC=UTC,
    #                                        DST=DST, day=day1)
    #     theta, phi, psi = st.get_engine_rotations(time, rpm, zenith_distance, boresight_angle)
    #     theta1, phi1, psi1 = st.get_engine_rotations(time1, rpm1, zenith_distance1, boresight_angle)
    #     theta2, phi2, psi2 = st.get_engine_rotations(time1, rpm1, zenith_distance, boresight_angle)
    #     x_fp_pol_angles, x_fp_pol_versors = st.get_full_fp_polarization_angles(
    #         './ScanningTools/fp_data/fp_psi.txt')
    #     n_horns = len(x_fp_pol_versors)
    #     fp_pol_pointings = st.get_fp_rotations(phi, theta, psi, x_fp_pol_versors, n_horns, time,
    #                                            n=n, cartesian=True) #rad
    #     fp_pol_pointings1 = st.get_fp_rotations(phi1, theta1, psi1, x_fp_pol_versors, n_horns,
    #                                             time1, n=n1, cartesian=True) #rad
    #     pol_ang_proj = st.get_polarization_angles(phi, theta, psi, x_fp_pol_versors, n_horns, time,
    #                                               n=n)
    #     pol_ang_proj1 = st.get_polarization_angles(phi1, theta1, psi1, x_fp_pol_versors, n_horns,
    #                                                time1, n=n1)
    #     pol_ang_proj2 = st.get_polarization_angles(phi2, theta2, psi2, x_fp_pol_versors, n_horns,
    #                                                time1, n=n1)
    #     pol_ang_proj_expected = np.concatenate((
    #         np.linspace(0, np.pi, sampling_rate * obs_t / 2 + 1),
    #         np.linspace(-np.pi, 0, sampling_rate * obs_t / 2 + 1)[1:-1]))
    #     self.assertTrue(np.allclose(np.arctan2(x_fp_pol_versors[..., 1], x_fp_pol_versors[..., 0]),
    #                                 x_fp_pol_angles))
    #     self.assertTrue(np.allclose(pol_ang_proj2[..., 0], x_fp_pol_angles))            
    #     self.assertTrue(np.allclose(pol_ang_proj, pol_ang_proj_expected))
    #     general_tests(x_fp_pol_versors, pol_ang_proj, fp_pol_pointings)
    #     general_tests(x_fp_pol_versors, pol_ang_proj1, fp_pol_pointings1)

        
    def test_get_scanning_strategy(self):

        def general_tests(packed_values):
            (x_fp, x_fp_pol_angles, n_horns, time, JD, theta, phi, psi, fp_pointings_spherical, Alt,
             Az, Dec, Ra, polarization_angles)  = packed_values
            self.assertTrue(np.allclose(x_fp.shape, (49, 3)))
            self.assertTrue(np.allclose(x_fp_pol_angles.shape, 49))
            self.assertEqual(n_horns, 49)
            self.assertTrue(np.allclose(time.shape, JD.shape))
            self.assertTrue(np.allclose(theta.shape, JD.shape))
            self.assertTrue(np.allclose(theta.shape, phi.shape))
            self.assertTrue(np.allclose(psi.shape, phi.shape))
            self.assertTrue(np.allclose(psi.shape, fp_pointings_spherical.shape[-2]))
            self.assertTrue(np.allclose(Alt.shape, fp_pointings_spherical.shape[:-1]))
            self.assertTrue(np.allclose(Alt.shape, Az.shape))
            self.assertTrue(np.allclose(Dec.shape, Az.shape))
            self.assertTrue(np.allclose(Dec.shape, Ra.shape))
            self.assertTrue(np.allclose(polarization_angles.shape, Ra.shape))
            self.assertTrue(np.allclose(time.shape, Ra.shape[-1]))
            self.assertTrue(np.allclose(fp_pointings_spherical.shape[-2], time.shape))

        obs_time = (0, 2, 0, 0, 0)
        sampling_rate = 2
        zenith_distance, polarization_angle = (10, 0)
        rpm = 5
        n1, n2 = (15, None)
        day1, day2 = (1, None)
        
        packed_values1 = st.get_scanning_strategy(
            obs_time, sampling_rate, zenith_distance, polarization_angle, rpm, n=n1, day=day2,
            LCT_start=(0, 0, 0), LCD_start=(1, 1, 2018), UTC=0, DST=0, LAT=np.array([28, 16, 24]),
            LONG=np.array([-16, 38, 32]), Height=2400,
            fp_theta_path='./ScanningTools/fp_data/fp_theta.txt',
            fp_phi_path='./ScanningTools/fp_data/fp_phi.txt',
            fp_psi_path='./ScanningTools/fp_data/fp_psi.txt')

        packed_values2 = st.get_scanning_strategy(
            obs_time, sampling_rate, zenith_distance, polarization_angle, rpm, n=n2, day=day1,
            LCT_start=(0, 0, 0), LCD_start=(1, 1, 2018), UTC=0, DST=0, LAT=np.array([28, 16, 24]),
            LONG=np.array([-16, 38, 32]), Height=2400,
            fp_theta_path='./ScanningTools/fp_data/fp_theta.txt',
            fp_phi_path='./ScanningTools/fp_data/fp_phi.txt',
            fp_psi_path='./ScanningTools/fp_data/fp_psi.txt')

        general_tests(packed_values1)
        general_tests(packed_values2)

        
if __name__ == '__main__':
    unittest.main()


