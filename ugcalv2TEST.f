C     ========================================================
C     Program ugcalv2ud to calculate upgraded tagger energy calibration.
C     ========================================================
C
C     27/3/13 copied from ugcalv2ub.f (same as a and c except name)
C     Added phenom corrn for 1.00-1.070T (840,855 and 883MeV) running
C     from D. Middleton's phenom corrns at 1.02517T see DM's
C     A2 internal Report 12/1 (4/2/13), fit coeffs from DM email 6/2/12 
C
C     20/10/10 Changed intercept in BCORR extrapolation below 0.5394498 T
C     to avoid a 'jump' at this field.
C
C     17/12/09 changed to neg slope in BCORR interpn 0.5394498 to 0.7545836 T
C     also 0.7545836 to 1.057T as Anal bk II p103, Ecal bkII p 64,58
C
C     'UNIVERSAL' version which interpolates BCORR between measured
C     values, extrapolates beyond, and applies phenom corrections in
C     ugcalv2c, ugcalv2c1895 and jugcalv2c for fields within +or- 1% 
C     of where these corrns were measured. Developed from ugcalv2c.f
C     8/12/09 bkII p 85, and Ecal analysis bk II p 101,103
C
C     Slight change to BCORR= 1.00032532 bkI p156 and
C     Tidied up messages to screen and to first output file 24/9/08
C     compared to prev version, now stored as ugcalv718208.f
C
C     Revised phenomenocogical correction frm Ecal anal p 171
C     Revised optimim BCORR for 1.834T, from Ecal analysis p 156
C     Phenomenological correction for 1.834 Tesla, 1508.0 MeV running
C     ugcalv2c.f jcm 14/2/08 
C     Optimum BCORR from April and Dec07 measts jcm 13/2/08 NOT YET A FN
C     of BNMR !!! DISABLED linear BCORR dependence on field.
C     Scint 33-39 widths changed to correct values ugcalv2.f jcm 11/12/07 
C     Reduced options version for CB users JCM 24/9/07 ugcalv1.f)
C     change to BCORR linear fn of field (NMR) 24/9/07 (ugcalv1.f
C     and ugcalv1gen.f) JCM
C     Change to BCORR 5/7/07 (ugcal.f) JCM
C
C     Changes for upgraded tagger JCM Feb 07 (modv7c.f)
C
C     Changes (from Ken Livingston's tagcalv6_cern) to get tagcalv7a
C     to compile on redhat  (not yet Scientific Linux !)
C     May 16,2006 done by JCM.
C
C     Ilya's changes to make this run on UNIX using CERN lib routine
C     done by JCM 18/8/00
C
C     JCM mod from GJM's TAGcalv7 to allow option of input of electron beam
C     energy (eg determined by MAMI people) rather than taking it from
C     the tagger field assuming perfect trajectory (including exit on axis)
C     through tagger. Input ebeamm, and get energies from ebeamm-... instead
C     of ebeam   14/6/99.
C
C     In this version the arrays [x_fpd],[y_fpd] contain the scintillator
C     positions (fixed) with respect to the "FPD" frame of reference (i.e.
C     a frame with x-axis along the rear straight edge of the Focal Plane
C     Detector). Three additional parameters (namely the FPD origin in the
C     "RADIATOR" frame of reference and the FPD x-axis gradient) are used to
C     transform the scintillator positions to the standard "RADIATOR" frame.
C
C     This version allows the actual (rather than nominal) scintillator widths
C     to define the energy bite of the overlapping scintillator pairs (and
C     triples).
C
      PROGRAM ugcalv2ud
C
      IMPLICIT NONE
C
      INTEGER NDIM,NP,JNP
      PARAMETER (NDIM=353,NP=353)
C
      REAL*4 BMAX
      PARAMETER (BMAX=1.4300)
C
      CHARACTER UPPER_CASE*40
      EXTERNAL UPPER_CASE
C
      REAL*4 X_FPD(NDIM),Y_FPD(NDIM),X(NDIM),Y(NDIM),P(NDIM),ANG(NDIM)
      REAL*4 DELE,KEMEAN,GRAD,X0,Y0,GRAD_DEF,X0_DEF,Y0_DEF,ANGLE
      REAL*4 PMEAN,DELP,theta_fpd(ndim),width(ndim),p1(ndim),p2(ndim)
      real*4 eg1_hi,eg1_lo,eg1_mean,eg2_hi,eg2_lo,eg2_mean
      real*4 eg_hi,eg_lo,eg_mean,ee_mean,dele1,dele2
      real*4 ang1(ndim),ang2(ndim),arm,pch,corrn
C
C     Added parameters for making BCORR a linear fn of NMR jcmmod 24/9/07
      real*4 bcorrm,bcorrc
C     end jcmmod
C
C     jcm mod 8/12/09 Add logical phenom to alter first output file header
      logical scint_plot,traj_plot,phenom			! diagnostic plots
C     end jcm mod 8/12/09
C      data    scint_plot,traj_plot/.false.,.false./ ! shifted to later jcm
C      data    arm/2000.0/                           ! 18/8/00
C
CCCCCC   jcm mod 14/6/99 add ebeamm
      REAL*4 ebeamm,EBEAM,BFIELD,SFAC,ANGI
CCCCCC   end jcm mod
      REAL*4 KEBEAM,PBEAM
C
      INTEGER I,J,K,L,IECL
      INTEGER EOF,iopt
C
      CHARACTER FILI*30,FILO*30,FILC*30,ANS*1
      CHARACTER HEAD(2,2,4)*40
C     jcm mod to set phenom false 8/12/09
      data    scint_plot,traj_plot,phenom/.false.,.false.,.false./
C     end jcm mod 8/12/09
      data    arm/2000.0/
C
      DATA HEAD/
     *' Goosy  Mean E_e  Bite   Mean E_g   ECL ',
     *'  HV  Station                           ',
     *' Chan.    (MeV)   (MeV)   (MeV)     out ',
     *'      Number                            ',
     *'  Goosy ch.   E_g_Low  E_g_High  E_g_Mea',
     *'n  E_g_bite    ECL   HV  Station        ',
     *'   firing      (MeV)     (MeV)     (MeV)',
     *'     (MeV)     out       Number         ',
     *'  Goosy ch.   E_g_Low  E_g_High  E_g_Mea',
     *'n  E_g_bite    ECL   HV  Station        ',
     *'   firing      (MeV)     (MeV)     (MeV)',
     *'     (MeV)     out       Number         ',
     *'                                        ',
     *'                                        ',
     *'                                        ',
     *'                                        '/
C
      character ecl(352)*1,hv(352)*2,sn(352)*2
      data ecl/16*'A',16*'B',16*'C',16*'D',16*'E',16*'F',16*'G',
     *         16*'H',16*'I',16*'J',16*'K',16*'L',16*'M',16*'N',
     *  16*'O',16*'P',16*'Q',16*'R',16*'S',16*'T',16*'U',16*'V'/
      data hv/'A','F','A','F','A','F','A','F','A','F','A','F',
     *        'A','F','A','F','A','F','A','F','A','F','A','F',
     *         56*'A',68*'B',68*'C',68*'D',36*'E',32*'E'''/
      data sn/'L','U','L','U','L','U','L','U','L','U','L','U',
     *        'L','U','L','U','L','U','L','U','L','U','L','U',
     *        328*'L'/
C
C     Double precision variables for NAG routine
      EXTERNAL RAYTRAK
      REAL*8 RELACC,ABSACC,PMIN,PMAX,PCAL,DELX
      INTEGER MAXCAL,IFAIL
      logical llm     ! Ilya's used for CERN function  JCM 18/8/00
C
      REAL*4 M0,KZ
C
C     Now the common parameters for the trajectory calculation
      REAL*4 EB,E0,S0X,W0,B0,A,LQ,SD,RHOB,R1,R2,PHIB,
     *ALPH1B,ALPH2B,PHI0,ALPH20,R2B
      COMMON/TPARMS/EB,E0,S0X,W0,B0,A,LQ,SD,RHOB,R1,R2,PHIB,
     *ALPH1B,ALPH2B,PHI0,ALPH20,R2B
C
      REAL*4 PB,P0,RHO0
      REAL*4 BHI,BLO,B,BCORR,BNMR
C
C     Next the common parameters for the minimisation calculation
      REAL*4 B1,X10,Y10,SIX,XSCINT,YSCINT,PHI,RHO
      INTEGER IBEND
      COMMON/MPARMS/B1,X10,Y10,SIX,XSCINT,YSCINT,PHI,RHO,IBEND
C     Default parameters for Mainz tagging spectrometer
C     assuming a 855.00Mev main beam energy
C
C
      DATA EB/855.511/,E0/458.311/,S0X/100.0/,W0/0.5/,B0/-4.4004/,A/5.3/,
C
C     JCM mod of SD,RHOB,R2 8/2/07 bk II p19, see also p39 and bkVI p 83 etc
C     but exact values not too important as BCORR is adjusted to best
C     fit the measurements.
C     *LQ/138.0/,SD/258.2/,RHOB/2800.0/,R1/183.8/,R2/-7995.0/,
     *LQ/138.0/,SD/283.0/,RHOB/2759.599464/,R1/183.8/,R2/-8021.0/,
     *PHIB/78.7636/,
C
C     JCM MOD of PHI0,ALPH20 8/2/07
C
C     *ALPH1B/16.7/,ALPH2B/-11.2364/,PHI0/75.43880/,ALPH20/-50.65460/,
     *ALPH1B/16.7/,ALPH2B/-11.2364/,PHI0/74.78526/,ALPH20/-51.16861/,
C     END JCM MOD
     *R2B/0.0/
C     Correction factor for non-homogeneous field profile. True field integral
C     = BCORR x homogeneous field integral (based on NMR value).
C  JCM mod of BCORR 24/9/08, 23/10/06, 12/2/07, 13/2/07, 5/7/07
C      DATA BCORR/1.00262/
C      DATA BCORR/1.01087/
C      DATA BCORR/1.00623/
C     Value for 1.508 GeV (1.834 T) deduced from first look at
C     April 07 calibn measts  - see 5/7/07 sheet
C      DATA BCORR/1.001299/
C     BCORR from combined April and Dec 07 measts 13/2/08 bkII p 145
C     and Ecal analysis bk p156 jcm
C      DATA BCORR/1.000295211/
C
C     change to value derived AFTER correcting Emami to 855.335 MeV in measts
C     bk I p 156
      DATA BCORR/1.00032532/
C
C     jcmmod to make bcorr a lin fn. of NMR
      data bcorrm/-0.010876/,bcorrc/1.021246/
C     END JCM MOD
C
C     ELECTRON REST ENERGY
      DATA M0/0.511/,KZ/29.979613E-02/
C
      REAL*4 RAD,ARG,MOM,ENERGY
C     DEFINE FUNCTION TO CONVERT FROM DEGREES TO RADIANS
C      RAD(ARG)=ARG*3.1415926/180.0
C     DEFINE FUNCTION TO CALULATE MOMENTUM FROM ENERGY
      MOM(ARG)=SQRT(ARG**2-0.511**2)
      ENERGY(ARG)=SQRT(ARG**2+0.511**2)
C
C     JCM mod 16/5/06
C     STUFF TO GET TRIG FUNCTIONS without using absoft (from Ken)
C
C
      REAL*4 COSD,SIND,ACOSD,ASIND,ATAND,THETAD,XD
      COSD(THETAD)=COS(THETAD*3.1415927/180.)
      SIND(THETAD)=SIN(THETAD*3.1415927/180.)
      ACOSD(XD)=180.*ACOS(XD)/3.1415927
      ASIND(XD)=180.*ASIN(XD)/3.1415927
      ATAND(XD)=180.*ATAN(XD)/3.1415927
C
C     end JCM mod
C
C     Next the coord locations of the scintillators. These values are
C     those measured in the "FPD" frame of reference (x-axis along the
C     rear straight edge of the FPD) and are fixed forever wrt this frame.
C     The scintillator positions in the "RADIATOR" frame (ie wrt the dipole,
C     radiator etc.) depend upon the FPD orientation within that frame - in
C     this case defined by the gradient of the FPD x-axis (GRAD) and the
C     position of the FPD origin (x0,y0) in the radiator frame.
C
      data (x_fpd(i),i=  1,120)/
     *  247.9970, 236.3300, 224.6630, 212.9960, 201.3290, 189.6600,
     *  177.9440, 166.1560, 154.2950, 142.3590, 130.3460, 118.2530,
     *  106.0800, 93.82600, 81.49000, 69.07200, 56.57300, 43.99600,
     *  31.34200, 18.61700, 5.826000,-7.024000,-19.92500,-32.86600,
     * -45.83700,-58.82800,-71.82700,-84.82500,-97.81300,-110.7860,
     * -123.7400,-136.6760,-149.5960,-162.5060,-175.4110,-188.3190,
     * -201.2370,-214.1700,-227.1140,-240.0650,-253.0210,-265.9820,
     * -278.9480,-291.9190,-304.8930,-317.8710,-330.8520,-343.8360,
     * -356.8230,-369.8130,-382.8040,-395.7970,-408.7910,-421.7870,
     * -434.7840,-447.7820,-460.7810,-473.7800,-486.7790,-499.7790,
     * -512.7790,-525.7790,-538.7790,-551.7790,-564.7780,-577.7770,
     * -590.7750,-603.7730,-616.7700,-629.7660,-642.7610,-655.7550,
     * -668.7490,-681.7420,-694.7330,-707.7240,-720.7130,-733.7010,
     * -746.6880,-759.6740,-772.6590,-785.6420,-798.6250,-811.6060,
     * -824.5850,-837.5640,-850.5410,-863.5180,-876.4930,-889.4670,
     * -902.4390,-915.4110,-928.3810,-941.3500,-954.3170,-967.2840,
     * -980.2500,-993.2140,-1006.178,-1019.140,-1032.101,-1045.061,
     * -1058.021,-1070.979,-1083.936,-1096.893,-1109.848,-1122.803,
     * -1135.757,-1148.710,-1161.663,-1174.614,-1187.565,-1200.515,
     * -1213.465,-1226.413,-1239.361,-1252.309,-1265.256,-1278.203/
      data (x_fpd(i),i=121,240)/
     * -1291.149,-1304.095,-1317.040,-1329.984,-1342.929,-1355.873,
     * -1368.817,-1381.761,-1394.704,-1407.647,-1420.590,-1433.532,
     * -1446.475,-1459.417,-1472.360,-1485.302,-1498.245,-1511.187,
     * -1524.129,-1537.071,-1550.014,-1562.957,-1575.900,-1588.842,
     * -1601.785,-1614.729,-1627.672,-1640.616,-1653.560,-1666.504,
     * -1679.448,-1692.393,-1705.338,-1718.284,-1731.230,-1744.177,
     * -1757.124,-1770.071,-1783.018,-1795.967,-1808.916,-1821.865,
     * -1834.815,-1847.765,-1860.716,-1873.667,-1886.620,-1899.572,
     * -1912.525,-1925.480,-1938.434,-1951.389,-1964.345,-1977.301,
     * -1990.258,-2003.216,-2016.175,-2029.134,-2042.094,-2055.055,
     * -2068.016,-2080.978,-2093.940,-2106.904,-2119.868,-2132.833,
     * -2145.799,-2158.765,-2171.732,-2184.701,-2197.669,-2210.639,
     * -2223.609,-2236.580,-2249.551,-2262.523,-2275.498,-2288.472,
     * -2301.447,-2314.422,-2327.399,-2340.376,-2353.354,-2366.333,
     * -2379.312,-2392.292,-2405.273,-2418.255,-2431.237,-2444.220,
     * -2457.204,-2470.188,-2483.173,-2496.159,-2509.145,-2522.132,
     * -2535.120,-2548.107,-2561.096,-2574.086,-2587.076,-2600.066,
     * -2613.057,-2626.049,-2639.041,-2652.033,-2665.027,-2678.020,
     * -2691.014,-2704.009,-2717.003,-2729.999,-2742.994,-2755.991,
     * -2768.987,-2781.984,-2794.982,-2807.979,-2820.977,-2833.975/
      data (x_fpd(i),i=241,353)/
     * -2846.973,-2859.971,-2872.970,-2885.969,-2898.968,-2911.967,
     * -2924.966,-2937.966,-2950.966,-2963.966,-2976.965,-2989.965,
     * -3002.965,-3015.966,-3028.966,-3041.966,-3054.966,-3067.966,
     * -3080.965,-3093.965,-3106.964,-3119.963,-3132.962,-3145.960,
     * -3158.958,-3171.956,-3184.954,-3197.952,-3210.949,-3223.945,
     * -3236.942,-3249.937,-3262.932,-3275.927,-3288.921,-3301.914,
     * -3314.906,-3327.898,-3340.889,-3353.879,-3366.869,-3379.858,
     * -3392.846,-3405.834,-3418.821,-3431.807,-3444.791,-3457.776,
     * -3470.759,-3483.740,-3496.721,-3509.701,-3522.680,-3535.657,
     * -3548.634,-3561.609,-3574.583,-3587.555,-3600.527,-3613.497,
     * -3626.466,-3639.433,-3652.399,-3665.364,-3678.326,-3691.288,
     * -3704.248,-3717.207,-3730.162,-3743.118,-3756.071,-3769.023,
     * -3781.973,-3794.920,-3807.866,-3820.811,-3833.753,-3846.693,
     * -3859.632,-3872.568,-3885.503,-3898.434,-3911.365,-3924.292,
     * -3937.218,-3950.141,-3963.063,-3975.981,-3988.898,-4001.812,
     * -4014.724,-4027.632,-4040.539,-4053.443,-4066.345,-4079.244,
     * -4092.140,-4105.034,-4117.925,-4130.813,-4143.698,-4156.581,
     * -4169.460,-4182.336,-4195.210,-4208.081,-4220.948,-4233.813,
     * -4246.675,-4259.534,-4272.390,-4285.241,-4298.091/
      data (y_fpd(i),i=  1,120)/
     * 6.714000,12.44800,18.18200,23.91700,29.65100,35.38200,
     * 41.01500,46.49700,51.81900,56.96900,61.93600,66.70800,
     * 71.27100,75.61000,79.71000,83.55700,87.13100,90.41700,
     * 93.39700,96.05500,98.37500,100.3420,101.9450,103.1770,
     * 104.0330,104.5170,104.6370,104.4120,103.8640,103.0270,
     * 101.9400,100.6510,99.21200,97.67900,96.11200,94.56900,
     * 93.11200,91.79700,90.59300,89.45900,88.38900,87.38700,
     * 86.44800,85.57400,84.75800,84.00200,83.30400,82.66300,
     * 82.07800,81.54500,81.06700,80.64100,80.26300,79.93600,
     * 79.65600,79.42200,79.23500,79.09200,78.99300,78.93700,
     * 78.92200,78.94700,79.01100,79.11600,79.25600,79.43400,
     * 79.64800,79.89700,80.17900,80.49400,80.84300,81.22300,
     * 81.63400,82.07400,82.54400,83.04200,83.56800,84.12100,
     * 84.70000,85.30500,85.93600,86.59000,87.26700,87.96800,
     * 88.69100,89.43600,90.20200,90.98800,91.79500,92.62000,
     * 93.46400,94.32700,95.20700,96.10400,97.01800,97.94800,
     * 98.89300,99.85400,100.8280,101.8180,102.8200,103.8350,
     * 104.8630,105.9040,106.9560,108.0190,109.0930,110.1770,
     * 111.2720,112.3750,113.4890,114.6100,115.7400,116.8780,
     * 118.0230,119.1760,120.3360,121.5010,122.6740,123.8510/
      data (y_fpd(i),i=121,240)/
     * 125.0340,126.2210,127.4130,128.6100,129.8110,131.0150,
     * 132.2240,133.4340,134.6470,135.8620,137.0800,138.2990,
     * 139.5200,140.7420,141.9650,143.1880,144.4110,145.6340,
     * 146.8590,148.0800,149.3020,150.5230,151.7420,152.9600,
     * 154.1750,155.3900,156.6010,157.8090,159.0160,160.2180,
     * 161.4180,162.6140,163.8050,164.9930,166.1760,167.3550,
     * 168.5300,169.7000,170.8630,172.0230,173.1760,174.3240,
     * 175.4650,176.6000,177.7290,178.8520,179.9680,181.0760,
     * 182.1780,183.2730,184.3590,185.4380,186.5090,187.5710,
     * 188.6270,189.6720,190.7100,191.7390,192.7580,193.7700,
     * 194.7710,195.7630,196.7450,197.7180,198.6810,199.6330,
     * 200.5760,201.5080,202.4300,203.3400,204.2400,205.1300,
     * 206.0070,206.8740,207.7290,208.5730,209.4050,210.2250,
     * 211.0340,211.8300,212.6130,213.3850,214.1440,214.8910,
     * 215.6250,216.3460,217.0530,217.7480,218.4290,219.0980,
     * 219.7530,220.3940,221.0220,221.6360,222.2360,222.8220,
     * 223.3940,223.9510,224.4940,225.0230,225.5370,226.0370,
     * 226.5220,226.9920,227.4470,227.8860,228.3110,228.7210,
     * 229.1150,229.4940,229.8570,230.2040,230.5370,230.8530,
     * 231.1530,231.4370,231.7050,231.9570,232.1930,232.4130/
      data (y_fpd(i),i=241,353)/
     * 232.6160,232.8020,232.9720,233.1260,233.2620,233.3820,
     * 233.4850,233.5720,233.6400,233.6920,233.7270,233.7440,
     * 233.7450,233.7270,233.6930,233.6400,233.5700,233.4820,
     * 233.3770,233.2540,233.1130,232.9540,232.7770,232.5820,
     * 232.3690,232.1390,231.8880,231.6210,231.3350,231.0300,
     * 230.7070,230.3650,230.0050,229.6270,229.2290,228.8120,
     * 228.3780,227.9240,227.4500,226.9580,226.4480,225.9180,
     * 225.3690,224.8000,224.2130,223.6060,222.9810,222.3350,
     * 221.6700,220.9860,220.2820,219.5590,218.8160,218.0540,
     * 217.2720,216.4700,215.6480,214.8070,213.9460,213.0640,
     * 212.1640,211.2440,210.3020,209.3410,208.3600,207.3600,
     * 206.3380,205.2970,204.2350,203.1540,202.0520,200.9290,
     * 199.7860,198.6230,197.4400,196.2370,195.0120,193.7670,
     * 192.5020,191.2160,189.9100,188.5830,187.2340,185.8660,
     * 184.4770,183.0670,181.6370,180.1860,178.7130,177.2210,
     * 175.7070,174.1730,172.6170,171.0410,169.4430,167.8250,
     * 166.1860,164.5250,162.8440,161.1420,159.4180,157.6730,
     * 155.9080,154.1220,152.3140,150.4840,148.6340,146.7620,
     * 144.8700,142.9560,141.0210,139.0640,137.0870/
C
C     Next the locating angle (in degrees) of the scintillators. These values
C     are those measured in the "FPD" frame of reference (x-axis along the
C     rear straight edge of the FPD) and are fixed forever wrt this frame.
C
      data (theta_fpd(i),i=  1,120)/
     *  133.6296,133.6290,133.6291,133.6291,133.6293,133.8052,
     *  134.4877,135.2440,136.0448,136.8943,137.7975,138.7551,
     *  139.7727,140.8513,141.9976,143.2102,144.4933,145.8470,
     *  147.2704,148.7554,150.3055,149.3918,148.5176,147.6798,
     *  146.8821,146.1218,145.3988,144.7126,144.0597,143.4384,
     *  142.8475,142.2822,141.7398,141.2183,140.7134,140.2221,
     *  139.7415,139.2690,138.8053,138.3518,137.9087,137.4748,
     *  137.0500,136.6333,136.2258,135.8265,135.4351,135.0516,
     *  134.6757,134.3066,133.9453,133.5906,133.2427,132.9011,
     *  132.5664,132.2371,131.9140,131.5971,131.2853,130.9799,
     *  130.6790,130.3843,130.0941,129.8094,129.5291,129.2540,
     *  128.9835,128.7176,128.4562,128.1993,127.9468,127.6982,
     *  127.4540,127.2136,126.9775,126.7448,126.5162,126.2914,
     *  126.0700,125.8522,125.6381,125.4272,125.2202,125.0160,
     *  124.8152,124.6175,124.4230,124.2315,124.0434,123.8579,
     *  123.6755,123.4961,123.3193,123.1454,122.9741,122.8060,
     *  122.6400,122.4766,122.3160,122.1579,122.0023,121.8489,
     *  121.6980,121.5498,121.4036,121.2597,121.1182,120.9790,
     *  120.8418,120.7068,120.5742,120.4434,120.3148,120.1883,
     *  120.0636,119.9411,119.8205,119.7018,119.5852,119.4702/
      data (theta_fpd(i),i=121,240)/
     *  119.3573,119.2461,119.1367,119.0292,118.9234,118.8194,
     *  118.7168,118.6164,118.5175,118.4201,118.3244,118.2304,
     *  118.1379,118.0470,117.9577,117.8699,117.7836,117.6989,
     *  117.6156,117.5337,117.4533,117.3744,117.2969,117.2207,
     *  117.1460,117.0725,117.0006,116.9299,116.8606,116.7925,
     *  116.7257,116.6602,116.5960,116.5330,116.4713,116.4107,
     *  116.3515,116.2934,116.2364,116.1807,116.1261,116.0727,
     *  116.0203,115.9692,115.9192,115.8701,115.8223,115.7755,
     *  115.7298,115.6850,115.6414,115.5989,115.5574,115.5168,
     *  115.4773,115.4387,115.4012,115.3646,115.3290,115.2944,
     *  115.2607,115.2279,115.1961,115.1652,115.1352,115.1062,
     *  115.0780,115.0507,115.0243,114.9987,114.9740,114.9503,
     *  114.9273,114.9052,114.8838,114.8634,114.8437,114.8249,
     *  114.8068,114.7896,114.7731,114.7575,114.7426,114.7284,
     *  114.7151,114.7025,114.6906,114.6794,114.6691,114.6594,
     *  114.6505,114.6423,114.6348,114.6280,114.6218,114.6164,
     *  114.6117,114.6077,114.6043,114.6016,114.5996,114.5982,
     *  114.5975,114.5974,114.5980,114.5992,114.6010,114.6035,
     *  114.6065,114.6103,114.6146,114.6195,114.6250,114.6311,
     *  114.6378,114.6451,114.6530,114.6615,114.6705,114.6801/
      data (theta_fpd(i),i=241,353)/
     *  114.6903,114.7009,114.7123,114.7241,114.7364,114.7494,
     *  114.7628,114.7768,114.7913,114.8064,114.8220,114.8380,
     *  114.8546,114.8717,114.8893,114.9074,114.9260,114.9451,
     *  114.9647,114.9848,115.0054,115.0265,115.0480,115.0700,
     *  115.0924,115.1154,115.1387,115.1626,115.1869,115.2116,
     *  115.2368,115.2625,115.2886,115.3151,115.3421,115.3695,
     *  115.3973,115.4255,115.4542,115.4833,115.5128,115.5427,
     *  115.5731,115.6038,115.6349,115.6665,115.6984,115.7308,
     *  115.7635,115.7966,115.8301,115.8640,115.8983,115.9329,
     *  115.9680,116.0033,116.0391,116.0753,116.1118,116.1487,
     *  116.1859,116.2235,116.2614,116.2998,116.3384,116.3774,
     *  116.4168,116.4565,116.4965,116.5368,116.5776,116.6186,
     *  116.6600,116.7017,116.7438,116.7861,116.8288,116.8718,
     *  116.9151,116.9588,117.0027,117.0470,117.0915,117.1364,
     *  117.1816,117.2271,117.2729,117.3190,117.3654,117.4121,
     *  117.4591,117.5064,117.5539,117.6018,117.6499,117.6984,
     *  117.7471,117.7961,117.8454,117.8949,117.9447,117.9948,
     *  118.0452,118.0958,118.1467,118.1979,118.2493,118.3010,
     *  118.3530,118.4052,118.4577,118.5105,118.5635/
C
C     Next the actual scintillator widths (in mm) employed in the
C     FPD construction. Note how the widths are batched where the
C     nominal variation is slow.
C
C     jcm mod 11/12/07 to correct widths of scints 33-39
C
      data (width(i),i=1,353)/
     *  32.607769,32.221680,31.852831,31.501040,31.164835,30.803890,
     *  30.347589,29.898945,29.464128,29.040533,28.626196,28.220430,
     *  27.821032,27.426767,27.035423,26.647329,26.259811,25.872179,
     *  25.483376,25.093826,24.701313,24.284260,23.821457,23.317320,
     *  22.778624,22.215502,21.637844,21.062616,20.501245,19.973793,
C     *  19.492910,19.075768,17*18.21,8*17.75,11*17.38,9*16.95,9*16.60,
     *  19.492910,19.075768,7*18.70,10*18.21,8*17.75,11*17.38,9*16.95,
     * 9*16.60,
     *  9*16.25,9*15.95,10*15.60,10*15.30,10*14.95,11*14.60,11*14.25,
     *  11*13.90,11*13.55,11*13.25,12*12.90,12*12.60,12*12.25,12*11.90,
     *  12*11.60,12*11.30,12*11.00,12*10.75,12*10.45,12*10.20,13*9.95,
     *  13*9.70,13*9.45,5*9.20/
C
C     Default values for FPD origin (in mm) and x-axis gradient.
C
      data x0_def,y0_def,grad_def/-707.134,616.522,-0.825231/
C
C
C     JCM mods format for 1500 MeV 16/5/06
C
C2005      format(3x,i3,3x,f7.3,3x,f5.3,3x,f7.3,3x,a1,i2,3x,a2,3x,a2,i3)
2005      format(3x,i3,3x,f8.3,3x,f8.3,3x,f8.3,3x,a1,i2,3x,a2,3x,a2,i3)
C
C2006      format(2x,i3,'  only',4x,f7.3,2x,f7.3,4x,f7.3,3x,f6.3,
C     1        5x,a1,i2,3x,a2,3x,a2,i3)
C	1        5x,a1,i2,3x,a2,3x,a2,i3) ! jcm
2006      format(2x,i3,'  only',4x,f8.3,2x,f8.3,4x,f8.3,3x,f8.3,
     1        5x,a1,i2,3x,a2,3x,a2,i3)
C
C2007      format(2x,i3,' & ',i3,4x,f7.3,2x,f7.3,4x,f7.3,3x,f6.3)
2007      format(2x,i3,' & ',i3,4x,f8.3,2x,f8.3,4x,f8.3,3x,f8.3)
C
C2008      format(2x,i3,' s+nd',5x,f7.3,2x,f7.3,4x,f7.3,3x,f6.3,
C     1        5x,a1,i2,3x,a2,3x,a2,i3)
C	1        5x,a1,i2,3x,a2,3x,a2,i3) ! jcm
2008      format(2x,i3,' s+nd',5x,f8.3,2x,f8.3,4x,f8.3,3x,f8.3,
     1        5x,a1,i2,3x,a2,3x,a2,i3)
C
C2009      format(2x,i3,' & ',i3,4x,f7.3,2x,f7.3,4x,f7.3,3x,f6.3,
C     1        ' < no geometric overlap')
C	1        ' < no geometric overlap') ! jcm
2009      format(2x,i3,' & ',i3,4x,f8.3,2x,f8.3,4x,f8.3,3x,f8.3,
     1        ' < no geometric overlap')
C
C     JCM MOD to fix bug in option 5 - should say 's+bnd'
C     rather than 'only' as in format 2006, so define format 2010
C
2010     format(2x,i3,' s+bnd',4x,f8.3,2x,f8.3,4x,f8.3,3x,f8.3,
     1        5x,a1,i2,3x,a2,3x,a2,i3)
C
C          end JCM mods
C
      WRITE(*,1011)
C     JCM mod 8/12/09 suppress write of BCORR to screen until later
C      WRITE (*,*) BCORR
1011  FORMAT(
     *'           Upgraded Tagger Calibration ugcalv2ud',/,
     *'     J.C.McGeorge 20/10/10,I.Anthony 7/2/92 G.J.Miller',//,
     *' Calculates energy calibration of the upgraded',/
     * 'Glasgow Mainz tagging spectrometer.',//,
     * 'USES INTERPOLATED FIELD CORRECTION FACTOR, BCORR',/,
     * 'UNKNOWN UNCERTAINTY AT FIELDS LESS THAN 0.5 TESLA'/
     * 'PHENOMEOLOGICAL ENERGY CORRECTION for effect of large-scale',/,
     * 'field non-uniformity derived from calibration measurements',/,
     * 'made April and December 07 BUT ONLY FOR fields within',/,
     * ' 1 percent of 1.89563, 1.834 and 1.443 Tesla where determined',//,
C     endjcm mod 8/12/09
     * ' ALSO for fields 1.000 - 1.070T from calibn measts July and Oct 2011: ',/,
     * ' see D. Middleton A2 Internal Report 2012/1',/,
C     *'     assuming scintillator centres to be limits of overlap.',/,
C     *' 2)- Photon energy limits for single and neighbouring double',/,
C     *'     o/p channels (separately). True scintillator geometry.',/,
     *' Outputs PHOTON energies & channel widths for ',/,
     *' single and neighbouring double channels combined',/,
     *' True scintillator geometry.',//,
     *' An output table in 5X3G format can also be produced suitable',/,
     *' for display by a graphics package',/,
     *' Input data required:-',/,
     *' Main beam (TOTAL) energy (MeV)',
     *' Tagger magnetic field from NMR (Tesla)',//,
     *' Assumes coordinate system origin at the radiator with +ve
     * Y axis in the',/,
     *' direction of the input electron beam. Scintillator locations
     * are ordered',/,
     *' from low to high electron momentum. Channels are then numbered
     * increasing',/,
     *' in this direction from channel 1 (first disc./coinc. card which
     * has a signal',/,
     *' output i.e. scint. station no. 2). All energies are TOTAL.',/
     *' ...............',/)
C
C     jcmmod 24/9/07 hard wired ANS to be option F
C      ANS=' '
      ANS='F'
C      DO WHILE (ANS.NE.'B'.AND.ANS.NE.'F')
C        WRITE(*,1008)
C1008    FORMAT(' Calculate calibration using:- main beam energy-(B) '
C     *  ,/,    '                               magnet field-(F) [A1]',$)
C        READ(*,1001)ANS
C        ANS=UPPER_CASE(ANS)
C      END DO
C     end jcmmod
      PB=MOM(EB)
C
C     jcm mod 24/9/07 removed questions not relevant in option F
C      IF(ANS.EQ.'B')THEN
C        WRITE(*,1006)
C1006    FORMAT(/,' Give desired main beam (TOTAL) energy (MeV) [G]:',$)
C
C       JCM mod format without absoft 16/5/06
C
C        READ(*,1007)EBEAM
C        READ(*,*)EBEAM
C
C       end JCM mod
C
C        PBEAM=MOM(EBEAM)
C        B1 = PBEAM/(KZ*RHOB) 		!Ideal field in magnet in tesla.
C        BNMR = B1/BCORR			!NMR field reading
C        SFAC=PBEAM/PB			!Scaling factor from standard setting
C
C      ELSE IF(ANS.EQ.'F')THEN
C     endjcm mod 24/9/07
        WRITE(*,1009)
1009    FORMAT(/,' Give magnet field NMR reading (Tesla) [G]:',$)
C
C
C     JCM mod 16/5/06
C     STUFF TO GET format without absoft
C
C        READ(*,1007)BNMR		!NMR field in magnet
        READ(*,*)BNMR		!NMR field in magnet
C
C     jcmmod 24/9/07 to make bcorr a linear fn of the NMR value
C     FOR OPTION F only !  TAKEN OUT FOR NOW jcm 14/2/08
C        bcorr=bcorrm*BNMR+bcorrc
C     end jcm mod
C
C     jcm mod 8/12/09 to interpolate BCORR value between measured points
C     ECAL analysis bkII p101,103 and bcorr1a chart 4
        IF(BNMR.GE.1.89563) BCORR=-0.014983223*BNMR+1.027677707
        IF((BNMR.GE.1.83400).AND.(BNMR.LT.1.89563)) BCORR=-0.017041392*BNMR+1.031579233
        IF((BNMR.GE.1.57037).AND.(BNMR.LT.1.83400)) BCORR=-0.019783283*BNMR+1.036607862
        IF((BNMR.GE.1.44300).AND.(BNMR.LT.1.57037)) BCORR=-0.024104679*BNMR+1.043394052
        IF((BNMR.GE.1.199669).AND.(BNMR.LT.1.44300)) BCORR=-0.002335444*BNMR+1.011981046
        IF((BNMR.GE.1.05700).AND.(BNMR.LT.1.199669)) BCORR=0.0012075924*BNMR+1.007730575
        IF((BNMR.GE.0.7545836).AND.(BNMR.LT.1.05700)) BCORR=-0.000301227*BNMR+1.009325397
        IF((BNMR.GE.0.5394498).AND.(BNMR.LT.0.7545836)) BCORR=-0.0025503059*BNMR+1.011022515
C        IF(BNMR.LT.0.5394498) BCORR=-0.000683625*BNMR+1.009839621
C     slight change  - see below
C     end jcm mod 8/12/09
C
C     jcm mod 20/10/10 to avoid 'jump' at 0.5394498 tesla
C
        IF(BNMR.LT.0.5394498) BCORR=-0.000683625*BNMR+1.010015534
C     end jcm mod 20/10/10
C
        B1 = BCORR * BNMR		!Ideal magnet field requested
        BFIELD = MOM(EB)/(KZ*RHOB)	!Standard ideal field 
        SFAC=B1/BFIELD			!Scaling factor from standard setting
        PBEAM=PB*SFAC
        EBEAM=ENERGY(PBEAM)
C
C     jcmmod 24/9/07 remove this endif 
C      ENDIF
C     end jcmmod 24/9/07
C
C     JCM mod 16/5/06
C     STUFF TO GET format without absoft
C
C1007  FORMAT(G)
C
C      end of JCM mod
C
1001  FORMAT(A)
CCCCCCCCCCCCCCCCC  jcm mod 14/6/99 to choose beam energy from tagger or MAMI
C
C     jcmmod 24/9/07 to remove this option in ugcalv1
C
      ANS='M'
C      ANS=' '
C      DO WHILE (ANS.NE.'T'.AND.ANS.NE.'M')
C        WRITE(*,991)
C991    FORMAT(' Take electron beam energy from:- tagger -(T) '
C     *  ,/,    '                               MAMI -(M) [A1]',$)
C        READ(*,1001)ANS
C        ANS=UPPER_CASE(ANS)
C      END DO
C      If(ans.eq.'T') then
C      ebeamm=ebeam
C      else if(ans.eq.'M') then
C     end jcmmod 24/9/07
C
      Write(*,992)
992    format(/,' Enter TOTAL (include rest mass) e beam
     *energy from MAMI ')
C
C     JCM mod - format without absoft 16/5/06
C
C      read(*,1007) ebeamm
      read(*,*) ebeamm
C     end JCM mod
C
C     jcmmod 24/9/07 remove this endif
C      endif
C     end jcmmod 24/9/07
C
      write(*,993) ebeamm
C
C     JCM mod format without absoft 16/5/06
C
C993    format(/,' Using total e beam energy = ', G)
993    format(/,' Using total e beam energy = ', F8.3)
C
C     END JCM MOD
C
C     jcm mod to 'hard wire' default FPD location 8/12/09
C1030  format(/,' Enter position of FPD origin (in mm) and gradient of',
C     *         ' FPD x-axis',/,' wrt radiator frame (X0,Y0,GRAD).',
C     *         ' 0.0,0.0,0.0 for default values : ')
1030  format(/,' Using standard FPD location origin (in mm) and gradient  : ')
C     end of jcm mod 8/12/09
C
C
C     JCM mod 16/5/06
C     STUFF TO GET format without absoft
C
C1031  format(3g)
C
C      end jcm mod
C
1032  format(' X0 = ',f9.3,' , Y0 = ',f9.3,' , GRAD = ',f9.6,/)
      write (*,1030)
C
C     JCM mod 16/5/06
C     STUFF TO GET format without absoft
C
C      read(*,1031)x0,y0,grad
C     jcm mod to suppress read in of FPD location and hard wire it8/12/09
C      read(*,*)x0,y0,grad
C
C     end jcm mod
C     suppress if and endif to hard wire default FPD location 8/12/09
C      if(x0.eq.0.0.and.y0.eq.0.0.and.grad.eq.0.0)then	! Use default values.
        x0=x0_def
        y0=y0_def
        grad=grad_def
        write(*,1032)x0,y0,grad
C      end if
C     end jcm mod 8/12/09
C
C
C     Scintillator positions in "RADIATOR" frame of reference.
C
      angle=atand(grad)
      do i=1,np
        x(i)=x_fpd(i)*cosd(angle)-y_fpd(i)*sind(angle)+x0
        y(i)=x_fpd(i)*sind(angle)+y_fpd(i)*cosd(angle)+y0
      end do
C
C
C     jcm mod 24/9/07 output option is hard wired, removed this list
C      WRITE(*,1012)
C1012  FORMAT(' Output table of:-',//,
C     *' (1)- Electron energy, photon energy and bite per channel',/,
C     *'      assuming scintillator centres to be limits of overlap.',//,
C     *' (2)- Photon energy limits for single and neighbouring double',/,
C     *'      o/p channels (separately). True scintillator geometry.',//,
C     *' (3)- Photon energy limits for single and neighbouring double',/,
C     *'      o/p channels (combined). True scintillator geometry.',//,
C     *' (4)- Photon energy limits for singles only',/,
C     *'      True scintillator geometry.',//,
C     *' (5)- Photon energy limits for single and both neighbouring',/,
C     *'      doubles. True scintillator geometry.',//,
C     *' Select option:',$)
C
C
C
C     JCM mod 16/5/06
C     STUFF TO GET format without absoft
C
C      READ(*,1013) IOPT
C      READ(*,*) IOPT
      IOPT=5
C     end jcm mod 24/9/07
C1013  FORMAT(I)
C      end JCM mod
C
      WRITE(*,1004)
1004  FORMAT(' Give filename for printed output [A30]:',$)
      READ(*,1001)FILO
      OPEN(UNIT=2,FILE=FILO,STATUS='NEW',FORM='FORMATTED',
     *ACCESS='SEQUENTIAL')
      WRITE(*,1014)
1014  FORMAT(' Give output filename for plot data [A30]
     * <cr> = none):',$)
      READ(*,1001)FILC
      IF(FILC.NE.' ')OPEN(UNIT=3,FILE=FILC,STATUS='NEW',
     *FORM='FORMATTED',ACCESS='SEQUENTIAL')
C
      if(scint_plot)then
       open(unit=7,file='scint',status='new',form='formatted',
     * access='sequential')
      end if
      if(traj_plot)then
       open(unit=8,file='traj',status='new',form='formatted',
     * access='sequential')
      end if
C
C
C     JCM mod 16/5/06
C     STUFF TO GET format without absoft
C
C
C1003  FORMAT(2G)
C
C     end JCM mod
      IBEND=-1
C     JCM MOD TO WRITE BCORR VALUE TO SCREEN
      WRITE(*,10001) BCORR
10001 FORMAT(' Using BCORR = ', F11.8)
C
C********************************************************************
C
C     FIRST SECTION CALCULATES MAIN BEAM TRAJECTORY
C
C********************************************************************
C
      SIX = S0X + LQ + SD
C
C
C
C******************************************************************
C
C     SECOND SECTION CALCULATES CENTRAL TRAJECTORY
C
C******************************************************************
C
C
C     CENTRAL RAY  MOMENTUM
      P0=SFAC*MOM(E0)
C
C     RHO0 IS RADIUS OF CENTRAL RAY IN MAGNET
      RHO0 = P0/(KZ*B1)
C
C     CALULATE INTERSECTION POINTS OF CENTRAL RAY WITH MAGNET FACE
C     EXIT FACE INTERSECTION (X10,Y10)
 
      X10 = IBEND * (RHO0 - RHO0*COSD(PHI0))
      Y10 = SIX  + RHO0*SIND(PHI0)
C
C     JCM MOD type out X10,Y10 for tests 8/2/07
C
      write(*,10151) X10,Y10,RHO0
10151 FORMAT(' (X10,Y10) RHO0 = ', 3F11.5)
C     END JCM MOD
C
C     JCM MOD ebeam changed to ebeamm everywhere from here on
C
C     jcm mod 8/12/09 to header wording and to include BCORR
C
      WRITE(2,1015)x0,y0,grad,EBEAMm,PBEAM
1015  FORMAT(' ugcalv2ud.f : Used interpolated field correction factor, BCORR ',/,
     *' (Uncertainty is unknown for fields less than 0.5 Tesla)',/,
     *' Channel widths include both neighbouring double hits',/,
     *' FPD origin (in mm) at (',f9.3,',',f9.3,')',/,
     *' and FPD x-axis gradient = ',f9.6,' in radiator frame.',/,
C     *' in radiator frame.',/,
     *' Main beam (TOTAL) energy = ',f9.4,' MeV,',
     *' mom fr NMR in opt F = ',f9.4,' MeV/c.')
      WRITE(2,1016)BNMR,BCORR,B1,iopt
1016  FORMAT(' NMR reading = ',f9.7, ' Interpolated BCORR=',f11.8,/,
     *'   (Equivalent uniform field = ',f9.7,' Tesla) IOPT = ',i1,/)
10161 FORMAT(' As field was within 1 percent of 1.443 T, used phemomenological ',/,
     * 'correction derived for 1204 MeV beam',/)
10162 FORMAT(' As field was within 1 percent of 1.834 T, used phemomenological ',/,
     * 'correction derived for 1508 MeV beam',/)
10163 FORMAT(' As field was within 1 percent of 1.89563 T, used phemomenological ',/,
     * 'correction derived for 1557 MeV beam',/)
C
10164 FORMAT(' As field was in 1.00-1.07T range, used phemomenological ',/,
     * 'EUF scaled correction, DM rept 2012/1 for 855 MeV MeV beam',/)
C
C     end jcm mod 8/12/09
C
C     OK read in coordinate positions
C     loop through and determine energy by minimizing distance
C     between trajectory and point
C     First set accuracies
C
1017      FORMAT(1H1)
1018      FORMAT(2(1X,2(A),/),2X,78('-'))
C
C
C     JCM mod 16/5/06
C     STUFF TO GET format without absoft
C
C1005      FORMAT(5X,3G)
C1021      FORMAT(5X,5G)
C1022      FORMAT(5X,2G)
C1002      FORMAT(
C     *    ' Failed to find minimum (IFAIL=',I1,') for point ',I4,/,
C     *    ' Coordinates (x,y) (',F8.2,',',F8.2,')',/,
C     *    ' No of iterations =',I4,/,
C     *    ' Current estimate of:- calibration momentum  =',G,/,
C     *    '                       x position shift   =',G,/,
C     *    ' Current limits on electron momentum',/,
C     *    '                              min ',G,/,
C     *    '                              max ',G,//,
C     *    ' Do you wish to retry (*R) or
C     *accept current estimate (A) [A1]')
C
1005      FORMAT(1X,F4.0,3X,F8.3,3X,F6.3)
C1021      FORMAT(5X,5G)
1022      FORMAT(5X,2F8.4)
1002      FORMAT(
     *    ' Failed to find minimum (IFAIL=',I1,') for point ',I4,/,
     *    ' Coordinates (x,y) (',F8.2,',',F8.2,')',/,
     *    ' No of iterations =',I4,/,
     *    ' Current estimate of:- calibration momentum  =',F8.4,/,
     *    '                       x position shift   =',F8.4,/,
     *    ' Current limits on electron momentum',/,
     *    '                              min ',F8.4,/,
     *    '                              max ',F8.4,//,
     *    ' Do you wish to retry (*R) or
     *accept current estimate (A) [A1]')
C
      DO I=1,NP
C
	XSCINT=X(I)		! Find electron trajectory through
        YSCINT=Y(I)		! scintillator centre.
        RELACC=0.00001
        ABSACC=0.001
CC        RELACC=0.0D0  ! what is this... change back to above jcm 18/8/00
CC        ABSACC=0.0D0
        MAXCAL=10000
        PMIN=0.06*PBEAM
        PMAX=PBEAM
        IFAIL=1
C     Ilya's accuracy stuff  change JCM 18/8/00
        relacc=relacc*pmin+absacc
C     JCM's guess at better values after reading the DMIN blurb 18/8/00
C     taken out 18/1/00
C        relacc=0.0000001
C        absacc=0.000001
C     end jcm mods
C200     CALL E04ABF(RAYTRAK,RELACC,ABSACC,PMIN,PMAX,MAXCAL,
C     *  PCAL,DELX,IFAIL)
200     CALL DMINFC(RAYTRAK,PMIN,PMAX,RELACC,ABSACC, ! Ilya's function
     &    PCAL,DELX,llm)                            ! from CERNlib
C                                                    ! used instead of E04ABF
C        IF(IFAIL.EQ.2) THEN
        IF( .NOT.llm) THEN   ! Ilya's check fit status
C         end of jcm mod
C
          WRITE(*,1002)IFAIL,I,X(I),Y(I),MAXCAL,PCAL,DELX,PMIN,PMAX
          READ(*,1001)ANS
          ANS=UPPER_CASE(ANS)
          IF(ANS.NE.'A') THEN
            MAXCAL=10000
            IFAIL=-1
            GOTO 200
          ENDIF
        ENDIF
        P(I)=PCAL			! Save in arrays the momentum and
	ANG(I)=PHI+90.0			! trajectory angle wrt x-axis.
C
C    *	Find electron trajectory through NEAR edge of scintillator.
C
	XSCINT=X(I)+0.5*width(i)*cosd(theta_fpd(i)+atand(grad))
        YSCINT=Y(I)+0.5*width(i)*sind(theta_fpd(i)+atand(grad))
        RELACC=0.00001
        ABSACC=0.001
CC        RELACC=0.0D0  what the ... did these do? jcm 18/8.00
CC        ABSACC=0.0D0    resurrect the values above
        MAXCAL=10000
        PMIN=0.06*PBEAM
        PMAX=PBEAM
        IFAIL=1
        relacc=relacc*pmin+absacc  ! Ilya's accuracy
C
C        relacc=0.0000001 ! jcm's better values from the blurb on dmin
C        absacc=0.000001 ! NOT put in 18/8/00
C
C201     CALL E04ABF(RAYTRAK,RELACC,ABSACC,PMIN,PMAX,MAXCAL,
C     *  PCAL,DELX,IFAIL)
201     CALL DMINFC(RAYTRAK,PMIN,PMAX,RELACC,ABSACC, ! Ilya's function
     &  PCAL,DELX,llm)                               ! from CERNlib
C                                                    ! insetead of E04ABF
C        IF(IFAIL.EQ.2) THEN
        IF( .NOT.llm) THEN
          WRITE(*,1002)IFAIL,I,X(I),Y(I),MAXCAL,PCAL,DELX,PMIN,PMAX
          READ(*,1001)ANS
          ANS=UPPER_CASE(ANS)
          IF(ANS.NE.'A') THEN
            MAXCAL=10000
            IFAIL=-1
            GOTO 201
          ENDIF
        ENDIF
        P1(I)=PCAL			! Save near-edge momentum in array P1,
	ANG1(I)=PHI+90.0		! trajectory angle wrt x-axis in ANG1.
	if(scint_plot)write(7,1022)xscint,yscint
	if(traj_plot)then
	 write(8,1022)xscint-arm*cosd(ang1(i)),yscint-arm*sind(ang1(i))
	 write(8,1022)xscint+arm*cosd(ang1(i)),yscint+arm*sind(ang1(i))
	end if
C
C
C    *	Find electron trajectory through FAR edge of scintillator.
C
	XSCINT=X(I)-0.5*width(i)*cosd(theta_fpd(i)+atand(grad))
        YSCINT=Y(I)-0.5*width(i)*sind(theta_fpd(i)+atand(grad))
        RELACC=0.00001
        ABSACC=0.001
CC        RELACC=0.0D0  ! what? change back to above values jcm 18/8/00
CC        ABSACC=0.0D0
        MAXCAL=10000
        PMIN=0.06*PBEAM
        PMAX=PBEAM
        relacc=relacc*pmin+absacc    ! Ilya's accuracy jcm 18/8/00
C
C        relacc=0.0000001  !  jcm's guess at accuracy from the dmon blurb
C        absacc=0.000001   !  NOT put in 18/8/00
        IFAIL=1
C202     CALL E04ABF(RAYTRAK,RELACC,ABSACC,PMIN,PMAX,MAXCAL,
C     *  PCAL,DELX,IFAIL)
202     CALL DMINFC(RAYTRAK,PMIN,PMAX,RELACC,ABSACC, ! Ilya's function
     &    PCAL,DELX,llm)                            ! from CERNlib
C        IF(IFAIL.EQ.2) THEN
        IF( .NOT.llm) THEN
          WRITE(*,1002)IFAIL,I,X(I),Y(I),MAXCAL,PCAL,DELX,PMIN,PMAX
          READ(*,1001)ANS
          ANS=UPPER_CASE(ANS)
          IF(ANS.NE.'A') THEN
            MAXCAL=10000
            IFAIL=-1
            GOTO 202
          ENDIF
        ENDIF
        P2(I)=PCAL			! Save far-edge momentum in array P2,
	ANG2(I)=PHI+90.0		! trajectory angle wrt x-axis in ANG2.
	if(scint_plot)write(7,1022)xscint,yscint
	if(traj_plot)then
	 write(8,1022)xscint+arm*cosd(ang2(i)),yscint+arm*sind(ang2(i))
	 write(8,1022)xscint-arm*cosd(ang2(i)),yscint-arm*sind(ang2(i))
	end if
C
      ENDDO
C
C     *	Write out calibration data.
C
	j=1
	do i=1,np
C
	 if(j.eq.1) then
	  if(i.ne.1)write(2,1017)
C   JCM fix for iopt =4 or 5 28/6/99
C         write(6,*) ' got to here ok i,j= ',i,j
            if(iopt.lt.4) then
            write(2,1018)((head(k,l,iopt),k=1,2),l=1,2)
            endif
            if(iopt.gt.3) then
            write(2,1018)((head(k,l,3),k=1,2),l=1,2)
            endif
	 endif
C
	 if(i.gt.1) then
	  angi = 0.5*(ang(i)+ang(i-1))
          iecl=mod(i-1,16)
	  if(iecl.eq.0)iecl=16
C
	  if(iopt.eq.1)then
C
C     *	Write out energy and energy bite for each logic output (Goosy) channel
C     *	based on rays through scintillator centres only.
C
	   ee_mean = 0.5 * (energy(p(i)) + energy(p(i-1)))
	   eg_mean = ebeamm-ee_mean
	   dele = energy(p(i))-energy(p(i-1))
	   write(2,2005)i-1,ee_mean,dele,eg_mean,ecl(i-1),
     1               iecl,hv(i-1),sn(i-1),(2*i-1)
C	1               iecl,hv(i-1),sn(i-1),(2*i-1)  !jcm
C
C     * CUPID file of Goosy ch.,E_e,dE_e,E_g,E_angle...etc as required
C
	   if(filc.ne.' ')then
	    write(3,1005)float(i-1),eg_mean,dele
C	    write(3,1021)float(i-1),ee_mean,dele,eg_mean,angi 
	   endif
C
	  else if(iopt.eq.2)then
C
C     *	Write out upper and lower photon energies for regions of double (only)
C     * scintillator overlap (1 goosy channel firing) and triple scintillator
C     * overlap (2 neighbouring goosy channels firing) using true scintillator 
C     * widths and geometries.
C
	   if(i.eq.2)then				! Region corresponding
	    eg1_hi=ebeamm-energy(p2(i))			! to single output
	   else						! channel firing
	    eg1_hi=ebeamm-energy(max(p1(i-2),p2(i)))	! only.
	   end if					!
	   if(i.eq.353)then				!
	    eg1_lo=ebeamm-energy(p1(i-1))		!
	   else						!
	    eg1_lo=ebeamm-energy(min(p2(i+1),p1(i-1)))	!
	   end if					!
	   eg1_mean=0.5*(eg1_hi+eg1_lo)			!
	   dele1=eg1_hi-eg1_lo				!
C
	   if(i.gt.2)then			! Region corresponding to two
	    eg2_hi=ebeamm-energy(p2(i))		! neighbouring channels
	    eg2_lo=ebeamm-energy(p1(i-2))	! firing.
	    dele2=eg2_hi-eg2_lo			!
	    if(dele2.lt.0.0)dele2=0.0		! < Check that triple scint.
	    eg2_mean=0.5*(eg2_hi+eg2_lo)	!   overlap is finite.
	   end if				!
C
	   if(i.gt.2)then
	    if(dele2.eq.0.0)then
             write(2,2009)i-2,i-1,eg2_lo,eg2_hi,eg2_mean,dele2
	    else
             write(2,2007)i-2,i-1,eg2_lo,eg2_hi,eg2_mean,dele2
	    end if
	   end if
	   write(2,2006)i-1,eg1_lo,eg1_hi,eg1_mean,dele1,ecl(i-1),
     1               iecl,hv(i-1),sn(i-1),(2*i-1)
C	1               iecl,hv(i-1),sn(i-1),(2*i-1) !jcm
C
C     * CUPID file of Goosy ch.,E_e,dE_e,E_g,E_angle...etc as required
C
	   if(filc.ne.' ')then
	    if(i.gt.2)write(3,1005)float(i)-1.5,eg2_mean,dele2
	    write(3,1005)float(i-1),eg1_mean,dele1
	   endif
	   j=j+1			! Extra increment on line counter
C
	  else if(iopt.eq.3)then
C
C     *	Write out upper and lower photon energies for COMBINED region
C     * containing both double (1 goosy channel firing) and triple scintillator
C     * overlap (2 neighbouring goosy channels firing) using true scintillator 
C     * widths and geometries. Neighbouring doubles are included with the
C     * preceding (i.e. lower) photon energy bin.
C
	   eg_hi=ebeamm-energy(p2(i))			! Region containing
	   if(i.eq.353)then				! both single and
	    eg_lo=ebeamm-energy(p1(i-1))		! neighbouring double
	   else						! output channels.
	    eg_lo=ebeamm-energy(min(p2(i+1),p1(i-1)))	!
CCCCC end of JCM mod ? in wrong place?
	   end if					!
	   eg_mean=0.5*(eg_hi+eg_lo)			!
	   dele=eg_hi-eg_lo				!
C
	   write(2,2008)i-1,eg_lo,eg_hi,eg_mean,dele,ecl(i-1),
     1               iecl,hv(i-1),sn(i-1),(2*i-1)
C	1               iecl,hv(i-1),sn(i-1),(2*i-1) !jcm
C
C     * CUPID file of Goosy ch.,E_e,dE_e,E_g,E_angle...etc as required
C
	   if(filc.ne.' ')then
	    write(3,1005)float(i-1),eg_mean,dele
	   endif
C
	  else if(iopt.eq.4)then
C
C     *	Write out upper and lower photon energies for regions of double (only)
C     * scintillator overlap (1 goosy channel firing) 
C     * using true scintillator 
C     * widths and geometries.
C
	   if(i.eq.2)then				! Region corresponding
	    eg1_hi=ebeamm-energy(p2(i))			! to single output
	   else						! channel firing
	    eg1_hi=ebeamm-energy(max(p1(i-2),p2(i)))	! only.
	   end if					!
	   if(i.eq.353)then				!
	    eg1_lo=ebeamm-energy(p1(i-1))		!
	   else						!
	    eg1_lo=ebeamm-energy(min(p2(i+1),p1(i-1)))	!
	   end if					!
	   eg1_mean=0.5*(eg1_hi+eg1_lo)			!
	   dele1=eg1_hi-eg1_lo				!
C
C	   if(i.gt.2)then			! Region corresponding to two
C	    eg2_hi=ebeamm-energy(p2(i))		! neighbouring channels
C	    eg2_lo=ebeamm-energy(p1(i-2))	! firing.
C	    dele2=eg2_hi-eg2_lo			!
C	    if(dele2.lt.0.0)dele2=0.0		! < Check that triple scint.
C	    eg2_mean=0.5*(eg2_hi+eg2_lo)	!   overlap is finite.
C	   end if				!
C
C	   if(i.gt.2)then
C	    if(dele2.eq.0.0)then
C             write(2,2009)i-2,i-1,eg2_lo,eg2_hi,eg2_mean,dele2
C	    else
C             write(2,2007)i-2,i-1,eg2_lo,eg2_hi,eg2_mean,dele2
C	    end if
C	   end if
	   write(2,2006)i-1,eg1_lo,eg1_hi,eg1_mean,dele1,ecl(i-1),
     1               iecl,hv(i-1),sn(i-1),(2*i-1)
C	1               iecl,hv(i-1),sn(i-1),(2*i-1) !jcm
C
C     * CUPID file of Goosy ch.,E_e,dE_e,E_g,E_angle...etc as required
C
	   if(filc.ne.' ')then
C	    if(i.gt.2)write(3,1005)float(i)-1.5,eg2_mean,dele2
	    write(3,1005)float(i-1),eg1_mean,dele1
	   endif
C	   j=j+1			! Extra increment on line counter
C
	  else if(iopt.eq.5)then
C
C     *	Write out upper and lower photon energies for regions of double 
C     * scintillator overlap (1 goosy channel alone or with neighbours firing)
C     * using true scintillator 
C     * widths and geometries.
C
C	   if(i.eq.2)then				! Region corresponding
C	    eg1_hi=ebeamm-energy(p2(i))			! to single output
C	   else						! channel firing
C	    eg1_hi=ebeamm-energy(max(p1(i-2),p2(i)))	! only.
C	   end if					!
C	   if(i.eq.353)then				!
C	    eg1_lo=ebeamm-energy(p1(i-1))		!
C	   else						!
C	    eg1_lo=ebeamm-energy(min(p2(i+1),p1(i-1)))	!
C	   end if					!
C	   eg1_mean=0.5*(eg1_hi+eg1_lo)			!
C	   dele1=eg1_hi-eg1_lo				!
C
           if(i.gt.1)then     ! Region corresponding to two
	    eg2_hi=ebeamm-energy(p2(i))		! neighbouring channels
	    eg2_lo=ebeamm-energy(p1(i-1))	! firing.
C
C     jcm mod 8/12/09 to apply phenom corrns within +or-1% of 1.89563,1.834 and 1.4430 T
C
            IF ((BNMR.GE.1.8767).AND.(BNMR.LE.1.9146)) THEN
C     PHENOMENOLOGICAL CORRECTION for 1.89563 Tesla, 1557 MeV 8/2/09 jcm
C     to correct measured deviations by quadratic fn bkII,p69
CC     to correct measured deviations see p 163, linear up to goosy 163
CC     quadratic 163-352, REVISED from p 171, 18/2/08
            pch=i-1.5
C            if(pch.le.163.0) corrn=0.005601388*pch-1.653349
C      if(pch.gt.163.0) corrn=pch*pch*6.528934E-05-pch*8.402245E-03-1.104
C            if(pch.le.205.0) corrn=0.007506531*pch - 1.774604
C      if(pch.gt.205.0) corrn=pch*pch*1.868485E-04-pch*6.894372E-02 +6.0488
            corrn=pch*pch*7.387262E-05-pch*5.614872E-03 -1.777049
C
CC     suppress write of corrn to screen 14/9/08
C      write(*,*) corrn,pch
C     write message to screen 8/12/09
            IF(i.eq.2) THEN 
               write(*,*) 'Using phenomenological correction derived for 1557 MeV beam'
               write(2,10163)
               ENDIF
C
            eg2_hi=eg2_hi+corrn
            eg2_lo=eg2_lo+corrn
            ENDIF
C     end phenomenological corrn jcm 14/2/08
C
            IF ((BNMR.GE.1.81566).AND.(BNMR.LE.1.85234)) THEN
C     PHENOMENOLOGICAL CORRECTION for 1.834 Tesla, 1508 MeV 14/2/08 jcm
C     to correct measured deviations see p 163, linear up to goosy 163
C     quadratic 163-352, REVISED from p 171, 18/2/08
            pch=i-1.5
C            if(pch.le.163.0) corrn=0.005601388*pch-1.653349
C      if(pch.gt.163.0) corrn=pch*pch*6.528934E-05-pch*8.402245E-03-1.104
            if(pch.le.205.0) corrn=0.007506531*pch - 1.774604
      if(pch.gt.205.0) corrn=pch*pch*1.868485E-04-pch*6.894372E-02 +6.0488
C
C     suppress write of corrn to screen 14/9/08
C      write(*,*) corrn,pch
C     jcm mod 8/12/09 write message to screen
      IF(i.eq.2) THEN
         write(*,*) 'Using phenomeological correction derived for 1508 MeV beam' 
         write(2,10162)
         ENDIF
C     end jcm mod 8/12/09
C
            eg2_hi=eg2_hi+corrn
            eg2_lo=eg2_lo+corrn
C     end phenomenological corrn jcm 14/2/08
            ENDIF
C
            IF ((BNMR.GE.1.42857).AND.(BNMR.LE.1.45743)) THEN
C     PHENOMENOLOGICAL CORRECTION for 1.443 Tesla, 1204 MeV 2/2/09 jcm
C     to correct measured deviations see bk II p56,57 with quad fn p61
C     
            pch=i-1.5
C
      corrn = -pch*pch*1.099716E-04 + pch*2.467585E-02 - 0.384447
C
C      write of corrn to screen 2/2/09
C      write(*,*) corrn,pch
C
C     jcm mod 8/12/09 write message to screen and first output file
      IF(i.eq.2) THEN
         write(*,*) 'Using phenomeological correction derived for 1204 MeV beam' 
      write(2,10161)
      ENDIF
C     end jcm mod 8/12/09
C
            eg2_hi=eg2_hi+corrn
            eg2_lo=eg2_lo+corrn
C     end phenomenological corrn jcm 2/2/09
               ENDIF
C     end of jcm mod 8/12/09
C
C
            IF ((BNMR.GE.1.00000).AND.(BNMR.LE.1.07000)) THEN
C     PHENOMENOLOGICAL CORRECTION for 1.00 - 1.070 Tesla, 840 - 883 MeV 28/3/13 jcm
C     from calibns done for 855 MeV mostly by D. Middleton A2 Int Rept 12/1 
C     coeffnts from DM email DM855_fit_details.txt 6/2/12
C
C     pch changed from 1204,1508,1557MeV above as DM fits are vs goosy not
C     the goosy-0.5 space of JCM's ROOT plots. 
            pch=i-1.0
            corrn=0.00000
            if(pch.le.143.0) corrn=0.920235
            if((pch.gt.143.0).and.(pch.le.246.0)) corrn=-pch*0.0155779+3.16701
            if((pch.gt.246.0).and.(pch.le.274.0)) corrn=-pch*0.0439358+10.1328
            if(pch.gt.274.0) corrn=+pch*0.0173007-6.70907
C     Scale with factor EUF/(EUF for 1.02157T, for which DM's corrns apply) 
            corrn=corrn*B1/1.034413579
C
C      write corrn to screen 14/9/08
C     jcm mod 8/12/09 write message to screen
      IF(i.eq.2) THEN
         write(*,*) 'Using PHENOMENOLOGICAL correction derived for 855 MeV beam' 
         write (*,*) ' Goosy ch, Raw correction, Corrn scaled by EUF '
         write(2,10164)
         ENDIF
      write(*,*) pch,1.034413579*corrn/B1,corrn
C     end jcm mod 2/4/13
C
            eg2_hi=eg2_hi+corrn
            eg2_lo=eg2_lo+corrn
C     end phenomenological corrn jcm 14/2/08
            ENDIF
C
C
	    dele2=eg2_hi-eg2_lo			!
	    if(dele2.lt.0.0)dele2=0.0		! < Check that triple scint.
	    eg2_mean=0.5*(eg2_hi+eg2_lo)	!   overlap is finite.
           end if                 !
C
C	   if(i.gt.2)then
C	    if(dele2.eq.0.0)then
C             write(2,2009)i-2,i-1,eg2_lo,eg2_hi,eg2_mean,dele2
C	    else
C             write(2,2007)i-2,i-1,eg2_lo,eg2_hi,eg2_mean,dele2
C	    end if
C          end if
	   write(2,2010)i-1,eg2_lo,eg2_hi,eg2_mean,dele2,ecl(i-1),
     1               iecl,hv(i-1),sn(i-1),(2*i-1)
C	1               iecl,hv(i-1),sn(i-1),(2*i-1) !jcm
C
C     * CUPID file of Goosy ch.,E_e,dE_e,E_g,E_angle...etc as required
C
	   if(filc.ne.' ')then
CC	    if(i.gt.2)write(3,1005)float(i)-1.5,eg2_mean,dele2
	    write(3,1005)float(i-1),eg2_mean,dele2
	   endif
C	   j=j+1			! Extra increment on line counter
C
C
	  end if
C
C
	 endif
	 j=j+1				! Increment line counter
	 if(j.ge.52)j=1			! New page
	enddo
C
	close(unit=1)
	close(unit=2)
	if(filc.ne.' ')close(unit=3)
	close(unit=7)
	stop
	end
C
C
C
C
C      SUBROUTINE RAYTRAK(MOMENTUM,YSHFT)
      real*8 function RAYTRAK(MOMENTUM) ! Ilya's function instead of subr      
C
C     SUBROUTINE RAYTRAK TO CALCULATE THE SHIFT IN X-COORD OF A
C     RAY EMERGING FROM A QD SPECTROMETER FROM A PRE-DEFINED POINT
C     X,Y WHEN THE Y-VALUE OF THE RAY AND THE Y-VALUE OF THE POINT
C     ARE EQUAL
C     CALCULATION TAKEN FROM PROGRAM TAGQD.FOR BY I. ANTHONY
C     8/5/91
C
      IMPLICIT NONE
C
      REAL*8 MOMENTUM,YSHFT
C
      INTEGER IOPT,IFLAG
C
      REAL*4 SLOPE,TA,TB,DSA,TEMP
      REAL*4 S1X
      REAL*4 X1,Y1,XI,YI
      REAL*4 E2,F2
      REAL*4 XT(2),YT(2)
C
C
      REAL*4 M0,KZ
C
C     Next the common parameters for the minimisation calculation
      REAL*4 B1,X10,Y10,SIX,XSCINT,YSCINT,PHI,RHO
      INTEGER IBEND
      COMMON/MPARMS/B1,X10,Y10,SIX,XSCINT,YSCINT,PHI,RHO,IBEND
      REAL*4 EB,E0,S0X,W0,B0,A,LQ,SD,RHOB,R1,R2,PHIB,
     *ALPH1B,ALPH2B,PHI0,ALPH20,R2B
      COMMON/TPARMS/EB,E0,S0X,W0,B0,A,LQ,SD,RHOB,R1,R2,PHIB,
     *ALPH1B,ALPH2B,PHI0,ALPH20,R2B
C
C     ELECTRON REST ENERGY
      REAL*4 RAD,ARG,MOM,ENERGY
      DATA M0/0.511/,KZ/29.979613E-02/
C      REAL*4 RAD,ARG,MOM,ENERGY  ! shift 2 line up jcm 18/8/00
C     DEFINE FUNCTION TO CONVERT FROM DEGREES TO RADIANS
C      RAD(ARG)=ARG*3.1415926/180.0
C     DEFINE FUNCTION TO CALULATE MOMENTUM FROM ENERGY
C      MOM(ARG)=SQRT(ARG**2-0.511**2)
C      ENERGY(ARG)=SQRT(ARG**2+0.511**2)
C
C
C     JCM mod 16/5/06
C     STUFF TO GET TRIG FUNCTIONS without using absoft (from Ken)
C
      REAL*4 COSD,SIND,ACOSD,ASIND,ATAND,THETAD,XD
      COSD(THETAD)=COS(THETAD*3.1415927/180.)
      SIND(THETAD)=SIN(THETAD*3.1415927/180.)
      ACOSD(XD)=180.*ACOS(XD)/3.1415927
      ASIND(XD)=180.*ASIN(XD)/3.1415927
      ATAND(XD)=180.*ATAN(XD)/3.1415927
C
C     end JCM mod
C
C
C
C******************************************************************
C
C     THIS SECTION CALCULATES TRAJECTORIES FOR NON- CENTRAL RAYS
C
C******************************************************************
C
C
C     RHO IS RADIUS OF BEND IN MAGNET
      RHO= MOMENTUM/(KZ*B1)
C
C
C     CHECK WHETHER EXIT FACE IS STRAIGHT OR CURVED
C
      IF(R2.EQ.0.0) THEN
C
C       EXIT FACE STRAIGHT
C
        IOPT=2
C
C       RAD2 IS SLOPE OF EXIT FACE(DEGREES TO HORIZONTAL)
C
        SLOPE=90.0*(1+IBEND)-IBEND*(PHI0-ALPH20)
        CALL INTRSEC(IBEND*RHO,SIX,RHO,X10,Y10,SLOPE,XT,YT
     *  ,IOPT,IFLAG)
C
C
      ELSE
C
C       CURVED FACE
C
        IOPT=1
C
C       (E2,F2) ARE COORDS. OF CENTRE OF EXIT CIRCULAR FACE
C
        E2= X10 - IBEND*R2*SIND(PHI0-ALPH20)
        F2= Y10 - R2*COSD(PHI0-ALPH20)
C
C     JCM MOD type out E2,F2 for tests 9/2/07
C
C       write(*,10152) E2,F2
10152 FORMAT(' (E2,F2) = ', 2F11.5)
C     END JCM MOD
C
C
C       DETERMINE INTERSECTION OF RAY WITH EXIT FACE
C
        CALL INTRSEC(E2,F2,R2,IBEND*RHO,SIX,RHO,XT,YT
     X  ,IOPT,IFLAG)
C
      ENDIF
C
C     Check that an intersection has been found with exit face
C     i.e. IFLAG=0
      IF(IFLAG.EQ.0) THEN
C
C       NOW CALC. DISTANCE BETWEEN EACH SOLUTION PT. AND THE CENTRAL RAY
C       INTERSECTION PT. ((X10,Y10)) TO DETERMINE CORRECT SOLN.
C
        TA=(XT(1)-X10)**2+(YT(1)-Y10)**2
        TB=(XT(2)-X10)**2+(YT(2)-Y10)**2
C       **********************
        IF(TA.GT.TB) THEN
          X1=XT(2)
          Y1=YT(2)
        ELSE
          X1=XT(1)
          Y1=YT(1)
        END IF
C       **********************
C
C       CALC. DEFLECTION ANGLE, PHI.
C
        PHI=ACOSD(1.0-IBEND*X1/RHO)
C
C
C
        XI=XSCINT
        S1X = (XI-X1)/COSD(90.0-IBEND*PHI)
        YI = Y1 + S1X*SIND(90-PHI)
C
C       SHIFT IN Y DIRECTION IS TAKEN AS SQUARE OF PHYSICAL SHIFT
C       TO ENSURE IT VARIES SMOOTHLY ABOUT 0
        YSHFT=(YI-YSCINT)**2
        RAYTRAK=yshft         ! Ilya's set return value
      ENDIF
      RETURN
      END
C
C
C -----------------------------------------------------------------------------
      CHARACTER*(*) FUNCTION UPPER_CASE(STRING)
C
C     This subroutine is used to change any characters in a string which
C     are lower case alphabetic characters to upper case. The string can
C     be of any length.
C
      IMPLICIT NONE
C
      CHARACTER*(*) STRING  ! this is a passed length string
      INTEGER LS,IC,I
C
      LS=LEN (STRING)
      UPPER_CASE=STRING
      DO I=1,LS
        IC=ICHAR(STRING(I:I))
        IF (IC .GE. 97 .AND. IC .LE. 122)UPPER_CASE(I:I)=CHAR(IC-32)
      END DO
      RETURN
      END
C
      SUBROUTINE INTRSEC(A,B,RAD1,C,D,RAD2,XP,YP,IOPT,ISUCCESS)
C
C     SUBROUTINE INTRSEC DETERMINES INTERSECTION POINTS OF:-
C     IOPT=1, TWO CIRCLES CENTRE (A,B) (C,D)
C                         RADIUS  RAD1  RAD2
C     IOPT=2, STRAIGHT LINE AND CIRCLE
C                         LINE   (X-C)/(Y-D)=TAN(RAD2), RAD2 IN DEGREES
C                         CIRCLE CENTRE (A,B), RADIUS RAD1
C
C     CALC. COEFFICIENTS FOR QUADRATIC EQN., WHICH GIVES INTERSECTION
C     CHOOSE CALC. OF P AND Q COEFFS. ACCORDING TO IOPT= 1,2
C
      REAL*4 XP(2),YP(2)
C
C
C
C     JCM mod 16/5/06
C     STUFF TO GET TRIG FUNCTIONS without using absoft (from Ken)
C
C
      REAL*4 COSD,SIND,ACOSD,ASIND,ATAND,THETAD,XD
      COSD(THETAD)=COS(THETAD*3.1415927/180.)
      SIND(THETAD)=SIN(THETAD*3.1415927/180.)
      ACOSD(XD)=180.*ACOS(XD)/3.1415927
      ASIND(XD)=180.*ASIN(XD)/3.1415927
      ATAND(XD)=180.*ATAN(XD)/3.1415927
C
C     end JCM mod
C
C     DEFINE TAN FUNCTION
C
      TAND(THETA)=SIND(THETA)/COSD(THETA)
      ISUCCESS=1
      IF(IOPT.EQ.1) THEN
C
C       INTERSECTION OF TWO CIRCLES
C
        PN= RAD1**2 - RAD2**2 - A**2 -B**2 +C**2 + D**2
        PD= 2.0*(C - A)
        P=PN/PD
        Q=(D-B)/(C-A)
      ELSE IF(IOPT.EQ.2) THEN
C
C       INTERSECTION OF STRAIGHT LINE AND CIRCLE
C
        IF(RAD2.EQ.90.0.OR.RAD2.EQ.270.0) THEN
C        
C         Line vertical
          P=C
          Q=0.0
        ELSE IF(RAD2.EQ.0.0.OR.RAD2.EQ.180.0)THEN
C
C         LINE HORIZONTAL
C
          P=D
          Q=0.0
        ELSE
          P=C - D/TAND(RAD2)
          Q=-1.0/TAND(RAD2)
        ENDIF
      ENDIF
C
C
C     NOW CALC. REMAINING QUADRATIC PARAMETERS
C
      IF(IOPT.EQ.2.AND.(RAD2.EQ.0.0.OR.RAD2.EQ.180.0)) THEN
        ALPHA=1.0
        BETA=-2.0*A
        GAMMA=A**2+(D-B)**2-RAD1**2
      ELSE
        ALPHA = 1.0+Q**2
        BETA  = 2.0*(Q*(P-A)+B)
        GAMMA = (P-A)**2+B**2-RAD1**2
      ENDIF
C
C     CALC SQRT(BETA**2 - 4*ALPHA*GAMMA) IN SOLN. OF QUADRATIC
C
      SQR=BETA**2-4.0*ALPHA*GAMMA
C     CHECK IF SOLNS. REAL
C
      IF(SQR.GE.0.0) THEN
C
C
C       CALC. SOLN.
C
        SQR=SQRT(SQR)
C
C       (XP(1),YP(1)),(XP(2),YP(2)) ARE TWO SOLNS.
C
        YP(1)=(BETA+SQR)/(2.0*ALPHA)
        YP(2)=(BETA-SQR)/(2.0*ALPHA)
        XP(1)=P-Q*YP(1)
        XP(2)=P-Q*YP(2)
C
C       IF LINE IS HORIZONTAL THEN SWAP X,Y OF SOLNS
C
C****************************
        IF((RAD2.EQ.0.0.OR.RAD2.EQ.180.0).AND.IOPT.EQ.2)THEN
          DO  J=1,2
            PT=YP(J)
            YP(J)=XP(J)
            XP(J)=PT
          END DO
        END IF
C****************************
        ISUCCESS=0
      ELSE
C
C       Type out of intersection error removed I.A. 6/3/86
C        TYPE 100,RAD1,RAD2,P,Q,A,B,C,SQR,IOPT
C100     FORMAT( ' INTERSECTION ERROR '/
C     X  ,' RAD1 = ',G,/,' RAD2 = ',G,/,' P,Q = ',2G,/,' A,B,C = ',3G,/,
C     X  ' SQR = ',G,/,' IOPT = ',I1)
      END IF
      RETURN
      END
C
C
