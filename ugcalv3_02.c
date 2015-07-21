/****************************************************************
 *                                                       	*
 *	Ugcalv3_01.c -- Ugcal in C-Programming              	*
 *                                                       	*
 *	Uses the TMinuit Minimizer from ROOT 
 *                                                       	*
 *      - Robert Materi 2015                             	*
 *                                                       	*
 ***************************************************************/
//==================================================================
//Program ugcalv3_01 to calculate upgraded tagger energy calibration
//==================================================================

// C++ libraries used for all programs in ugcalv3_01.c
#include <cstdio>
#include <cmath>
#include <string>
#include <cstdlib>
#include <iostream>

#define PI 3.1415927

// Global variables that are required in multiple functions
// Trajectory global variables
const double EB = 855.511;
const double E0 = 458.311;
const double S0X = 100.0;
const double LQ = 138.0;
const double SD = 283.0;
const double Rho_B = 2759.599464;
const double R2 = -8021.0;
const double Phi_0 = 74.78526;
const double Alpha_20 = -51.16861;

// Minimizer global variables
const int Bend = -1;

/********************************************************************
*
* These are a collection of small functions which are here to help
* with calculations above to make the computer code easier to read
* Robert 06/10/15
*
********************************************************************/

// Calculates the Momentum and returns the result as a double
double MomentumCal(double arg)
{
double result;
double M_0 = 0.511;
result = sqrt(arg * arg - M_0 * M_0);
return result;
}

/*******************************************************************/

// Calculates the Energy and returns the result as a double
double EnergyCal(double arg)
{
double result;
double M_0 = 0.511;
result = sqrt(arg * arg + M_0 * M_0);
return result;
}

/*******************************************************************/

// These are sine, cosine, tangent, arcsine, arccosine and
// arctangent functions which instead of using radians like the
// "cmath" library does, uses degrees

double cos_deg(double angle)
{
double result;
angle *= PI/180.0;
result = cos(angle);
return result;
}

double sin_deg(double angle)
{
double result;
angle *= PI/180.0;
result = sin(angle);
return result;
}

double tan_deg(double angle)
{
double result;
result = sin_deg(angle)/cos_deg(angle);
return result;
}

double acos_deg(double arg)
{
double angle;
angle = acos(arg);
angle *= 180.0/PI;
return angle;
}

double asin_deg(double arg)
{
double angle;
angle = asin(arg);
angle *= 180.0/PI;
return angle;
}

double atan_deg(double arg)
{
double angle;
angle = atan(arg);
angle *= 180.0/PI;
return angle;
}

/*******************************************************************/

// A struct for the output of the intersection function to place
// data in so that it can return its four entries
struct Intersection_Point
{
	double XP[2];
	double YP[2];
};

/********************************************************************
*
* The Intersection function determines the intersection points of :
* Option = 1, Two Circles centre (A,B) (C,D)
*                       Radius  Rad1  Rad2
* Option = 2, Straight Line and Circle
*               Line (X-C)/(Y-D) = tan(Rad2), Rad2 in Degrees
*               Circle Centre (A, B), Radius Rad1
*
* Calculates coefficients for quadratic equation which gives
* intersection. Choose calculation of P & Q coefficients according
* to Option = 1,2
*
********************************************************************/

Intersection_Point Intersection(double A, double B, double Rad1, double C, double D, double Rad2, int Option, int Success)
{

// Allows access to our struct so that the Intersection function
// can assign values to its components for future use
Intersection_Point Result;

Success = 1;
double PN, PD, P, Q, PT;
double XP[2], YP[2];

// Quadratic equation: Alpha*x^2 + Beta*x + Gamma = 0
double Alpha, Beta, Gamma, Discriminant;

if (Option == 1)
{
        // Intersection of two circles
        PN = Rad1 * Rad1 - Rad2 * Rad2 - A * A - B * B + C * C
                + D * D;
        PD = 2.0 * (C - A);
        P = PN / PD;
        Q = (D - B) / (C - A);
}
else if(Option == 2)
{
        // Intersection of straight line and circle
        if (Rad2 == 90.0 || Rad2 == 270.0)
        {
                // Line Vertical
                P = C;
                Q = 0.0;
        }
        else if (Rad2 == 0.0 || Rad2 == 180.0)
        {
                // Line Horizontal
                P = D;
                Q = 0.0;
        }
        else
        {
                P = C - D / tan_deg(Rad2);
                Q = -1.0 / tan_deg(Rad2);
        }
}

// Now Calculate remaining quadratic parameters
if (Option == 2 && (Rad2 == 0.0 || Rad2 == 180.0))
{
        Alpha = 1.0;
        Beta = -2.0 * A;
        Gamma = A * A + (D - B) * (D - B) - Rad1 * Rad1;
}
else
{
        Alpha = 1.0 + Q * Q;
        Beta = 2.0 * (Q * (P - A) + B);
        Gamma = (P - A) * (P - A) + B * B - Rad1 * Rad1;
}

// Calculate Square Root of (Beta^2 - 4*Alpha*Gamma) in solution of
// quadratic
Discriminant = Beta * Beta - 4.0 * Alpha * Gamma;

// Check if real solutions
if (Discriminant >= 0.0)
{
        // Calculate Solution
        Discriminant = sqrt(Discriminant);
        // (XP(1), YP(1)), (XP(2), YP(2)) are two solutions
        YP[0] = (Beta + Discriminant) / (2.0 * Alpha);
        YP[1] = (Beta - Discriminant) / (2.0 * Alpha);
        XP[0] = P - Q * YP[0];
        XP[1] = P - Q * YP[1];

        // If line is horizontal, swap X, Y of solutions
        if((Rad2 == 0.0 || Rad2 == 180.0) && Option == 2)
        {
                for (int J = 0; J < 2; J++)
                {
                        PT = YP[J];
                        YP[J] = XP[J];
                        XP[J] = PT;
                }
        }
        Success = 0;
}

// These are the intersection points begin stored into our struct
// so that they can be used in our RayTracking function
Result.XP[0] = XP[0];
Result.XP[1] = XP[1];
Result.YP[0] = YP[0];
Result.YP[1] = YP[1];
return Result;
}

/*******************************************************************/

// This is a struct that will be used for holding common variables 
// in the main ugcalv3 function and the RayTrak function 

struct RayTrak_Vars
{
double B1;
double X_10;
double Y_10;
double SIX;
double XScint;
double YScint;
int Bend;
};

/********************************************************************
*
* Raytrak calculates the shift in x-coordinates of a ray emerging
* from a QD Spectrometer from a per-defined point x,y when the
* y-value of the point are equal. Calculation taken from program
* 'TAGQD.FOR' by I.Anthony 8/5/91
*
********************************************************************/

double RayTrak(double Momentum)
{

// Variables
int Option, Flag;

double Slope, TA, TB, DSA, Temp;
double S1X;
double X1, Y1, XI, YI;

// Coordinates of the center of the exit circular face
double E2, F2;

double M_0 = 0.511, KZ = 0.29979613;

// common parameters for the minimisation calculation
double B1, X10, Y10, SIX, XScint, YScint, Phi, Rho;
int Bend;

// Our result, what this function will output
double YShift;

// So we have access to our struct
Intersection_Point Point;

// This Section Calculates Trajectories for Non-Central Rays
//
// Rho is the radius of bend in magnet
Rho = Momentum/(KZ * B1);

// Check whether exit face is straight or curved
if (R2 == 0.0)
{
        // Exit Face Straight
        Option = 2;

        // Rad2 is slope of exit face (degrees to horizontal)
        Slope = 90.0 * (1 + Bend) - Bend * (Phi_0 - Alpha_20);
	Point = Intersection(Bend*Rho, SIX, Rho, X10, Y10, Slope, Option, Flag);
}
else
{
	// Curved Face
	Option = 1;

	// (E2,F2) are coordinates of centre of exit circular face
	E2 = X10 - Bend * R2 * sin_deg(Phi_0 - Alpha_20);
	F2 = Y10 - R2 * cos_deg(Phi_0 - Alpha_20);

	Point =  Intersection(E2, F2, R2, Bend*Rho, SIX, Rho, Option, Flag);
}
if (Flag == 0)
{
        // Calculate distance between each solution point and the
        // central ray intersection point (X10,Y10) to determine
        // the correct solution
        TA = (Result.XP[0]-X10) * (Result.XP[0]-X10) + (Result.YP[0]-Y10) * (Result.YP[0]-(Y10));
        TB = (Result.XP[1]-X10) * (Result.XP[1]-X10) + (Result.YP[1]-Y10) * (Result.YP[1]-(Y10));
        if(TA > TB)
        {
                X1 = Result.XP[1];
                Y1 = Result.YP[1];
        }
        else
        {
                X1 = Result.XP[0];
                Y1 = Result.YP[0];
        }

// Calculates deflection angle Phi
Phi = acos_deg(1.0 - Bend * X1 / Rho);

XI = XScint;
S1X = (XI - X1)/cos_deg(90.0 - Bend * Phi);
YI = Y1 + S1X * sin_deg(90.0 - Phi);

// Shift in y-direction is taken as square of physical shift to 
// ensure it varies smoothly about 0
YShift = (YI - YScint) * (YI - YScint);
return YShift;
}
}

//===================================================================
//
// ugcalv3_01() is the ugcal function for C
//
//===================================================================

void ugcalv3_01()
{

// Variables
const int NDim = 353;
const int NP = 353; 

double X[NDim], Y[NDim], P[NDim], Angle[NDim];
double P1[NDim], P2[NDim];
double X_0, Y_0, Gradient, Angle_deg;
double X_0_Def = -707.134, Y_0_Def = 616.522, Gradient_Def = -0.825231;

// Energies of the photons
double PhoEner1_high, PhoEner1_low, PhoEner1_mean;
double PhoEner2_high, PhoEner2_low, PhoEner2_mean;
double PhoEner_high, PhoEner_low, PhoEner_mean;

// Energy of the electron
double ElecEner_mean;

// Changes in appropriate photon energy
double Delta_E, Delta_E1, Delta_E2;

double Angle1[NDim], Angle2[NDim], pch, corrn;

// Magnetic field corrections
double bcorrm = -0.010876, bcorrc = 1.021246;

//
double Scaling_Factor, AngleI, PBeam; 

// loop index variables
int I, J, K, L;

int ECL;
int Option = 5;

// File output name can be up to 31 characters long (this includes
// the .txt at the end. This is due to the first 31 characters 
// needing to be unique (can't have duplicate names)
char File_Output[31];

// Used when an error occurs and the program needs input from the 
// operator to determine how to proceed. Will store either an 'A'
// for accepting the circumstances or a 'R' which will retry the
// calculation
char ANS;

// A form of labelling for the scintillators so one knows which
// area the scintillator is located (One can find the letters
// corresponding to which scintillator on the back of the 
// scintillators in the A2 hall)
char ecl[352] = {
    'A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A',
    'B','B','B','B','B','B','B','B','B','B','B','B','B','B','B','B',
    'C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
    'D','D','D','D','D','D','D','D','D','D','D','D','D','D','D','D',
    'E','E','E','E','E','E','E','E','E','E','E','E','E','E','E','E',
    'F','F','F','F','F','F','F','F','F','F','F','F','F','F','F','F',
    'G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G',
    'H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H',
    'I','I','I','I','I','I','I','I','I','I','I','I','I','I','I','I',
    'J','J','J','J','J','J','J','J','J','J','J','J','J','J','J','J',
    'K','K','K','K','K','K','K','K','K','K','K','K','K','K','K','K',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O',
    'P','P','P','P','P','P','P','P','P','P','P','P','P','P','P','P',
    'Q','Q','Q','Q','Q','Q','Q','Q','Q','Q','Q','Q','Q','Q','Q','Q',
    'R','R','R','R','R','R','R','R','R','R','R','R','R','R','R','R',
    'S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','S',
    'T','T','T','T','T','T','T','T','T','T','T','T','T','T','T','T',
    'U','U','U','U','U','U','U','U','U','U','U','U','U','U','U','U',
    'V','V','V','V','V','V','V','V','V','V','V','V','V','V','V','V'};

char hv[352] = {
    'A','F','A','F','A','F','A','F','A','F','A','F','A','F','A','F',
    'A','F','A','F','A','F','A','F','A','A','A','A','A','A','A','A',
    'A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A',
    'A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A',
    'A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A',
    'B','B','B','B','B','B','B','B','B','B','B','B','B','B','B','B',
    'B','B','B','B','B','B','B','B','B','B','B','B','B','B','B','B',
    'B','B','B','B','B','B','B','B','B','B','B','B','B','B','B','B',
    'B','B','B','B','B','B','B','B','B','B','B','B','B','B','B','B',
    'B','B','B','B','C','C','C','C','C','C','C','C','C','C','C','C',
    'C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
    'C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
    'C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
    'C','C','C','C','C','C','C','C','D','D','D','D','D','D','D','D',
    'D','D','D','D','D','D','D','D','D','D','D','D','D','D','D','D',
    'D','D','D','D','D','D','D','D','D','D','D','D','D','D','D','D',
    'D','D','D','D','D','D','D','D','D','D','D','D','D','D','D','D',
    'D','D','D','D','D','D','D','D','D','D','D','D','E','E','E','E',
    'E','E','E','E','E','E','E','E','E','E','E','E','E','E','E','E',
    'E','E','E','E','E','E','E','E','E','E','E','E','E','E','E','E',
    'E','E','E','E','E','E','E','E','E','E','E','E','E','E','E','E',
    'E','E','E','E','E','E','E','E','E','E','E','E','E','E','E','E'};

char sn[352] = {
    'L','U','L','U','L','U','L','U','L','U','L','U','L','U','L','U',
    'L','U','L','U','L','U','L','U','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L',
    'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','L'};

//==================================================================
// Variables for NAG routine
//
// The relative accuracy that the minimilisation calculation should
// have 
double Relative_Accuracy = 0.00001;

// The absolute accuracy that the minimilisation calculation can be
double Absolute_Accuracy = 0.001;

// The momentum calculation variables
double PMin, PMax, PCal, Delta_X;

// The maximum number of minimilisation calls that will be called
// (Often this number should not be exceeded)
int Max_Call;

int Fail = 1;
bool llm;

double M_0 = 0.511, KZ = 0.29979613;

//==================================================================
// Parameters for the trajectory calculation
double PB, P0, Rho_0;

// BCORR is the Magnetic field correction constant
double BCORR = 1.00032532;

// The magnetic field reading from the NMR
double BNMR;

// Parameters for the minimisation calculation
double B1, X_10, Y_10, SIX, XScint, YScint, Phi, Rho;
int Bend = -1;
RayTrak_Vars raytrakValues;
raytrakValues.Bend = Bend;

double Momentum;
 
// Next the coord locations of the scintillators. these values are
// those measured in the "FPD" frame of reference (x-axis along the
// rear straight edge of the FPD) and are fixed forever wrt to this 
// frame. The scintillator positions in the "RADIATOR" frame (ie
// wrt the dipole, radiator etc.) depend upon the FPD orientation
// within that frame - in this case defined by the gradient of the
// FPD x-axis (GRAD) and the position of the FPD origin (x0, y0) in
// the radiator frame.
double X_FPD[] = {
       247.9970, 236.3300, 224.6630, 212.9960, 201.3290, 189.6600,
       177.9440, 166.1560, 154.2950, 142.3590, 130.3460, 118.2530,
       106.0800, 93.82600, 81.49000, 69.07200, 56.57300, 43.99600,
       31.34200, 18.61700, 5.826000,-7.024000,-19.92500,-32.86600,
      -45.83700,-58.82800,-71.82700,-84.82500,-97.81300,-110.7860,
      -123.7400,-136.6760,-149.5960,-162.5060,-175.4110,-188.3190,
      -201.2370,-214.1700,-227.1140,-240.0650,-253.0210,-265.9820,
      -278.9480,-291.9190,-304.8930,-317.8710,-330.8520,-343.8360,
      -356.8230,-369.8130,-382.8040,-395.7970,-408.7910,-421.7870,
      -434.7840,-447.7820,-460.7810,-473.7800,-486.7790,-499.7790,
      -512.7790,-525.7790,-538.7790,-551.7790,-564.7780,-577.7770,
      -590.7750,-603.7730,-616.7700,-629.7660,-642.7610,-655.7550,
      -668.7490,-681.7420,-694.7330,-707.7240,-720.7130,-733.7010,
      -746.6880,-759.6740,-772.6590,-785.6420,-798.6250,-811.6060,
      -824.5850,-837.5640,-850.5410,-863.5180,-876.4930,-889.4670,
      -902.4390,-915.4110,-928.3810,-941.3500,-954.3170,-967.2840,
      -980.2500,-993.2140,-1006.178,-1019.140,-1032.101,-1045.061,
      -1058.021,-1070.979,-1083.936,-1096.893,-1109.848,-1122.803,
      -1135.757,-1148.710,-1161.663,-1174.614,-1187.565,-1200.515,
      -1213.465,-1226.413,-1239.361,-1252.309,-1265.256,-1278.203,
      -1291.149,-1304.095,-1317.040,-1329.984,-1342.929,-1355.873,
      -1368.817,-1381.761,-1394.704,-1407.647,-1420.590,-1433.532,
      -1446.475,-1459.417,-1472.360,-1485.302,-1498.245,-1511.187,
      -1524.129,-1537.071,-1550.014,-1562.957,-1575.900,-1588.842,
      -1601.785,-1614.729,-1627.672,-1640.616,-1653.560,-1666.504,
      -1679.448,-1692.393,-1705.338,-1718.284,-1731.230,-1744.177,
      -1757.124,-1770.071,-1783.018,-1795.967,-1808.916,-1821.865,
      -1834.815,-1847.765,-1860.716,-1873.667,-1886.620,-1899.572,
      -1912.525,-1925.480,-1938.434,-1951.389,-1964.345,-1977.301,
      -1990.258,-2003.216,-2016.175,-2029.134,-2042.094,-2055.055,
      -2068.016,-2080.978,-2093.940,-2106.904,-2119.868,-2132.833,
      -2145.799,-2158.765,-2171.732,-2184.701,-2197.669,-2210.639,
      -2223.609,-2236.580,-2249.551,-2262.523,-2275.498,-2288.472,
      -2301.447,-2314.422,-2327.399,-2340.376,-2353.354,-2366.333,
      -2379.312,-2392.292,-2405.273,-2418.255,-2431.237,-2444.220,
      -2457.204,-2470.188,-2483.173,-2496.159,-2509.145,-2522.132,
      -2535.120,-2548.107,-2561.096,-2574.086,-2587.076,-2600.066,
      -2613.057,-2626.049,-2639.041,-2652.033,-2665.027,-2678.020,
      -2691.014,-2704.009,-2717.003,-2729.999,-2742.994,-2755.991,
      -2768.987,-2781.984,-2794.982,-2807.979,-2820.977,-2833.975,
      -2846.973,-2859.971,-2872.970,-2885.969,-2898.968,-2911.967,
      -2924.966,-2937.966,-2950.966,-2963.966,-2976.965,-2989.965,
      -3002.965,-3015.966,-3028.966,-3041.966,-3054.966,-3067.966,
      -3080.965,-3093.965,-3106.964,-3119.963,-3132.962,-3145.960,
      -3158.958,-3171.956,-3184.954,-3197.952,-3210.949,-3223.945,
      -3236.942,-3249.937,-3262.932,-3275.927,-3288.921,-3301.914,
      -3314.906,-3327.898,-3340.889,-3353.879,-3366.869,-3379.858,
      -3392.846,-3405.834,-3418.821,-3431.807,-3444.791,-3457.776,
      -3470.759,-3483.740,-3496.721,-3509.701,-3522.680,-3535.657,
      -3548.634,-3561.609,-3574.583,-3587.555,-3600.527,-3613.497,
      -3626.466,-3639.433,-3652.399,-3665.364,-3678.326,-3691.288,
      -3704.248,-3717.207,-3730.162,-3743.118,-3756.071,-3769.023,
      -3781.973,-3794.920,-3807.866,-3820.811,-3833.753,-3846.693,
      -3859.632,-3872.568,-3885.503,-3898.434,-3911.365,-3924.292,
      -3937.218,-3950.141,-3963.063,-3975.981,-3988.898,-4001.812,
      -4014.724,-4027.632,-4040.539,-4053.443,-4066.345,-4079.244,
      -4092.140,-4105.034,-4117.925,-4130.813,-4143.698,-4156.581,
      -4169.460,-4182.336,-4195.210,-4208.081,-4220.948,-4233.813,
      -4246.675,-4259.534,-4272.390,-4285.241,-4298.091};

double Y_FPD[] = {
      6.714000,12.44800,18.18200,23.91700,29.65100,35.38200,
      41.01500,46.49700,51.81900,56.96900,61.93600,66.70800,
      71.27100,75.61000,79.71000,83.55700,87.13100,90.41700,
      93.39700,96.05500,98.37500,100.3420,101.9450,103.1770,
      104.0330,104.5170,104.6370,104.4120,103.8640,103.0270,
      101.9400,100.6510,99.21200,97.67900,96.11200,94.56900,
      93.11200,91.79700,90.59300,89.45900,88.38900,87.38700,
      86.44800,85.57400,84.75800,84.00200,83.30400,82.66300,
      82.07800,81.54500,81.06700,80.64100,80.26300,79.93600,
      79.65600,79.42200,79.23500,79.09200,78.99300,78.93700,
      78.92200,78.94700,79.01100,79.11600,79.25600,79.43400,
      79.64800,79.89700,80.17900,80.49400,80.84300,81.22300,
      81.63400,82.07400,82.54400,83.04200,83.56800,84.12100,
      84.70000,85.30500,85.93600,86.59000,87.26700,87.96800,
      88.69100,89.43600,90.20200,90.98800,91.79500,92.62000,
      93.46400,94.32700,95.20700,96.10400,97.01800,97.94800,
      98.89300,99.85400,100.8280,101.8180,102.8200,103.8350,
      104.8630,105.9040,106.9560,108.0190,109.0930,110.1770,
      111.2720,112.3750,113.4890,114.6100,115.7400,116.8780,
      118.0230,119.1760,120.3360,121.5010,122.6740,123.8510,
      125.0340,126.2210,127.4130,128.6100,129.8110,131.0150,
      132.2240,133.4340,134.6470,135.8620,137.0800,138.2990,
      139.5200,140.7420,141.9650,143.1880,144.4110,145.6340,
      146.8590,148.0800,149.3020,150.5230,151.7420,152.9600,
      154.1750,155.3900,156.6010,157.8090,159.0160,160.2180,
      161.4180,162.6140,163.8050,164.9930,166.1760,167.3550,
      168.5300,169.7000,170.8630,172.0230,173.1760,174.3240,
      175.4650,176.6000,177.7290,178.8520,179.9680,181.0760,
      182.1780,183.2730,184.3590,185.4380,186.5090,187.5710,
      188.6270,189.6720,190.7100,191.7390,192.7580,193.7700,
      194.7710,195.7630,196.7450,197.7180,198.6810,199.6330,
      200.5760,201.5080,202.4300,203.3400,204.2400,205.1300,
      206.0070,206.8740,207.7290,208.5730,209.4050,210.2250,
      211.0340,211.8300,212.6130,213.3850,214.1440,214.8910,
      215.6250,216.3460,217.0530,217.7480,218.4290,219.0980,
      219.7530,220.3940,221.0220,221.6360,222.2360,222.8220,
      223.3940,223.9510,224.4940,225.0230,225.5370,226.0370,
      226.5220,226.9920,227.4470,227.8860,228.3110,228.7210,
      229.1150,229.4940,229.8570,230.2040,230.5370,230.8530,
      231.1530,231.4370,231.7050,231.9570,232.1930,232.4130,
      232.6160,232.8020,232.9720,233.1260,233.2620,233.3820,
      233.4850,233.5720,233.6400,233.6920,233.7270,233.7440,
      233.7450,233.7270,233.6930,233.6400,233.5700,233.4820,
      233.3770,233.2540,233.1130,232.9540,232.7770,232.5820,
      232.3690,232.1390,231.8880,231.6210,231.3350,231.0300,
      230.7070,230.3650,230.0050,229.6270,229.2290,228.8120,
      228.3780,227.9240,227.4500,226.9580,226.4480,225.9180,
      225.3690,224.8000,224.2130,223.6060,222.9810,222.3350,
      221.6700,220.9860,220.2820,219.5590,218.8160,218.0540,
      217.2720,216.4700,215.6480,214.8070,213.9460,213.0640,
      212.1640,211.2440,210.3020,209.3410,208.3600,207.3600,
      206.3380,205.2970,204.2350,203.1540,202.0520,200.9290,
      199.7860,198.6230,197.4400,196.2370,195.0120,193.7670,
      192.5020,191.2160,189.9100,188.5830,187.2340,185.8660,
      184.4770,183.0670,181.6370,180.1860,178.7130,177.2210,
      175.7070,174.1730,172.6170,171.0410,169.4430,167.8250,
      166.1860,164.5250,162.8440,161.1420,159.4180,157.6730,
      155.9080,154.1220,152.3140,150.4840,148.6340,146.7620,
      144.8700,142.9560,141.0210,139.0640,137.0870};

// Next the loacting angle (in degrees) of the scintillators. These
// values are those measured in the "FPD" frame of reference (x-axis
// along the rear straight edge of the FPD) and are fixed forever
// wrt this frame.
double Theta_FPD[] = {
       133.6296,133.6290,133.6291,133.6291,133.6293,133.8052,
       134.4877,135.2440,136.0448,136.8943,137.7975,138.7551,
       139.7727,140.8513,141.9976,143.2102,144.4933,145.8470,
       147.2704,148.7554,150.3055,149.3918,148.5176,147.6798,
       146.8821,146.1218,145.3988,144.7126,144.0597,143.4384,
       142.8475,142.2822,141.7398,141.2183,140.7134,140.2221,
       139.7415,139.2690,138.8053,138.3518,137.9087,137.4748,
       137.0500,136.6333,136.2258,135.8265,135.4351,135.0516,
       134.6757,134.3066,133.9453,133.5906,133.2427,132.9011,
       132.5664,132.2371,131.9140,131.5971,131.2853,130.9799,
       130.6790,130.3843,130.0941,129.8094,129.5291,129.2540,
       128.9835,128.7176,128.4562,128.1993,127.9468,127.6982,
       127.4540,127.2136,126.9775,126.7448,126.5162,126.2914,
       126.0700,125.8522,125.6381,125.4272,125.2202,125.0160,
       124.8152,124.6175,124.4230,124.2315,124.0434,123.8579,
       123.6755,123.4961,123.3193,123.1454,122.9741,122.8060,
       122.6400,122.4766,122.3160,122.1579,122.0023,121.8489,
       121.6980,121.5498,121.4036,121.2597,121.1182,120.9790,
       120.8418,120.7068,120.5742,120.4434,120.3148,120.1883,
       120.0636,119.9411,119.8205,119.7018,119.5852,119.4702,
       119.3573,119.2461,119.1367,119.0292,118.9234,118.8194,
       118.7168,118.6164,118.5175,118.4201,118.3244,118.2304,
       118.1379,118.0470,117.9577,117.8699,117.7836,117.6989,
       117.6156,117.5337,117.4533,117.3744,117.2969,117.2207,
       117.1460,117.0725,117.0006,116.9299,116.8606,116.7925,
       116.7257,116.6602,116.5960,116.5330,116.4713,116.4107,
       116.3515,116.2934,116.2364,116.1807,116.1261,116.0727,
       116.0203,115.9692,115.9192,115.8701,115.8223,115.7755,
       115.7298,115.6850,115.6414,115.5989,115.5574,115.5168,
       115.4773,115.4387,115.4012,115.3646,115.3290,115.2944,
       115.2607,115.2279,115.1961,115.1652,115.1352,115.1062,
       115.0780,115.0507,115.0243,114.9987,114.9740,114.9503,
       114.9273,114.9052,114.8838,114.8634,114.8437,114.8249,
       114.8068,114.7896,114.7731,114.7575,114.7426,114.7284,
       114.7151,114.7025,114.6906,114.6794,114.6691,114.6594,
       114.6505,114.6423,114.6348,114.6280,114.6218,114.6164,
       114.6117,114.6077,114.6043,114.6016,114.5996,114.5982,
       114.5975,114.5974,114.5980,114.5992,114.6010,114.6035,
       114.6065,114.6103,114.6146,114.6195,114.6250,114.6311,
       114.6378,114.6451,114.6530,114.6615,114.6705,114.6801,
       114.6903,114.7009,114.7123,114.7241,114.7364,114.7494,
       114.7628,114.7768,114.7913,114.8064,114.8220,114.8380,
       114.8546,114.8717,114.8893,114.9074,114.9260,114.9451,
       114.9647,114.9848,115.0054,115.0265,115.0480,115.0700,
       115.0924,115.1154,115.1387,115.1626,115.1869,115.2116,
       115.2368,115.2625,115.2886,115.3151,115.3421,115.3695,
       115.3973,115.4255,115.4542,115.4833,115.5128,115.5427,
       115.5731,115.6038,115.6349,115.6665,115.6984,115.7308,
       115.7635,115.7966,115.8301,115.8640,115.8983,115.9329,
       115.9680,116.0033,116.0391,116.0753,116.1118,116.1487,
       116.1859,116.2235,116.2614,116.2998,116.3384,116.3774,
       116.4168,116.4565,116.4965,116.5368,116.5776,116.6186,
       116.6600,116.7017,116.7438,116.7861,116.8288,116.8718,
       116.9151,116.9588,117.0027,117.0470,117.0915,117.1364,
       117.1816,117.2271,117.2729,117.3190,117.3654,117.4121,
       117.4591,117.5064,117.5539,117.6018,117.6499,117.6984,
       117.7471,117.7961,117.8454,117.8949,117.9447,117.9948,
       118.0452,118.0958,118.1467,118.1979,118.2493,118.3010,
       118.3530,118.4052,118.4577,118.5105,118.5635};


// Next the actual scintillator widths (in mm) employed in the FPD
// construction. Note how the widths are batched where the nominal
// variation is slow.
double Width[] = {
       32.607769,32.221680,31.852831,31.501040,31.164835,30.803890,
       30.347589,29.898945,29.464128,29.040533,28.626196,28.220430,
       27.821032,27.426767,27.035423,26.647329,26.259811,25.872179,
       25.483376,25.093826,24.701313,24.284260,23.821457,23.317320,
       22.778624,22.215502,21.637844,21.062616,20.501245,19.973793,
       19.492910,19.075768,18.70    ,18.70    ,18.70    ,18.70    ,
       18.70    ,18.70    ,18.70    ,18.21    ,18.21    ,18.21    ,
       18.21    ,18.21    ,18.21    ,18.21    ,18.21    ,18.21    ,
       18.21    ,17.75    ,17.75    ,17.75    ,17.75    ,17.75    ,
       17.75    ,17.75    ,17.75    ,17.38    ,17.38    ,17.38    ,
       17.38    ,17.38    ,17.38    ,17.38    ,17.38    ,17.38    ,
       17.38    ,17.38    ,16.95    ,16.95    ,16.95    ,16.95    ,
       16.95    ,16.95    ,16.95    ,16.95    ,16.95    ,16.60    ,
       16.60    ,16.60    ,16.60    ,16.60    ,16.60    ,16.60    ,
       16.60    ,16.60    ,16.60    ,16.25    ,16.25    ,16.25    ,
       16.25    ,16.25    ,16.25    ,16.25    ,16.25    ,16.25    ,
       15.95    ,15.95    ,15.95    ,15.95    ,15.95    ,15.95    ,
       15.95    ,15.95    ,15.95    ,15.60    ,15.60    ,15.60    ,
       15.60    ,15.60    ,15.60    ,15.60    ,15.60    ,15.60    ,
       15.60    ,15.30    ,15.30    ,15.30    ,15.30    ,15.30    ,
       15.30    ,15.30    ,15.30    ,15.30    ,15.30    ,14.95    ,
       14.95    ,14.95    ,14.95    ,14.95    ,14.95    ,14.95    ,
       14.95    ,14.95    ,14.95    ,14.60    ,14.60    ,14.60    ,
       14.60    ,14.60    ,14.60    ,14.60    ,14.60    ,14.60    ,
       14.60    ,14.60    ,14.25    ,14.25    ,14.25    ,14.25    ,
       14.25    ,14.25    ,14.25    ,14.25    ,14.25    ,14.25    ,
       14.25    ,13.90    ,13.90    ,13.90    ,13.90    ,13.90    ,
       13.90    ,13.90    ,13.90    ,13.90    ,13.90    ,13.90    ,
       13.55    ,13.55    ,13.55    ,13.55    ,13.55    ,13.55    ,
       13.55    ,13.55    ,13.55    ,13.55    ,13.55    ,13.25    ,
       13.25    ,13.25    ,13.25    ,13.25    ,13.25    ,13.25    ,
       13.25    ,13.25    ,13.25    ,12.90    ,12.90    ,12.90    ,
       12.90    ,12.90    ,12.90    ,12.90    ,12.90    ,12.90    ,
       12.90    ,12.90    ,12.90    ,12.60    ,12.60    ,12.60    ,
       12.60    ,12.60    ,12.60    ,12.60    ,12.60    ,12.60    ,
       12.60    ,12.60    ,12.60    ,12.25    ,12.25    ,12.25    ,
       12.25    ,12.25    ,12.25    ,12.25    ,12.25    ,12.25    ,
       12.25    ,12.25    ,12.25    ,11.9     ,11.9     ,11.9     ,
       11.9     ,11.9     ,11.9     ,11.9     ,11.9     ,11.9     ,
       11.9     ,11.9     ,11.9     ,11.60    ,11.60    ,11.60    ,
       11.60    ,11.60    ,11.60    ,11.60    ,11.60    ,11.60    ,
       11.60    ,11.60    ,11.60    ,11.30    ,11.30    ,11.30    ,
       11.30    ,11.30    ,11.30    ,11.30    ,11.30    ,11.30    ,
       11.30    ,11.30    ,11.30    ,11.00    ,11.00    ,11.00    ,
       11.00    ,11.00    ,11.00    ,11.00    ,11.00    ,11.00    ,
       11.00    ,11.00    ,11.00    ,10.75    ,10.75    ,10.75    ,
       10.75    ,10.75    ,10.75    ,10.75    ,10.75    ,10.75    ,
       10.75    ,10.75    ,10.75    ,10.45    ,10.45    ,10.45    ,
       10.45    ,10.45    ,10.45    ,10.45    ,10.45    ,10.45    ,
       10.45    ,10.45    ,10.45    ,10.20    ,10.20    ,10.20    ,
       10.20    ,10.20    ,10.20    ,10.20    ,10.20    ,10.20    ,
       10.20    ,10.20    ,10.20    , 9.95    , 9.95    , 9.95    ,
        9.95    , 9.95    , 9.95    , 9.95    , 9.95    , 9.95    ,
        9.95    , 9.95    , 9.95    , 9.95    , 9.70    , 9.70    ,
        9.70    , 9.70    , 9.70    , 9.70    , 9.70    , 9.70    ,
        9.70    , 9.70    , 9.70    , 9.70    , 9.70    , 9.45    ,
        9.45    , 9.45    , 9.45    , 9.45    , 9.45    , 9.45    ,
        9.45    , 9.45    , 9.45    , 9.45    , 9.45    , 9.45    ,
        9.20    , 9.20    , 9.20    , 9.20    , 9.20    };


X_0 = X_0_Def;
Y_0 = Y_0_Def;
Gradient = Gradient_Def;

// Scintillator positions in "RADIATOR" frame of reference
Angle_deg = atan_deg(Gradient);
int index;
for (index = 0; index < NP; index++)
{
	X[index] = X_FPD[index] * cos_deg(Angle_deg)
	- Y_FPD[index] * sin_deg(Angle_deg) + X_0;
	Y[index] = X_FPD[index] * sin_deg(Angle_deg)
	+ Y_FPD[index] * cos_deg(Angle_deg) + Y_0;
}

// Standard ideal magnetic field (Tesla)
double BField;

// Average Energy of the beam (MeV)
double Beam_Energy;

// Screen Output
printf("\t\t Upgraded Tagger Calibration ugcalv2ud\n");
printf("\tJ.C. McGeorge 20/10/10, I. Anthony 7/2/92, G.J.Miller\n\n");
printf("Calculates energy calibration of the upgraded Glasgow Mainz\n");
printf("tagging spectrometer.\n\n");
printf("USES INTERPOLATED FIELD CORRECTION FACTOR, BCORR\n");
printf("UNKNOWN UNCERTAINTY AT FIELDS LESS THAN 0.5 TESLA\n");
printf("PHENOMEOLOGICAL ENERGY CORRECTION for effect of large-scale\n");
printf("field non-uniformity derived from calibration measurements\n");
printf("made April and December 07 BUT ONLY FOR fields within\n");
printf("1 percent of 1.89563, 1.834 and 1.443 Tesla where determined\n\n");
printf("ALSO for fields 1.000 - 1.070T from calibration measurements");
printf("July and Oct 2011: see D. Middleton A2 Internal Report 2012/1\n");
printf("Outputs PHOTON energies & channel widths for\n");
printf("single and neighbouring double channels combined\n");
printf("True scintillator geometry.\n\n");
printf("An output table in 5X3G format can also be produced suitable\n");
printf("for display by a graphics package\n");
printf("Input data required:-\n");
printf("Main beam (TOTAL) energy (MeV)\n");
printf("Tagger magnetic field from NMR (Tesla)\n\n");
printf("Assumes coordinate system origin at the radiator with +ve Y axis\n");
printf("in the direction of the input electron beam. Scintillator locations\n");
printf("are ordered from low to high electron momentum. Channels are then\n");
printf("numbered increasing in this direction from channel 1 (first\n");
printf("disc./coinc. card which has a signal output i.e scint. station\n");
printf("no. 2). All energies are TOTAL.\n");
printf("***************\n\n\n");
printf("Give magnet field NMR reading (Tesla) [G]:\n");
scanf("%lf", &BNMR);

if (BNMR >= 1.89563)
{
	BCORR = -0.014983223 * BNMR + 1.027677707;
}
else if (BNMR >= 1.83400 && BNMR < 1.89563)
{
	BCORR = -0.017041392 * BNMR + 1.031579233;
}
else if (BNMR >= 1.57037 && BNMR < 1.83400)
{
	BCORR = -0.019783283 * BNMR + 1.036607862;
}
else if (BNMR >= 1.44300 && BNMR < 1.57037)
{
	BCORR = -0.024104679 * BNMR + 1.043394052;
}
else if (BNMR >= 1.199669 && BNMR < 1.44300)
{
	BCORR = -0.002335444 * BNMR + 1.011981046;
}
else if (BNMR >= 1.05700 && BNMR < 1.199669) 	// No Negative?
{
	BCORR = 0.0012075924 * BNMR + 1.007730575;
}
else if (BNMR >= 0.7545836 && BNMR < 1.05700)
{
	BCORR = -0.000301227 * BNMR + 1.009325397;
}
else if (BNMR >= 0.5394498 && BNMR < 0.7545836)
{
	BCORR = -0.0025503059 * BNMR + 1.011022515;
}
else
{
	BCORR = -0.000683625 * BNMR + 1.010015534;
}

PB = MomentumCal(EB);
B1 = BCORR * BNMR;           // Ideal magnet field requested
BField = PB / (KZ * RhoB);   // Standard Ideal Field
Scaling_Factor = B1 / BField;// Scaling factor from standard setting
raytrakValues.B1 = B1;
	
PBeam = PB * Scaling_Factor;

printf("Enter TOTAL (include rest mass) ebeam energy from MAMI (MeV):\n");
scanf("%lf", &Beam_Energy);
printf("Using total ebeam energy = %lf\n\n", Beam_Energy);
printf("Using standard FPD location origin (in mm) and gradient:\n");
printf("X_0 = %9.3lf, Y_0 = %9.3lf, Gradient = %9.6lf\n\n", X_0, Y_0, Gradient);
printf("Give filename for printed output [A30]:");
scanf("%s", File_Output);
printf("Using BCORR = %11.8lf\n", BCORR);


/********************************************************************
*
*	FIRST SECTION CALCULATES MAIN BEAM TRAJECTORY
*
********************************************************************/

SIX = S0X + LQ + SD;
raytrakValues.SIX = SIX;

/********************************************************************
*
*	SECOND SECTION CALCULATES CENTRAL TRAJECTORY
*
********************************************************************/

// Central Ray Momentum
Momentum = MomentumCal(E0);
P0 = Scaling_Factor * Momentum;

// Rho_0 is the radius of the central ray in the magnet
Rho_0 = P0/(KZ*B1);

// Calculate intersection points of central ray with magnet face
// exit face intersection (X_10, Y_10)
X_10 = Bend * (Rho_0 - Rho_0 * cos_deg(Phi_0));
raytrakValues.X_10 = X_10;
Y_10 = SIX + Rho_0 * sin_deg(Phi_0);
raytrakValues.Y_10 = Y_10;

printf("(X_10, Y_10) Rho_0 = %11.5lf %11.5lf %11.5lf\n", X_10, Y_10, Rho_0);

// Text File Output
FILE *fp;
fp = fopen(File_Output, "w");
fprintf(fp,"ugcalv3_01.c : Used interpolated field correction factor, BCORR\n");
fprintf(fp,"(Uncertainty is unknown for fields less than 0.5 Tesla)\n");
fprintf(fp,"Channel widths include both neighbouring double hits\n");
fprintf(fp,"FPD origin (in mm) at (%9.3lf, %9.3lf)\n", X_0, Y_0);
fprintf(fp,"and FPD x-axis gradient = %9.6lf in radiator frame.\n", Gradient);
fprintf(fp,"Main beam (TOTAL) energy = %9.4lf MeV,\n", Beam_Energy);
fprintf(fp,"mom fr NMR in opt F = %9.4lf MeV/c.\n", PBeam);
fprintf(fp,"NMR reading = %9.7lf\n", BNMR);
fprintf(fp,"Interpolated BCORR = %11.8lf\n", BCORR);
fprintf(fp,"(Equivalent uniform field = %9.7lf Tesla)\n", B1);
fprintf(fp,"Option = %d\n", Option);

for (I = 0; I < NP; I++)
{
   XScint = X[I]; // Find electron trajectory through
   YScint = Y[I]; // scintillator centre
   Max_Call = 10000;
   PMin = 0.06 * PBeam;
   PMax = PBeam;

   Fail = 1;
   Relative_Accuracy *= PMin + Absolute_Accuracy;

   // Cernlib stuff

   P[I] = PCal; 		// Save in arrays the momentum and
   Angle[I] = Phi + 90.0;	// trajectory angle wrt x-axis

   // Find electron trajectory through NEAR edge of scintillator

   XScint = X[I] + 0.5 * Width[I] * cos((Theta_FPD[I] + atan(Gradient)) * PI / 180);
   YScint = Y[I] + 0.5 * Width[I] * sin((Theta_FPD[I] + atan(Gradient)) * PI / 180);
   Relative_Accuracy = 0.00001;
   Absolute_Accuracy = 0.001;
   Max_Call = 10000;
   PMin = 0.06 * PBeam;
   PMax = PBeam;
   Fail = 1;
   Relative_Accuracy *= PMin + Absolute_Accuracy;

   // cernlib stuff

   P1[I] = PCal;		// Save near-edge momentum in array
   Angle1[I] = Phi + 90.0;	// P1, trajectory angle wrt x-axis in
				// Angle1.
	
   // Find electron trajectory through FAR edge of scintillator

   XScint = X[I] - 0.5 * Width[I] * cos((Theta_FPD[I] + atan(Gradient)) * PI / 180);
   YScint = Y[I] - 0.5 * Width[I] * sin((Theta_FPD[I] + atan(Gradient)) * PI / 180);
   Relative_Accuracy = 0.00001;
   Absolute_Accuracy = 0.001;
   Max_Call = 10000;
   PMin = 0.06 * PBeam;
   PMax = PBeam;
   Relative_Accuracy *= PMin + Absolute_Accuracy;
   Fail = 1;

   // cernlib
	
   P2[I] = PCal;		// Save far-edge momentum in array P2
   Angle2[I] = Phi + 90.0;	// Trajectory angle wrt x-axis in 
				// Angle2
}

// Calibration Data
J = 0;
for (I = 0; I < NP; I++)
{
   if (J == 1)
   {
      if (I != 0)
      {
         if (Option < 4)
	 {
	    fprintf(fp,"\n\n  Goosy    Mean E_e    Bite     Mean E_g      ECL      HV     Station\n");
            fprintf(fp,"  Chan.      (MeV)     (MeV)     (MeV)        out             Number\n");
            fprintf(fp,"**************************************************************************************\n");
	 }
	 else if (Option > 3)
	 {	
            fprintf(fp,"\n\n   Goosy ch.     E_g_Low    E_g_High    E_g_Bite       ECL      HV     Station\n");
            fprintf(fp,"    firing        (MeV)       (MeV)       (MeV)        out             number\n");
            fprintf(fp,"**************************************************************************************\n");
	 }
      }
   }
   if (I > 0)
   {
      AngleI = 0.5 * (Angle[I] + Angle[I - 1]);
      ECL = (I - 1) % 16;
      if (ECL == 0)
      {
         ECL = 16;
      }
      if (Option == 1)
      {
         // Write out energy and energy bite for each logic output
         // (Goosy) channel based on rays through scintillator
	 // centres only
				
	 ElecEner_mean = 0.5 * (EnergyCal(P[I]) + EnergyCal(P[I - 1])); 
         PhoEner_mean =  Beam_Energy - ElecEner_mean;
	 Delta_E = EnergyCal(P[I]) - EnergyCal(P[I - 1]);

	 // Writing information to file
	 fprintf(fp,"%6d %11.3f %11.3f %11.3f %4c %2d %5c %5c %d\n", I - 1, ElecEner_mean, Delta_E, PhoEner_mean, ecl[I-1], ECL, hv[I-1], sn[I-1],2*I-1);
      }
      else if (Option == 2)
      {
	 // Write out upper and lower photon energies for regions of
	 // double (only) scintillator overlap (1 goosy channel
	 // firing) and triple scintillator overlap (2 neighbouring
	 // goosy channels firing) using true scintillator widths
	 // and geometries

	 // Region corresponding to single output channel firing only
         if (I == 1)
         {
            PhoEner1_high = Beam_Energy - EnergyCal(P2[I]);
         }
         else
	 {
            PhoEner1_high = Beam_Energy - EnergyCal((fmax(P1[I - 2], P2[I])));
         }
         if (I == 352)
         {
	    PhoEner1_low = Beam_Energy - EnergyCal(P1[I - 1]);
	 }
	 else
	 {
	    PhoEner1_low = Beam_Energy - EnergyCal(fmin(P2[I + 1], P1[I - 1]));
	 }
	 PhoEner1_mean = 0.5 * (PhoEner1_high + PhoEner1_low);
	 Delta_E1 = PhoEner1_high - PhoEner1_low;

	 // Region corresponding to two neighbouring channels firing
	 if (I > 1)
	 {
	    PhoEner2_high = Beam_Energy - EnergyCal(P2[I]);
	    PhoEner2_low = Beam_Energy - EnergyCal(P2[I - 2]);
	    Delta_E2 = PhoEner2_high - PhoEner2_low;
	    
	    // Check that triple scintillator overlap is finite
	    if (Delta_E2 < 0.0)
	    {
	       Delta_E2 = 0.0;
            }
	 PhoEner2_mean = 0.5 * (PhoEner2_high + PhoEner2_low);
	 }
	 if (I > 1)
	 {
	    if(Delta_E2 == 0.0)
	    {
	       fprintf(fp,"%5d & %3d %12.3f %10.3f %12.3f %11.3f < no geometric overlap\n", I-2, I-1, PhoEner2_low, PhoEner2_high, PhoEner2_mean, Delta_E2);
	    }
	    else
	    {
	       fprintf(fp,"%5d & %3d %12.3f %10.3f %12.3f %11.3f\n", I-2, I-1, PhoEner2_low, PhoEner2_high, PhoEner2_mean, Delta_E2);
	    }
	 }
	 fprintf(fp,"%5d only %12.3f %10.3f %12.3f %11.3f %6c%2d %5c %5c%3d\n", I-1, PhoEner1_low, PhoEner1_high, PhoEner1_mean, Delta_E1, ecl[I-1], ECL, hv[I-1],sn[I-1],2*I-1);
	 J++;
      } 
      else if (Option == 3)
      {
	 // Region containing both single and double output channels
	 PhoEner_high = Beam_Energy - EnergyCal(P2[I]);
	 if (I == 352)
	 {
	    PhoEner_low = Beam_Energy - EnergyCal(P1[I - 1]);
	 }
	 else
	 {
	    PhoEner_low = Beam_Energy - EnergyCal(fmin(P2[I + 1], P1[I - 1]));
	 }
	 PhoEner_mean = 0.5 * (PhoEner_high + PhoEner_low);
	 Delta_E = PhoEner_high - PhoEner_low;
	 fprintf(fp,"%5d only%12.3f %10.3f %12.3f %11.3f %6c%2d %5c %5c%3d\n", I-1, PhoEner_low, PhoEner_high, PhoEner_mean, Delta_E, ecl[I-1], ECL, hv[I-1], sn[I-1], 2*I-1);
      }
      else if (Option == 4)
      {
	 // Write out upper and lower photon energies for regions of
	 // double (only) scintillator overlap (1 goosy channel
	 // firing) using true scintillator widths and geometries

	 // Region corresponding to single output channel firing only
	 if (I == 1)
	 {
	    PhoEner1_high = Beam_Energy - EnergyCal(P2[I]);
	 }
	 else
	 {
	    PhoEner1_high = Beam_Energy - EnergyCal(fmax(P1[I - 2], P2[I]));
	 }
	 if (I == 352)
	 {
	    PhoEner1_low = Beam_Energy - EnergyCal(P1[I - 1]);
	 }
	 else
	 {
	    PhoEner1_low = Beam_Energy - EnergyCal(fmin(P2[I + 1], P1[I - 1]));
	 }
	 PhoEner1_mean = 0.5 * (PhoEner1_high + PhoEner1_low);
	 Delta_E1 = PhoEner1_high - PhoEner1_low;
	 fprintf(fp,"%5d only%12.3f %10.3f %12.3f %11.3f %6c%2d %5c %5c%3d\n", I, PhoEner1_low, PhoEner1_high, PhoEner1_mean, Delta_E1, ecl[I-1], ECL, hv[I-1], sn[I-1], 2*I+1);
      }
      else if (Option == 5)
      {
	 // Write out upper and lower photon energies for regions of
	 // double scintillator overlap (1 goosy channel alone or
	 // with neightbours firing) using true scintillator widths
	 // and geometries
	 if (I == 0)
	 {
	    // Region corresponding to two neightbouring channels
	    // firing
	    PhoEner2_high = Beam_Energy - EnergyCal(P2[I]);
	    PhoEner2_low = Beam_Energy - EnergyCal(P1[I - 1]);
	    if ((BNMR >= 1.8767) && (BNMR <= 1.9146))
	    {
	       // PHENOMENOLOGICAL CORRECTION for 1.89563 Tesla, 
	       // 1557 MeV 8/2/09 jcm to correct measured deviations
	       // by quadratic fn bkII,p69 to correct measured
	       // deviations see p163, linear up to goosy 163
	       // quadratic 163-352, REVISED from p171, 18/2/08
	       pch = I - 0.5;
	       corrn = pch * pch * 7.387262E-05 - pch * 5.614872E-03 - 1.777049;
	       if (I == 1)
	       {
		  printf("Using phenomenological correction derived for 1557 MeV beam\n");
		  fprintf(fp,"As field was within 1 percent of 1.89563T, used phemomenological correction derived for 1557 MeV beam\n");
	       } 
	       PhoEner2_high = PhoEner2_high + corrn;
	       PhoEner2_low = PhoEner2_low + corrn;
	    }
	    if ((BNMR >= 1.81566) && (BNMR <= 1.85234))
	    {
	       // PHENOMENOLOGICAL CORRECTION for 1.834 Tesla, 1508
	       // MeV 14/2/08 jcm to correct measured deviations see
	       // p163, linear up to goosy 163 quadratic 163-352,
	       // REVISED from p171, 18/2/08
	       pch = I - 0.5;
	       if(pch <= 205.0)
	       {
		  corrn = 0.007506531 * pch - 1.774604;
	       }
	       else
	       {
		  corrn = pch * pch * 1.868485E-04 - pch * 6.894372E-02 + 6.0488;
	       }
	       if (I == 1)
	       {
		  printf("Using phenomenological correction derived for 1508 MeV beam\n");
		  fprintf(fp,"As field was within 1 percent of 1.834T, used phemomenological correction derived for 1508 MeV beam\n");
	       }
	    }
	    if ((BNMR >= 1.42857) && (BNMR <= 1.45743))
	    {
	       // PHENOMONOLOGICAL CORRECTION for 1.443 Tesla, 1204
	       // MeV 2/2/09 jcm to correct measured deviations see 
	       // bk II p56,57 with quad fn p61
	       pch = I - 0.5;
	       corrn = -pch * pch * 1.099716E-04 + pch * 2.2467585E-02 - 0.384447;
	       if(I == 1)
	       {
		  printf("Using phenomenological correction derived for 1204 MeV beam\n");
		  fprintf(fp,"As field was within 1 percent of 1.443T, used phemomenological correction derived for 1204 MeV beam\n");
	       }
	       PhoEner2_high = PhoEner2_high + corrn;
	       PhoEner2_low = PhoEner2_low + corrn;
	    }
	    if ((BNMR >= 1.00000) && (BNMR <= 1.07000))
	    {
	       // PHENOMENOLOGICAL CORRECTION for 1.00 - 1.070 Tesla,
	       // 840 - 883 MeV 28/3/13 jcm from calibrations done
               // for 885 MeV mostly by D. Middleton A2 Int Rept
               // 12/1 coefficients from DM email
               // DM855_fit_details.txt 6/2/12

	       // pch changed from 1204, 1508, 1557 MeV above as DM
	       // fits are vs goosy not the goosy - 0.5 space of
               // JCM's ROOT plots.
	       pch = I;
	       corrn = 0.00000;
	       if (pch <= 143.0)
	       {
		  corrn = 0.920235;
	       }
	       else if ((pch > 143.0) && (pch <= 246.0))
	       {
		  corrn = -pch * 0.0155779 + 9.16701;
	       }
	       else if ((pch > 246.0) && (pch <= 274.0))
	       {
		  corrn = -pch * 0.0439358 + 10.1328;
	       }
	       else // (pch > 274.0)
	       {
		  corrn = pch * 0.0173007 - 6.70907;
	       }
	       // Scale with factor EUF/(EUF for 1.02157T, for which
	       // DM's corrns apply)
	       corrn = corrn * B1 / 1.034413579;
	       if (I == 2)
	       {
		  printf("Using PHENOMENOLOGICAL correction derived for 855 MeV beam\n");
		  printf(" Goosy ch, Raw correction, Corrn scaled by EUF\n ");
		  fprintf(fp,"As field was in 1.00-1.07T range, used phemomenological EUF scaled correction, DM rept 2012/1 for 855 MeV beam\n");
	       }
	       // unsure if should be here now
	       // printf("pch, 1.034413579 * corrn / B1, corrn");
	       PhoEner2_high = PhoEner2_high + corrn;
	       PhoEner1_low = PhoEner2_low + corrn;
	    }
	    Delta_E2 = PhoEner2_high - PhoEner2_low;
	    
	    // Check that triple scintialltor overlap is finite
	    if (Delta_E2 < 0.0)
	    {
	       Delta_E2 = 0.0;
	    }
	    PhoEner2_mean = 0.5 * (PhoEner2_high + PhoEner2_low);
	 }
	 fprintf(fp,"%5d s+bnd%12.3f %10.3f %12.3f %11.3f %6c%2d %5c %5c%3d\n", I-1, PhoEner2_low, PhoEner2_high, PhoEner2_mean, Delta_E2, ecl[I-1], ECL, hv[I-1], sn[I-1], 2*I-1);
      }   
   }	 
   J = J + 1; 	// Increment Line Counter
   if (J >= 52)
   {
      J = 0;
   }
}
int result = fclose(fp);
}// End of ugcalv3_01


/*******************************************************************/

int main()
{
ugcalv3_01();
return 0;
}
