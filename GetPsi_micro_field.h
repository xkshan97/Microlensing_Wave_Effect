#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <time.h>
#include <functional>
#include <algorithm>
#include <ctime>
#include <vector>
#include <thread>
#include <string>
#include <random>
#include "../spline.h"

using std::sin;
using std::cos;
using std::tan;
using std::exp; 
using std::pow;
using std::log;
using std::sqrt;
using std::atan;
using std::acos;

using namespace std;


#define PI acos(-1)

double* Preparation4CreatPhiKappaStar(double kappa, double gamma, double kappaStar, double coeffi, int PrecisionFactor);
double* CreatMicroLens(double SkyLimitX, double SkyLimitY, int NStar);
double MicroCreatPsi(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double X1, double X2, double *MicroLensCoorXY, int NStar, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal);
long* NearFieldMicroIndex(double SkyLimitX, double SkyLimitY, int NStar, vector<vector<vector<long>>>& X1X2SetL1IncludeMicroIndex, double ResolutionL1, long X1L1Length, long X2L1Length, double X1L2Tmp, double X2L2Tmp);
double* FarFieldAlphaAndCenterPotential(double* MicroLensCoorXY, int NStar, long* ResNearFieldMicroIndex, double ResolutionL2, double X1L2Tmp, double X2L2Tmp, double* MassSample, double AveMassTotal);
double GetMinusPhi(double X1, double X2, double a1, double a2, double b1, double b2);
double GetMinusPhiDiffX1(double X1, double X2, double a1, double a2, double b1, double b2);
double GetMinusPhiDiffX2(double X1, double X2, double a1, double a2, double b1, double b2);
double MacroMicroAndMinusSheet(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double X1, double X2, double kappaStar, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal);
double* MacroMicroPartialTau(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal);
double MacroMicroFracK(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double InitialResolution, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal);
double MacroMicroDeltaTauInPixel(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double InitialResolution, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal);
int MacroMicroLayersNumber(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double InitialResolution, double EpsilonRes, double TimeDelayEsp, double coeffi, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal);
double* MacroMicroTheoryOut(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double Resolution, double coeffi, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal);