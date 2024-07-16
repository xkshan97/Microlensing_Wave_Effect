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
#include <numeric>
#include "../spline.h"
#include <random>
#include <string>
#include"./GetPsi_micro_field.h" //包含头文件
#include"./CosMoCal.h" //包含头文件
#include"./SampleMethod/RejectAndAcceptSample.h" //包含头文件

using std::sin;
using std::cos;
using std::tan;
using std::exp; 
using std::pow;
using std::log;
using std::sqrt;
using std::atan;



using namespace std;

#define PI acos(-1)

double MainDiffraction(double kappa, double gamma, double kappaStar_in, double LensRedshift, double SourceRedshift, int ThreadCount);