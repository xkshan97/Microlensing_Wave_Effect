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
#include "../spline.h"
#include <random>
#include <string>


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


double HubbleParam(double z)
{
    double Omegam = 0.30966;
    double Omegav = 0.69034;
    double H0 = 67.66;
    return H0 * sqrt(Omegam * pow(1 + z, 3) + Omegav);
}

double AngularDiameterDis(double zStart, double zEnd)
{
    vector<double> zset;
    double zstep = 0.0001;
    double Integrad = 0;
    double lightspeed = 2.9979246 * pow(10,5);
    for(double zi = zStart; zi <= zEnd; zi += zstep)
    {
        zset.push_back(zi);
        Integrad += 1/HubbleParam(zi);
        // cout << zi << endl;
    }
    return Integrad * zstep * lightspeed / (1 + zEnd);
}





// int main()
// {
//     double DL = AngularDiameterDis(0, 1);
//     cout << DL << endl;
// }