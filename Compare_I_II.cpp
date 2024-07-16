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
#include"./GetPsi_micro_field.h" //包含头文件
#include"./CosMoCal.h" //包含头文件
#include"./SampleMethod/RejectAndAcceptSample.h" //包含头文件
#include"./Micro_field_adaptive.h" //包含头文件

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

int main()
{
    //读取透镜的配置
    for(int subImageIndex = 0; subImageIndex < 50; subImageIndex ++ )
    {
        MainDiffraction(0.5967729088673934, 0.5967729088673934, 1.1 * 0.2786377498329086, 0.5, subImageIndex, 250);
        MainDiffraction(0.4038037265674532, 0.4038037265674532, 1.1 * 0.27675318987558506, 0.5, subImageIndex, 250);
    }
}