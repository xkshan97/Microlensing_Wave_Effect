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


//This is for multi micro lenses, one need to set the \kappa_*.
double* Preparation4CreatPhiKappaStar(double kappa, double gamma, double kappaStar, double coeffi, int PrecisionFactor)
{
    int MaxNStar = pow(10, 6)/2;
    double RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)});
    while(RNeed2Cal<50) //防止太小，这里设置的最小是50。
    {
        PrecisionFactor += 1;
        RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}); 
    }
    while(RNeed2Cal>900) //这里写的是900，为了留出来200的空隙（为的是与以前读数据的算法一样）
    {
        PrecisionFactor -= 1;
        RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}); 
    }
    // while((4 * RNeed2Cal * kappaStar / PI > pow(10, 6)) && (PrecisionFactor > 50)) //防止微透镜数目太多，这里选的是100万。
    // {
    //     PrecisionFactor -= 1;
    //     RNeed2Cal = PrecisionFactor*2*sqrt(kappaStar/PI)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}); 
    // }
    // cout << "SNR = 60 exceed sky map" << endl;
    // cout << "NOW SNR = " << PrecisionFactor << endl;
    double TimeNeed2CalMax = coeffi/2*min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)})*pow(PrecisionFactor*1.134*sqrt(kappaStar)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}),2);
    // cout << "TimeNeed2CalMax = " << TimeNeed2CalMax << endl;
    double MagnificationR = 1 - kappa + gamma;
    double MagnificationT = 1 - kappa - gamma;
    double SkyLimit;
    //Minimum和Maximum
    if(((MagnificationR > 0) && (MagnificationT > 0)) || ((MagnificationR < 0) && (MagnificationT < 0)))
    {
        double* PreOutPut = new double[6]; //这个是用来存后面得到的最大值和最小值的。
        double X10 = sqrt(abs(2*TimeNeed2CalMax/coeffi/MagnificationR));
        double X20 = sqrt(abs(2*TimeNeed2CalMax/coeffi/MagnificationT));
        SkyLimit = max({X10,X20}); 
        // if(SkyLimit<50)
        // {
        //     SkyLimit = 50;
        // }
        // cout << TmpX1 << TmpX2 << SkyLimit << endl;
        double AreaS = 4 * pow(SkyLimit, 2); //   其中0.1为类似于软化因子的东西，Skylimit是统计dA/dt的范围。撒恒星的范围要比他大0.1。
        int NStar = kappaStar * AreaS / PI;
        double kappaStarNew = NStar * PI / AreaS;
        //下面是为了防止微透镜数目过多，设置的最大容忍微透镜数目。
        while((NStar > MaxNStar) && (PrecisionFactor>20))
        {
            PrecisionFactor -= 1;
            TimeNeed2CalMax = coeffi/2*min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)})*pow(PrecisionFactor*1.134*sqrt(kappaStar)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}),2);
            X10 = sqrt(abs(2*TimeNeed2CalMax/coeffi/MagnificationR));
            X20 = sqrt(abs(2*TimeNeed2CalMax/coeffi/MagnificationT));
            SkyLimit = max({X10,X20}); 
            // if(SkyLimit<50)
            // {
            //     SkyLimit = 50;
            // }
            // cout << TmpX1 << TmpX2 << SkyLimit << endl;
            AreaS = 4 * pow(SkyLimit, 2); //   其中0.1为类似于软化因子的东西，Skylimit是统计dA/dt的范围。撒恒星的范围要比他大0.1。
            NStar = kappaStar * AreaS / PI;
            kappaStarNew = NStar * PI / AreaS;

        }
        cout << "NOW SNR = " << PrecisionFactor << endl;
        // cout << "Nstar = " << NStar << endl;
        // cout << "SkyLimit = " << SkyLimit << endl;
        PreOutPut[0] = SkyLimit;
        PreOutPut[1] = NStar;
        PreOutPut[2] = TimeNeed2CalMax;
        PreOutPut[3] = X20;
        PreOutPut[4] = X10;
        PreOutPut[5] = kappaStarNew;
        return PreOutPut;

    }
    //Saddle
    else if((MagnificationR > 0) && (MagnificationT < 0))
    {
        double* PreOutPut = new double[10]; //这个是用来存后面得到的最大值和最小值的。
        double SkyLimitX;
        double SkyLimitY;
        double X20, X10, X2010, X20Prime, X1020, X10Prime, TestX20X10, TestX20X10Max, X10New, X20New;
        double EllipticalA = sqrt(1/MagnificationR/coeffi);
        double EllipticalB = sqrt(-1/MagnificationT/coeffi); //FIXME如果不是magnificationT<0的话，这里需要改。
        double Epsilon1 =50;//50;//30;
        double Epsilon2 = EllipticalB/EllipticalA*Epsilon1; 
        // cout << "EllipticalA = " << EllipticalA << " " << "EllipticalB = " << EllipticalB << endl;
        // cout << "Epsilon for x20 = " << Epsilon1 << " " << "Epsilon for x10 = " << Epsilon2  << endl;

        X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
        X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
        //X20Prime = EllipticalB/EllipticalA*X2010;
        X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
        X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
        //X10Prime = EllipticalA/EllipticalB*X1020;
        SkyLimitY = X1020;//max({X20Prime, X1020});
        SkyLimitX = X2010;//max({X2010,X10Prime});
        TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
        TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 

        
        // cout << "X20 = " << X20 << " X2010 = " << X2010 << " X10 = " << X10 << " X1020 = " << X1020 << endl;
        // cout << "SkyLimitX = " << SkyLimitX << " SkyLimitY = " << SkyLimitY << endl;
        // cout << "SkyLimitY/SkyLimitX = " << SkyLimitY/SkyLimitX << " X20/X10 = " << X20/X10 << " EllipticalB/EllipticalA = " << EllipticalB/EllipticalA << endl;
        

        while((TestX20X10 < 10)&&(Epsilon1>=1)) //防止太小
        {

            Epsilon1 -= 0.5;
            Epsilon2 = EllipticalB/EllipticalA*Epsilon1;
            X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
            X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
            //X20Prime = EllipticalB/EllipticalA*X2010;
            X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
            X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
            //X10Prime = EllipticalA/EllipticalB*X1020;
            SkyLimitY = X1020;//max({X20Prime, X1020});
            SkyLimitX = X2010;//max({X2010,X10Prime});
            TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
            TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 
        }
        while(TestX20X10Max > 1000)
        {

            Epsilon1 += 0.5;
            Epsilon2 = EllipticalB/EllipticalA*Epsilon1;
            X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
            X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
            //X20Prime = EllipticalB/EllipticalA*X2010;
            X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
            X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
            //X10Prime = EllipticalA/EllipticalB*X1020;
            SkyLimitY = X1020;//max({X20Prime, X1020});
            SkyLimitX = X2010;//max({X2010,X10Prime});
            TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
            TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 
        }

        // cout << "After judge TestX20X10Max > 1000" << endl;
        // cout << "X20 = " << X20 << " X2010 = " << X2010 << " X10 = " << X10 << " X1020 = " << X1020 << endl;
        // cout << "SkyLimitX = " << SkyLimitX << " SkyLimitY = " << SkyLimitY << endl;
        // cout << "SkyLimitY/SkyLimitX = " << SkyLimitY/SkyLimitX << " X20/X10 = " << X20/X10 << " EllipticalB/EllipticalA = " << EllipticalB/EllipticalA << endl;
        X10New = X10;//pow(-2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*SkyLimitY*SkyLimitY,0.5);
        X20New = X20;//pow(-2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*SkyLimitX*SkyLimitX,0.5); 
        SkyLimit = max({SkyLimitX, SkyLimitY});
        // cout << "Final Epsilon" << endl;
        // cout << "Epsilon for X20 = " << Epsilon1 << " Epsilon for X10 = " << Epsilon2 << endl;

        double AreaS = 4 * SkyLimit * SkyLimit; //   其中0.1为类似于软化因子的东西，Skylimit是统计dA/dt的范围。撒恒星的范围要比他大0.1。
        int NStar = kappaStar * AreaS / PI;
        double kappaStarNew = NStar * PI / AreaS;
        //下面是为了防止微透镜数目过多，设置的最大容忍微透镜数目。
        while((NStar > MaxNStar) && (PrecisionFactor > 20))
        {
            PrecisionFactor -= 1;
            TimeNeed2CalMax = coeffi/2*min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)})*pow(PrecisionFactor*1.134*sqrt(kappaStar)/min({abs(1 - kappa + gamma), abs(1 - kappa - gamma)}),2);

            X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
            X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
            //X20Prime = EllipticalB/EllipticalA*X2010;
            X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
            X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
            //X10Prime = EllipticalA/EllipticalB*X1020;
            SkyLimitY = X1020;//max({X20Prime, X1020});
            SkyLimitX = X2010;//max({X2010,X10Prime});
            TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
            TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 

            
            // cout << "X20 = " << X20 << " X2010 = " << X2010 << " X10 = " << X10 << " X1020 = " << X1020 << endl;
            // cout << "SkyLimitX = " << SkyLimitX << " SkyLimitY = " << SkyLimitY << endl;
            // cout << "SkyLimitY/SkyLimitX = " << SkyLimitY/SkyLimitX << " X20/X10 = " << X20/X10 << " EllipticalB/EllipticalA = " << EllipticalB/EllipticalA << endl;
            

            while((TestX20X10 < 10)&&(Epsilon1>=1)) //防止太小
            {

                Epsilon1 -= 0.5;
                Epsilon2 = EllipticalB/EllipticalA*Epsilon1;
                X20 = EllipticalB*(-pow(Epsilon1,2) + 2 * pow(EllipticalA,2) * TimeNeed2CalMax)/(2 * EllipticalA * Epsilon1);
                X2010 = pow(2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*X20*X20,0.5);
                //X20Prime = EllipticalB/EllipticalA*X2010;
                X10 = - EllipticalA * Epsilon2 / 2 / EllipticalB + EllipticalA * EllipticalB * TimeNeed2CalMax / Epsilon2;  
                X1020 = pow(2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*X10*X10,0.5);
                //X10Prime = EllipticalA/EllipticalB*X1020;
                SkyLimitY = X1020;//max({X20Prime, X1020});
                SkyLimitX = X2010;//max({X2010,X10Prime});
                TestX20X10 = min({SkyLimitX, SkyLimitY, X20, X10});
                TestX20X10Max = max({SkyLimitX, SkyLimitY, X20, X10}); 
            }

            // cout << "After judge TestX20X10Max > 1000" << endl;
            // cout << "X20 = " << X20 << " X2010 = " << X2010 << " X10 = " << X10 << " X1020 = " << X1020 << endl;
            // cout << "SkyLimitX = " << SkyLimitX << " SkyLimitY = " << SkyLimitY << endl;
            // cout << "SkyLimitY/SkyLimitX = " << SkyLimitY/SkyLimitX << " X20/X10 = " << X20/X10 << " EllipticalB/EllipticalA = " << EllipticalB/EllipticalA << endl;
            X10New = X10;//pow(-2 * pow(EllipticalA,2) * TimeNeed2CalMax + pow(EllipticalA,2)/pow(EllipticalB,2)*SkyLimitY*SkyLimitY,0.5);
            X20New = X20;//pow(-2 * pow(EllipticalB,2) * TimeNeed2CalMax + pow(EllipticalB,2)/pow(EllipticalA,2)*SkyLimitX*SkyLimitX,0.5); 
            SkyLimit = max({SkyLimitX, SkyLimitY});
            // cout << "Final Epsilon" << endl;
            // cout << "Epsilon for X20 = " << Epsilon1 << " Epsilon for X10 = " << Epsilon2 << endl;

            AreaS = 4 * SkyLimit * SkyLimit; //   其中0.1为类似于软化因子的东西，Skylimit是统计dA/dt的范围。撒恒星的范围要比他大0.1。
            NStar = kappaStar * AreaS / PI;
            kappaStarNew = NStar * PI / AreaS;

        }
        cout << "NOW SNR = " << PrecisionFactor << endl;
        // cout << "Nstar = " << NStar << endl;
        PreOutPut[0] = SkyLimit;
        PreOutPut[1] = NStar;
        PreOutPut[2] = TimeNeed2CalMax;
        PreOutPut[3] = SkyLimitX;
        PreOutPut[4] = SkyLimitY;
        PreOutPut[5] = kappaStarNew;
        PreOutPut[6] = X10New;
        PreOutPut[7] = X20New; 
        PreOutPut[8] = Epsilon1;
        PreOutPut[9] = Epsilon2;
        return PreOutPut;
    }
    //未知
    else
    {
        cout << "Unknown Macrolensing Type" << endl;
        return 0;
    }
    
}


//微透镜场的
//微透镜坐标
double* CreatMicroLens(double SkyLimitX, double SkyLimitY, int NStar)
{
    double* MicroLensCoorXY = new double [NStar*2];
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> disx(-SkyLimitX, SkyLimitX);
    std::uniform_real_distribution<> disy(-SkyLimitY, SkyLimitY);
    for(long StarNum = 0; StarNum < NStar; StarNum ++ )
    {
        MicroLensCoorXY[2*StarNum] = disx(gen);
        MicroLensCoorXY[2*StarNum+1] = disy(gen);
        // MicroLensCoorXY[2*StarNum] = 0.1;//cos(PI/3);//disx(gen);
        // MicroLensCoorXY[2*StarNum+1] = 0.0;//sin(PI/3);//disy(gen);
        // MicroLensCoorXY[2*StarNum] = disx(gen);
        // MicroLensCoorXY[2*StarNum+1] = disy(gen);
        // cout << "x = " << MicroLensCoorX[StarNum] << " y = " << MicroLensCoorY[StarNum] << endl;
    }
    return MicroLensCoorXY;
    

}

//二级网格double X1L2Tmp, double X2L2Tmp处的近场星坐标。
long* NearFieldMicroIndex(double SkyLimitX, double SkyLimitY, int NStar, vector<vector<vector<long>>>& X1X2SetL1IncludeMicroIndex, double ResolutionL1, long X1L1Length, long X2L1Length, double X1L2Tmp, double X2L2Tmp)
{
    long X1L2TmpAtX1SetL1Index = (X1L2Tmp + SkyLimitX) / ResolutionL1;
    long X2L2TmpAtX2SetL1Index = (X2L2Tmp + SkyLimitY) / ResolutionL1;
    // cout << "X1L2TmpAtX1SetL1Index = " << X1L2TmpAtX1SetL1Index << endl;
    // cout << "X2L2TmpAtX2SetL1Index = " << X2L2TmpAtX2SetL1Index << endl;
    vector <long> ResNearFieldMicroIndexVec = {};
    for(long X1L1TmpIndex = (0 < X1L2TmpAtX1SetL1Index - 1 ? X1L2TmpAtX1SetL1Index - 1 : 0); X1L1TmpIndex <= (X1L1Length - 1 < X1L2TmpAtX1SetL1Index + 1 ? X1L1Length - 1 : X1L2TmpAtX1SetL1Index + 1); X1L1TmpIndex ++ )
    {
        for(long X2L1TmpIndex = (0 < X2L2TmpAtX2SetL1Index - 1 ? X2L2TmpAtX2SetL1Index - 1 : 0); X2L1TmpIndex <= (X2L1Length - 1 < X2L2TmpAtX2SetL1Index + 1 ? X2L1Length - 1 : X2L2TmpAtX2SetL1Index + 1) ; X2L1TmpIndex ++ )
        { 
            for(long testi = 0; testi < NStar && X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex][testi] != -1; testi ++ )
            {   /*将一级网格内保存的微透镜坐标取出，合并到周围九个大盒子中。*/
                ResNearFieldMicroIndexVec.push_back(X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex][testi]);
            }
        }
    }
    // for(long MicroTmpI = 0; MicroTmpI < NStar; MicroTmpI ++ )
    // {
    //     long TmpMicroAtX1SetL1 = MicroAtX1SetL1[MicroTmpI];
    //     long TmpMicroAtX2SetL1 = MicroAtX2SetL1[MicroTmpI];
    //     if((TmpMicroAtX1SetL1 >= X1L2TmpAtX1SetL1Index - 1) && (TmpMicroAtX1SetL1 <= X1L2TmpAtX1SetL1Index + 1) && (TmpMicroAtX2SetL1 >= X2L2TmpAtX2SetL1Index - 1) && (TmpMicroAtX2SetL1 <= X2L2TmpAtX2SetL1Index + 1))
    //     {
    //         ResNearFieldMicroIndexVec.push_back(MicroTmpI);
    //     }
    // }
    long ResNearFieldMicroIndexVecLength = ResNearFieldMicroIndexVec.size();
    // cout << "ResNearFieldMicroIndexVecLength = " << ResNearFieldMicroIndexVecLength << endl;
    long* ResNearFieldMicroIndex = new long [ResNearFieldMicroIndexVecLength + 1];
    ResNearFieldMicroIndex[0] = ResNearFieldMicroIndexVecLength; //第一个元素储存这个数组中近场星的个数。 
    for(long ResNearFieldMicroIndexVecIndex = 1; ResNearFieldMicroIndexVecIndex <= ResNearFieldMicroIndexVecLength; ResNearFieldMicroIndexVecIndex ++ )
    {
        ResNearFieldMicroIndex[ResNearFieldMicroIndexVecIndex] = ResNearFieldMicroIndexVec[ResNearFieldMicroIndexVecIndex - 1];
    }
    return ResNearFieldMicroIndex; 
}


/*计算远场在二级网格四个顶点的偏转角以及中心点的引力势*/
double* FarFieldAlphaAndCenterPotential(double* MicroLensCoorXY, int NStar, long* ResNearFieldMicroIndex, double ResolutionL2, double X1L2Tmp, double X2L2Tmp, double* MassSample, double AveMassTotal)
{
    double* ResFarFieldAlphaAndCenterPotential = new double [9];
    //中心点引力势
    double Psi = 0;
    //四个顶点的偏转角，从右上角A B C D逆时针旋转，与xuechun的PhD thesis一样。
    double AlphaA1 = 0;
    double AlphaA2 = 0;
    double AlphaB1 = 0;
    double AlphaB2 = 0; 
    double AlphaC1 = 0;
    double AlphaC2 = 0;
    double AlphaD1 = 0;
    double AlphaD2 = 0;

    double X1L2Low = X1L2Tmp - ResolutionL2 / 2;
    double X1L2Up = X1L2Tmp + ResolutionL2 / 2;
    double X2L2Low = X2L2Tmp - ResolutionL2 / 2;
    double X2L2Up = X2L2Tmp + ResolutionL2 / 2;
    
    
    //下面计算所有微透镜的贡献。
    for(long MicroTmpI = 0; MicroTmpI < NStar; MicroTmpI ++ )
    {
        double MicroX1 = MicroLensCoorXY[2*MicroTmpI];
        double MicroX2 = MicroLensCoorXY[2*MicroTmpI + 1];
        double denominatorA = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Up, 2);
        double denominatorB = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Up, 2);
        double denominatorC = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Low, 2);
        double denominatorD = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Low, 2);
        AlphaA1 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Up) / denominatorA;
        AlphaA2 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Up) / denominatorA;
        AlphaB1 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Low) / denominatorB;
        AlphaB2 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Up) / denominatorB;
        AlphaC1 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Low) / denominatorC;
        AlphaC2 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Low) / denominatorC;
        AlphaD1 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Up) / denominatorD;
        AlphaD2 += MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Low) / denominatorD;
        // double denominatorA = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Low, 2);
        // double denominatorB = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Low, 2);
        // double denominatorC = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Up, 2);
        // double denominatorD = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Up, 2);
        // AlphaA1 += (MicroX1 - X1L2Low) / denominatorA;
        // AlphaA2 += (MicroX2 - X2L2Low) / denominatorA;
        // AlphaB1 += (MicroX1 - X1L2Up) / denominatorB;
        // AlphaB2 += (MicroX2 - X2L2Low) / denominatorB;
        // AlphaC1 += (MicroX1 - X1L2Up) / denominatorC;
        // AlphaC2 += (MicroX2 - X2L2Up) / denominatorC;
        // AlphaD1 += (MicroX1 - X1L2Low) / denominatorD;
        // AlphaD2 += (MicroX2 - X2L2Up) / denominatorD;
        
        Psi += 0.5 * MassSample[MicroTmpI] / AveMassTotal * log(pow(MicroX1 - X1L2Tmp, 2) + pow(MicroX2 - X2L2Tmp, 2)); 
        // cout << Psi << endl;
    }
    //下面减去近场微透镜的贡献；
    for(long NearTmpI = 0; NearTmpI < ResNearFieldMicroIndex[0]; NearTmpI ++ )
    {
        long MicroTmpI = ResNearFieldMicroIndex[NearTmpI + 1];
        double MicroX1 = MicroLensCoorXY[2*MicroTmpI];
        double MicroX2 = MicroLensCoorXY[2*MicroTmpI + 1];
        double denominatorA = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Up, 2);
        double denominatorB = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Up, 2);
        double denominatorC = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Low, 2);
        double denominatorD = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Low, 2);
        AlphaA1 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Up) / denominatorA;
        AlphaA2 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Up) / denominatorA;
        AlphaB1 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Low) / denominatorB;
        AlphaB2 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Up) / denominatorB;
        AlphaC1 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Low) / denominatorC;
        AlphaC2 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Low) / denominatorC;
        AlphaD1 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX1 + X1L2Up) / denominatorD;
        AlphaD2 -= MassSample[MicroTmpI] / AveMassTotal * (- MicroX2 + X2L2Low) / denominatorD;
        // double denominatorA = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Low, 2);
        // double denominatorB = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Low, 2);
        // double denominatorC = pow(MicroX1 - X1L2Up, 2) + pow(MicroX2 - X2L2Up, 2);
        // double denominatorD = pow(MicroX1 - X1L2Low, 2) + pow(MicroX2 - X2L2Up, 2);
        // AlphaA1 -= (MicroX1 - X1L2Low) / denominatorA;
        // AlphaA2 -= (MicroX2 - X2L2Low) / denominatorA;
        // AlphaB1 -= (MicroX1 - X1L2Up) / denominatorB;
        // AlphaB2 -= (MicroX2 - X2L2Low) / denominatorB;
        // AlphaC1 -= (MicroX1 - X1L2Up) / denominatorC;
        // AlphaC2 -= (MicroX2 - X2L2Up) / denominatorC;
        // AlphaD1 -= (MicroX1 - X1L2Low) / denominatorD;
        // AlphaD2 -= (MicroX2 - X2L2Up) / denominatorD;
        
        Psi -= 0.5 * MassSample[MicroTmpI] / AveMassTotal * log(pow(MicroX1 - X1L2Tmp, 2) + pow(MicroX2 - X2L2Tmp, 2)); 
    }
    ResFarFieldAlphaAndCenterPotential[0] = Psi;
    ResFarFieldAlphaAndCenterPotential[1] = AlphaA1;
    ResFarFieldAlphaAndCenterPotential[2] = AlphaA2;
    ResFarFieldAlphaAndCenterPotential[3] = AlphaB1;
    ResFarFieldAlphaAndCenterPotential[4] = AlphaB2;
    ResFarFieldAlphaAndCenterPotential[5] = AlphaC1;
    ResFarFieldAlphaAndCenterPotential[6] = AlphaC2;
    ResFarFieldAlphaAndCenterPotential[7] = AlphaD1;
    ResFarFieldAlphaAndCenterPotential[8] = AlphaD2;
    return ResFarFieldAlphaAndCenterPotential;
    


}


//纯微透镜部分，计算引力势
double MicroCreatPsi(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double X1, double X2, double *MicroLensCoorXY, int NStar, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal)
{
    //先找到X1L3Tmp和X2L3Tmp所在的二级网格位置
    long X1L3TmpAtX1SetL2Index = (X1 + SkyLimitX) / ResolutionL2;
    long X2L3TmpAtX2SetL2Index = (X2 + SkyLimitY) / ResolutionL2;
    // cout << "X1L3TmpAtX1SetL2Index = " << X1L3TmpAtX1SetL2Index << endl;
    // cout << "X2L3TmpAtX2SetL2Index = " << X2L3TmpAtX2SetL2Index << endl;
    // cout << "X1SetL2 = " << X1SetL2[X1L3TmpAtX1SetL2Index] << endl;
    // cout << "X2SetL2 = " << X2SetL2[X2L3TmpAtX2SetL2Index] << endl;
    
    double Psi = 0.5 * ((1 - kappa + gamma) * X1 * X1 + (1 - kappa - gamma) * X2 * X2);
    // cout << "Galaxy = " << Psi << endl;
    // Psi = 0;
        
    //下面先加上近场微透镜的贡献；
    for(long NearTmpI = 0; NearTmpI < ResNearFieldMicroIndexSum[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][0]; NearTmpI ++ )
    {
        long MicroTmpI = ResNearFieldMicroIndexSum[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][NearTmpI + 1];
        double MicroX1 = MicroLensCoorXY[2*MicroTmpI];
        double MicroX2 = MicroLensCoorXY[2*MicroTmpI + 1];
        double denominator = pow(MicroX1 - X1, 2) + pow(MicroX2 - X2, 2);
        Psi -= 0.5 * MassSample[MicroTmpI] / AveMassTotal * log(denominator);
    }
    // cout << "Near micro = " << Psi << endl;
    // Psi = 0;
    //下面再加上远场微透镜的贡献
    
    //读取远场引力势和偏转角
    double PsiFar_0 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][0];
    double AlphaA1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][1];
    double AlphaA2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][2];
    double AlphaB1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][3];
    double AlphaB2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][4];
    double AlphaC1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][5];
    double AlphaC2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][6];
    double AlphaD1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][7];
    double AlphaD2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][8];

    double PsiFar_1, PsiFar_2, PsiFar_11, PsiFar_12, PsiFar_22, PsiFar_111, PsiFar_112, PsiFar_122, PsiFar_222, PsiFar_1111, PsiFar_1112, PsiFar_1122, PsiFar_1222, PsiFar_2222;
    double Delta_1 = X1 - X1SetL2[X1L3TmpAtX1SetL2Index];
    double Delta_2 = X2 - X2SetL2[X2L3TmpAtX2SetL2Index];
    double Delta_1_2 = pow(Delta_1, 2);
    double Delta_1_3 = pow(Delta_1, 3);
    double Delta_1_4 = pow(Delta_1, 4);
    double Delta_2_2 = pow(Delta_2, 2);
    double Delta_2_3 = pow(Delta_2, 3);
    double Delta_2_4 = pow(Delta_2, 4);

    if((Delta_1 > ResolutionL2) || (Delta_2 > ResolutionL2))
    {
        cout << "Error !!!!!!!!!! L3 not in L2" << endl;
    }
    // cout << "Delta 1 = " << Delta_1 << endl;
    // cout << "Delta 2 = " << Delta_2 << endl;
    

    PsiFar_1 = (double) 1/4. * (AlphaA1 + AlphaB1 + AlphaC1 + AlphaD1);
    PsiFar_2 = (double) 1/4 * (AlphaA2 + AlphaB2 + AlphaC2 + AlphaD2);
    PsiFar_11 = (double) 1/(4 * ResolutionL2) * (AlphaA1 - AlphaA2 - AlphaB1 - AlphaB2 - AlphaC1 + AlphaC2 + AlphaD1 + AlphaD2);
    PsiFar_22 = - PsiFar_11;
    PsiFar_12 = (double) 1/(4 * ResolutionL2) * (AlphaA1 + AlphaA2 + AlphaB1 - AlphaB2 - AlphaC1 - AlphaC2 - AlphaD1 + AlphaD2);
    PsiFar_111 = (double) 1/(ResolutionL2 * ResolutionL2) * (AlphaB2 - AlphaA2 + AlphaD2 - AlphaC2);
    PsiFar_122 = - PsiFar_111;
    PsiFar_112 = (double) 1/(ResolutionL2 * ResolutionL2) * (AlphaA1 - AlphaB1 + AlphaC1 - AlphaD1);
    PsiFar_222 = - PsiFar_112;
    PsiFar_1111 = (double) 3/(ResolutionL2 * ResolutionL2 * ResolutionL2) * ( - AlphaA1 - AlphaA2 + AlphaB1 - AlphaB2 + AlphaC1 + AlphaC2 - AlphaD1 + AlphaD2);
    PsiFar_1122 = - PsiFar_1111;
    PsiFar_2222 = PsiFar_1111;
    PsiFar_1112 = (double) 3/(ResolutionL2 * ResolutionL2 * ResolutionL2) * (AlphaA1 - AlphaA2 + AlphaB1 + AlphaB2 - AlphaC1 + AlphaC2 - AlphaD1 - AlphaD2);
    PsiFar_1222 = - PsiFar_1112;
    double PsiFar = (double) PsiFar_0 + Delta_1 * PsiFar_1 + Delta_2 * PsiFar_2 + 1/2. * PsiFar_11 * Delta_1_2 + PsiFar_12 * Delta_1 * Delta_2 + 1/2. * PsiFar_22 * Delta_2_2 \
                    + 1/6. * PsiFar_111 * Delta_1_3 + 1/2. * PsiFar_112 * Delta_1_2 * Delta_2 + 1/2. * PsiFar_122 * Delta_1 * Delta_2_2 + 1/6. * PsiFar_222 * Delta_2_3 \
                    + 1/24. * PsiFar_1111 * Delta_1_4 + 1/6. * PsiFar_1112 * Delta_1_3 * Delta_2 + 1/4. * PsiFar_1122 * Delta_1_2 * Delta_2_2 + 1/6. * PsiFar_1222 * Delta_1 * Delta_2_3 + 1/24. * PsiFar_2222 * Delta_2_4;
    Psi -= PsiFar;
    // cout << "Far micro = " << Psi << endl;
    // cout << "Psi0 = " << PsiFar_0 << endl;
    // cout << "AlphaA1 = " << AlphaA1 << endl;
    // cout << "AlphaA2 = " << AlphaA2 << endl;
    // cout << "AlphaB1 = " << AlphaB1 << endl;
    // cout << "AlphaB2 = " << AlphaB2 << endl;
    // cout << "AlphaC1 = " << AlphaC1 << endl;
    // cout << "AlphaC2 = " << AlphaC2 << endl;
    // cout << "AlphaD1 = " << AlphaD1 << endl;
    // cout << "AlphaD2 = " << AlphaD2 << endl;
    // cout << "PsiFar_1 = " << PsiFar_1 << endl;
    // cout << "PsiFar_2 = " << PsiFar_2 << endl;
    // cout << "PsiFar_11 = " << PsiFar_11 << endl;
    // cout << "PsiFar_12 = " << PsiFar_12 << endl;
    // cout << "PsiFar_111 = " << PsiFar_111 << endl;
    // cout << "PsiFar_112 = " << PsiFar_112 << endl;
    // cout << "PsiFar_1111 = " << PsiFar_1111 << endl;
    // cout << "PsiFar_1112 = " << PsiFar_1112 << endl;

    return Psi;

}



// //纯微透镜部分
// double MicroCreatPsi(double kappa, double gamma, double X1, double X2, double *MicroLensCoorXY, int NStar)
// {
//     //创建二维的数组用来存放势能的数据，并初始化为0；
//     double Psi = 0.5 * ((1 - kappa + gamma) * X1 * X1 + (1 - kappa - gamma) * X2 * X2);
        
//     for(long k = 0; k < NStar; k++)
//     {
//         double MicroX1 = MicroLensCoorXY[2*k];
//         double MicroX2 = MicroLensCoorXY[2*k + 1];
//         Psi -= 0.5 * log(pow(MicroX1 - X1, 2) + pow(MicroX2 - X2, 2));
//     }
           
//     return Psi;
// }



//负质量片
//a1最小值，a2最大值，b1最小值，b2最大值。
double GetMinusPhi(double X1, double X2, double a1, double a2, double b1, double b2)
{
    double p1, p2, p3, p4, p5, p6, p7, p8, p9;
    /*软化因子*/
    // a1 = a1 - 0.1;
    // a2 = a2 + 0.1;
    // b1 = b1 - 0.1;
    // b2 = b2 + 0.1;
    p1=(X2-b2)*(X1-a2)*log((X1-a2)*(X1-a2)+(X2-b2)*(X2-b2));
    p2=-(X2-b2)*(X1-a1)*log((X1-a1)*(X1-a1)+(X2-b2)*(X2-b2));
    p3=-(X2-b1)*(X1-a2)*log((X1-a2)*(X1-a2)+(X2-b1)*(X2-b1));
    p4=(X2-b1)*(X1-a1)*log((X1-a1)*(X1-a1)+(X2-b1)*(X2-b1));
    p5=(X2-b2)*(X2-b2)*(atan((X1-a2)/(X2-b2))-atan((X1-a1)/(X2-b2)));
    p6=-(X2-b1)*(X2-b1)*(atan((X1-a2)/(X2-b1))-atan((X1-a1)/(X2-b1)));
    p7=(X1-a2)*(X1-a2)*(atan((X2-b2)/(X1-a2))-atan((X2-b1)/(X1-a2)));
    p8=-(X1-a1)*(X1-a1)*(atan((X2-b2)/(X1-a1))-atan((X2-b1)/(X1-a1)));
    p9=3.*(b2-b1)*(a1-a2);
    return p1+p2+p3+p4+p5+p6+p7+p8+p9;
}

//负质量片对X1坐标求导
double GetMinusPhiDiffX1(double X1, double X2, double a1, double a2, double b1, double b2)
{
    double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14;
    p1 = -2*(a1 - X1)*atan((b1 - X2)/(a1 - X1));
    p2 = 2*(a2 - X1)*atan((b1 - X2)/(a2 - X1));
    p3 = 2*a1*atan((b2 - X2)/(a1 - X1));
    p4 = -2*X1*atan((b2 - X2)/(a1 - X1));
    p5 = -2*a2*atan((b2 - X2)/(a2 - X1));
    p6 = 2*X1*atan((b2 - X2)/(a2 - X1));
    p7 = -b1*log(pow(a1 - X1,2) + pow(b1 - X2,2));
    p8 = X2*log(pow(a1 - X1,2) + pow(b1 - X2,2));
    p9 = b1*log(pow(a2 - X1,2) + pow(b1 - X2,2));
    p10 = -X2*log(pow(a2 - X1,2) + pow(b1 - X2,2));
    p11 = b2*log(pow(a1 - X1,2) + pow(b2 - X2,2));
    p12 = -X2*log(pow(a1 - X1,2) + pow(b2 - X2,2));
    p13 = -b2*log(pow(a2 - X1,2) + pow(b2 - X2,2));
    p14 = X2*log(pow(a2 - X1,2) + pow(b2 - X2,2));
    return p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14;
}

//负质量片对X2坐标求导
double GetMinusPhiDiffX2(double X1, double X2, double a1, double a2, double b1, double b2)
{
    double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14;
    p1 = -2*(b1 - X2)*atan((a1 - X1)/(b1 - X2));
    p2 = 2*(b1 - X2)*atan((a2 - X1)/(b1 - X2));
    p3 = 2*b2*atan((a1 - X1)/(b2 - X2));
    p4 = -2*X2*atan((a1 - X1)/(b2 - X2));
    p5 = -2*b2*atan((a2 - X1)/(b2 - X2));
    p6 = 2*X2*atan((a2 - X1)/(b2 - X2));
    p7 = -a1*log(pow(a1 - X1,2) + pow(b1 - X2,2));
    p8 = X1*log(pow(a1 - X1,2) + pow(b1 - X2,2));
    p9 = a2*log(pow(a2 - X1,2) + pow(b1 - X2,2));
    p10 = -X1*log(pow(a2 - X1,2) + pow(b1 - X2,2));
    p11 = a1*log(pow(a1 - X1,2) + pow(b2 - X2,2));
    p12 = -X1*log(pow(a1 - X1,2) + pow(b2 - X2,2));
    p13 = -a2*log(pow(a2 - X1,2) + pow(b2 - X2,2));
    p14 = X1*log(pow(a2 - X1,2) + pow(b2 - X2,2));
    return p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14;
}

//强透镜加微透镜场加负质量片的引力势
double MacroMicroAndMinusSheet(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double X1, double X2, double kappaStar, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal)
{
    double Psi = MicroCreatPsi(SkyLimitX, SkyLimitY, kappa, gamma, X1, X2, MicroLensCoorXY, NStar, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
    // printf("Psi in MacroMicroAndMinusSheet = %.10f\n", Psi);
    double MinusSheet = kappaStar / 2 / PI * GetMinusPhi(X1, X2, a1, a2, b1, b2); 
    // printf("MinusSheet in MacroMicroAndMinusSheet = %.10f\n", MinusSheet);
    return Psi + MinusSheet;
}

// //强透镜加微透镜加负质量片的求导部分
// double* MacroMicroPartialTau(double kappa, double gamma, double kappaStar, double X1, double X2, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2)
// {
//     double* Res = new double [2]();
    
    
//     double PartialX1 = (1 - kappa + gamma) * X1;
//     double PartialX2 = (1 - kappa - gamma) * X2;

//     for(long k = 0; k < NStar; k++ )
//     {
//         double MicroX1 = MicroLensCoorXY[2*k];
//         double MicroX2 = MicroLensCoorXY[2*k + 1];
//         double denominator = pow(MicroX1 - X1, 2) + pow(MicroX2 - X2, 2);
//         PartialX1 += (MicroX1 - X1) / denominator;
//         PartialX2 += (MicroX2 - X2) / denominator;
//     }
    

//     Res[0] = PartialX1 + kappaStar / 2 / PI * GetMinusPhiDiffX1(X1, X2, a1, a2, b1, b2);
//     Res[1] = PartialX2 + kappaStar / 2 / PI * GetMinusPhiDiffX2(X1, X2, a1, a2, b1, b2);
//     return Res;
// }



//强透镜加微透镜加负质量片，计算偏转角 (强透镜加微透镜加负质量片的求导部分)
double* MacroMicroPartialTau(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal)
{
    //先找到X1L3Tmp和X2L3Tmp所在的二级网格位置
    // cout << "X1 = " << X1 << "X2 = " << X2 << endl;
    long X1L3TmpAtX1SetL2Index = (X1 + SkyLimitX) / ResolutionL2;
    long X2L3TmpAtX2SetL2Index = (X2 + SkyLimitY) / ResolutionL2;
    // cout << "X1L3TmpAtX1SetL2Index = " << X1L3TmpAtX1SetL2Index << endl;
    // cout << "X2L3TmpAtX2SetL2Index = " << X2L3TmpAtX2SetL2Index << endl;
    // cout << "X1SetL2 = " << X1SetL2[X1L3TmpAtX1SetL2Index] << endl;
    // cout << "X2SetL2 = " << X2SetL2[X2L3TmpAtX2SetL2Index] << endl;
    
    double AlphaX1 = (1 - kappa + gamma) * X1;
    double AlphaX2 = (1 - kappa - gamma) * X2;
    // cout << "Galaxy = " << Psi << endl;
    // Psi = 0;
        
    //下面先加上近场微透镜的贡献；
    for(long NearTmpI = 0; NearTmpI < ResNearFieldMicroIndexSum[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][0]; NearTmpI ++ )
    {
        long MicroTmpI = ResNearFieldMicroIndexSum[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][NearTmpI + 1];
        double MicroX1 = MicroLensCoorXY[2*MicroTmpI];
        double MicroX2 = MicroLensCoorXY[2*MicroTmpI + 1];
        double denominator = pow(MicroX1 - X1, 2) + pow(MicroX2 - X2, 2);
        AlphaX1 -= MassSample[MicroTmpI] / AveMassTotal * (X1 - MicroX1) / denominator;
        AlphaX2 -= MassSample[MicroTmpI] / AveMassTotal * (X2 - MicroX2) / denominator;
    }
    // cout << "Near micro = " << Psi << endl;
    // Psi = 0;
    //下面再加上远场微透镜的贡献
    
    //读取远场引力势和偏转角
    double PsiFar_0 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][0];
    double AlphaA1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][1];
    double AlphaA2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][2];
    double AlphaB1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][3];
    double AlphaB2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][4];
    double AlphaC1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][5];
    double AlphaC2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][6];
    double AlphaD1 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][7];
    double AlphaD2 = ResFarFieldAlphaAndCenterPotential[X1L3TmpAtX1SetL2Index][X2L3TmpAtX2SetL2Index][8];

    double PsiFar_1, PsiFar_2, PsiFar_11, PsiFar_12, PsiFar_22, PsiFar_111, PsiFar_112, PsiFar_122, PsiFar_222, PsiFar_1111, PsiFar_1112, PsiFar_1122, PsiFar_1222, PsiFar_2222;
    double Delta_1 = X1 - X1SetL2[X1L3TmpAtX1SetL2Index];
    double Delta_2 = X2 - X2SetL2[X2L3TmpAtX2SetL2Index];
    double Delta_1_2 = pow(Delta_1, 2);
    double Delta_1_3 = pow(Delta_1, 3);
    double Delta_2_2 = pow(Delta_2, 2);
    double Delta_2_3 = pow(Delta_2, 3);
    

    if((Delta_1 > ResolutionL2) || (Delta_2 > ResolutionL2))
    {
        cout << "Error !!!!!!!!!! L3 not in L2" << endl;
    }
    // cout << "Delta 1 = " << Delta_1 << endl;
    // cout << "Delta 2 = " << Delta_2 << endl;
    

    PsiFar_1 = (double) 1/4. * (AlphaA1 + AlphaB1 + AlphaC1 + AlphaD1);
    PsiFar_2 = (double) 1/4 * (AlphaA2 + AlphaB2 + AlphaC2 + AlphaD2);
    PsiFar_11 = (double) 1/(4 * ResolutionL2) * (AlphaA1 - AlphaA2 - AlphaB1 - AlphaB2 - AlphaC1 + AlphaC2 + AlphaD1 + AlphaD2);
    PsiFar_22 = - PsiFar_11;
    PsiFar_12 = (double) 1/(4 * ResolutionL2) * (AlphaA1 + AlphaA2 + AlphaB1 - AlphaB2 - AlphaC1 - AlphaC2 - AlphaD1 + AlphaD2);
    PsiFar_111 = (double) 1/(ResolutionL2 * ResolutionL2) * (AlphaB2 - AlphaA2 + AlphaD2 - AlphaC2);
    PsiFar_122 = - PsiFar_111;
    PsiFar_112 = (double) 1/(ResolutionL2 * ResolutionL2) * (AlphaA1 - AlphaB1 + AlphaC1 - AlphaD1);
    PsiFar_222 = - PsiFar_112;
    PsiFar_1111 = (double) 3/(ResolutionL2 * ResolutionL2 * ResolutionL2) * ( - AlphaA1 - AlphaA2 + AlphaB1 - AlphaB2 + AlphaC1 + AlphaC2 - AlphaD1 + AlphaD2);
    PsiFar_1122 = - PsiFar_1111;
    PsiFar_2222 = PsiFar_1111;
    PsiFar_1112 = (double) 3/(ResolutionL2 * ResolutionL2 * ResolutionL2) * (AlphaA1 - AlphaA2 + AlphaB1 + AlphaB2 - AlphaC1 + AlphaC2 - AlphaD1 - AlphaD2);
    PsiFar_1222 = - PsiFar_1112;
    double AlphaFarX1 = (double) PsiFar_1 + PsiFar_11 * Delta_1 + PsiFar_12 * Delta_2 \
                        + 1/2. * PsiFar_111 * (Delta_1_2 - Delta_2_2) + PsiFar_112 * Delta_1 * Delta_2 \
                        + 1/6. * PsiFar_1111 * (Delta_1_3 - 3 * Delta_1 * Delta_2_2) + 1/6. * PsiFar_1112 * (3 * Delta_1_2 * Delta_2 - Delta_2_3);
    double AlphaFarX2 = (double) PsiFar_2 + PsiFar_22 * Delta_2 + PsiFar_12 * Delta_1 \
                        + 1/2. * PsiFar_222 * (Delta_2_2 - Delta_1_2) + PsiFar_122 * Delta_1 * Delta_2 \
                        + 1/6. * PsiFar_2222 * (Delta_2_3 - 3 * Delta_2 * Delta_1_2) + 1/6. * PsiFar_1222 * (3 * Delta_2_2 * Delta_1 - Delta_1_3);
    AlphaX1 -= AlphaFarX1;
    AlphaX2 -= AlphaFarX2;
    // printf("AlphaX1 in MacroMicroPartialTau = %.10f\n", AlphaX1);
    // printf("AlphaX2 in MacroMicroPartialTau = %.10f\n", AlphaX2);
    //加上负质量片
    AlphaX1 += kappaStar / 2 / PI * GetMinusPhiDiffX1(X1, X2, a1, a2, b1, b2);
    AlphaX2 += kappaStar / 2 / PI * GetMinusPhiDiffX2(X1, X2, a1, a2, b1, b2);
    double* ResAlpha = new double [2];
    ResAlpha[0] = AlphaX1;
    ResAlpha[1] = AlphaX2;
    // cout << "Far micro = " << Psi << endl;
    // cout << "Psi0 = " << PsiFar_0 << endl;
    // cout << "AlphaA1 = " << AlphaA1 << endl;
    // cout << "AlphaA2 = " << AlphaA2 << endl;
    // cout << "AlphaB1 = " << AlphaB1 << endl;
    // cout << "AlphaB2 = " << AlphaB2 << endl;
    // cout << "AlphaC1 = " << AlphaC1 << endl;
    // cout << "AlphaC2 = " << AlphaC2 << endl;
    // cout << "AlphaD1 = " << AlphaD1 << endl;
    // cout << "AlphaD2 = " << AlphaD2 << endl;
    // cout << "PsiFar_1 = " << PsiFar_1 << endl;
    // cout << "PsiFar_2 = " << PsiFar_2 << endl;
    // cout << "PsiFar_11 = " << PsiFar_11 << endl;
    // cout << "PsiFar_12 = " << PsiFar_12 << endl;
    // cout << "PsiFar_111 = " << PsiFar_111 << endl;
    // cout << "PsiFar_112 = " << PsiFar_112 << endl;
    // cout << "PsiFar_1111 = " << PsiFar_1111 << endl;
    // cout << "PsiFar_1112 = " << PsiFar_1112 << endl;

    return ResAlpha;

}


//获得像素内梯度误差比
double MacroMicroFracK(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double InitialResolution, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal)
{
    double* X1Five = new double [5]();
    double* X2Five = new double [5]();
    X1Five[0] = X1 - InitialResolution / 2;
    X1Five[1] = X1 + InitialResolution / 2;
    X1Five[2] = X1 + InitialResolution / 2;
    X1Five[3] = X1 - InitialResolution / 2;
    X1Five[4] = X1;

    X2Five[0] = X2 - InitialResolution / 2;
    X2Five[1] = X2 - InitialResolution / 2;
    X2Five[2] = X2 + InitialResolution / 2;
    X2Five[3] = X2 + InitialResolution / 2;
    X2Five[4] = X2;

   
    double* GradRes = new double [5]();
    double* PartialResX1 = new double [5]();
    double* PartialResX2 = new double [5]();
    double* GradResDiff = new double [10]();
    double AveGrad = 0;
    for(int i = 0; i < 5; i ++ )
    {
        double* PartialResTmp = MacroMicroPartialTau(SkyLimitX, SkyLimitY, kappa, gamma, kappaStar, X1Five[i], X2Five[i], MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
        GradRes[i] = sqrt(PartialResTmp[0] * PartialResTmp[0] + PartialResTmp[1] * PartialResTmp[1]);
        PartialResX1[i] = PartialResTmp[0];
        PartialResX2[i] = PartialResTmp[1];
        AveGrad += GradRes[i] / 5;
        delete [] PartialResTmp;
    }
    int TmpIndex = 0;
    for(int i = 0; i < 4; i ++ )
    {
        for(int j = i + 1; j < 5; j ++ )
        {
            GradResDiff[TmpIndex] = sqrt((PartialResX1[j] - PartialResX1[i]) * (PartialResX1[j] - PartialResX1[i]) + (PartialResX2[j] - PartialResX2[i]) * (PartialResX2[j] - PartialResX2[i]));
            TmpIndex += 1;
        }

    }
    double MaxGradDiff = 0;
   
    for(int i = 0; i < 10; i ++ )
    {
        if(GradResDiff[i] > MaxGradDiff)
        {
            MaxGradDiff = GradResDiff[i];
        }
    }
   
    double Fraction = MaxGradDiff / AveGrad;
    // cout << "Average = " << AvePartial << endl;
    // cout << "Maximum = " << MaxPartialDiff << endl;
    delete[] X1Five;
    delete[] X2Five;
    delete[] GradRes;
    delete[] PartialResX1;
    delete[] PartialResX2;
    delete[] GradResDiff;
    return Fraction;

}



//获得像素内最大的时间差
double MacroMicroDeltaTauInPixel(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double InitialResolution, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal)
{
    double* X1Five = new double [5]();
    double* X2Five = new double [5]();
    X1Five[0] = X1 - InitialResolution / 2;
    X1Five[1] = X1 + InitialResolution / 2;
    X1Five[2] = X1 + InitialResolution / 2;
    X1Five[3] = X1 - InitialResolution / 2;
    X1Five[4] = X1;

    X2Five[0] = X2 - InitialResolution / 2;
    X2Five[1] = X2 - InitialResolution / 2;
    X2Five[2] = X2 + InitialResolution / 2;
    X2Five[3] = X2 + InitialResolution / 2;
    X2Five[4] = X2;
    

    double* TauTmp = new double [5]();
  
    for(int i = 0; i < 5; i ++ )
    {
        TauTmp[i] = MacroMicroAndMinusSheet(SkyLimitX, SkyLimitY, kappa, gamma, X1Five[i], X2Five[i], kappaStar, MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
    }
    double MaxDeltaTau = 0;
    for(int i = 0; i < 4; i ++ )
    {
        for(int j = i + 1; j < 5; j ++ )
        {
            double TauTmpTmp = abs(TauTmp[i] - TauTmp[j]);
            if(TauTmpTmp > MaxDeltaTau)
            {
                MaxDeltaTau = TauTmpTmp;
            }
        }

    }
   
    delete[] X1Five;
    delete[] X2Five;
    delete[] TauTmp;
   

    return MaxDeltaTau;

}


//获得加密层数
int MacroMicroLayersNumber(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double InitialResolution, double EpsilonRes, double TimeDelayEsp, double coeffi, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal)
{
    int LayersMax = 7;
    double Fraction = MacroMicroFracK(SkyLimitX, SkyLimitY, kappa, gamma, kappaStar, X1, X2, InitialResolution, MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
    double MaxDeltaTau = MacroMicroDeltaTauInPixel(SkyLimitX, SkyLimitY, kappa, gamma, kappaStar, X1, X2, InitialResolution, MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
    // cout << "Sgn = " << Sgn << endl;
    //本来想粗略估计一个，但是发现效果并不好，而且还有bug。
    // double NablaTauX1 = 2 * (1 - kappa + gamma) * X1;
    // double NablaTauX2 = 2 * (1 - kappa - gamma) * X2;
    // double NablaNablaTauX1 = 2 * (1 - kappa + gamma);
    // double NablaNablaTauX2 = 2 * (1 - kappa - gamma);
    // double TheoryDeltaL = EpsilonRes * sqrt(NablaTauX1 * NablaTauX1 + NablaTauX2 * NablaTauX2) / sqrt(NablaNablaTauX1 * NablaNablaTauX1 + NablaNablaTauX2 * NablaNablaTauX2);
    // int NumLayers = InitialResolution / TheoryDeltaL;
    // // cout << "NumLayers = " << NumLayers << endl; 
    // if(NumLayers >= LayersMax)
    // {
    //     NumLayers = LayersMax;
    // }
    // if(NumLayers != 0)
    // {
    //     NumLayers -= 1;
    // }
     
    int NumLayers = 0;
    double X1Low = X1 - InitialResolution / 2;
    double X1High = X1 + InitialResolution / 2;
    double X2Low = X2 - InitialResolution / 2;
    double X2High = X2 + InitialResolution / 2;
    double X1MaxFrac, X2MaxFrac;
    // for(int NumLayersTmp = 0; NumLayersTmp < LayersMax; NumLayersTmp ++ )
    // while(((Fraction > EpsilonRes)&&(NumLayers <= LayersMax)))
    while(((Fraction > EpsilonRes)&&(MaxDeltaTau * coeffi > 0.2 * TimeDelayEsp)&&(NumLayers < LayersMax)))
    {
        // if((Fraction > EpsilonRes)&&(MaxDeltaTau * coeffi > 0.2 * TimeDelayEsp))
        // {
        double* X1Center = new double [(int) pow(4, NumLayers)]();
        double* X2Center = new double [(int) pow(4, NumLayers)]();
        int IndexX1Circle = 0;
        int IndexX2Circle = 0;
        int SumIndex = 0;
        while(SumIndex < pow(4, NumLayers))
        {
            if(IndexX1Circle < pow(2, NumLayers))
            {
                X1Center[SumIndex] = X1Low + InitialResolution / pow(2, NumLayers + 1) + IndexX1Circle * InitialResolution / pow(2, NumLayers);
                X2Center[SumIndex] = X2Low + InitialResolution / pow(2, NumLayers + 1) + IndexX2Circle * InitialResolution / pow(2, NumLayers);
                IndexX1Circle += 1;
                SumIndex += 1;
            }
            else
            {
                IndexX1Circle = 0;
                IndexX2Circle += 1;
            }
        }
        double FractionMax = 0;
        double MaxDeltaTauMax = 0;
        for(int i = 0; i < pow(4, NumLayers); i ++ )
        {
            double FractionTmp = MacroMicroFracK(SkyLimitX, SkyLimitY, kappa, gamma, kappaStar, X1Center[i], X2Center[i], InitialResolution / pow(2, NumLayers), MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
            double MaxDeltaTauTmp = MacroMicroDeltaTauInPixel(SkyLimitX, SkyLimitY, kappa, gamma, kappaStar, X1Center[i], X2Center[i], InitialResolution / pow(2, NumLayers), MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
            if(FractionTmp > FractionMax)
            {
                FractionMax = FractionTmp;
                X1MaxFrac = X1Center[i];
                X2MaxFrac = X2Center[i];
            }
            if(MaxDeltaTauTmp > MaxDeltaTauMax)
            {
                MaxDeltaTauMax = MaxDeltaTauTmp;
            }
            
            
        }
        
        Fraction = FractionMax;
        MaxDeltaTau = MaxDeltaTauMax;
        
        // cout << "Fraction = " << Fraction << endl;
        NumLayers += 1;
        // cout << "NumLayers: " << NumLayers << endl;
        delete[] X1Center;
        delete[] X2Center;
        // }
        // else
        // {
        //     break;
        // }

        
        
    }
    if((NumLayers != 0) && (NumLayers != LayersMax))
    {
        NumLayers -= 1;
    }

    // if(NumLayers == LayersMax)
    // {
        // cout << "Num Layers May not be enough!!!!!!!!!" << endl;
        // cout << "X1 :" << X1 << endl;
        // cout << "X2 :" << X2 << endl;
        // cout << "Fraction :" << Fraction << endl;
        // cout << "Sgn :" << Sgn << endl;
        // cout << "Max Frac X1 :"<< X1MaxFrac << endl;
        // cout << "Max Frac X2 :"<< X2MaxFrac << endl;
        // printf("Sgn Wrong X1 :%2.12lf\n", X1SgnWrong);
        // printf("Sgn Wrong X2 :%2.12lf\n", X2SgnWrong);
    // }
    return NumLayers;
}


//梯度近似算法的输出
double* MacroMicroTheoryOut(double SkyLimitX, double SkyLimitY, double kappa, double gamma, double kappaStar, double X1, double X2, double Resolution, double coeffi, double *MicroLensCoorXY, int NStar, double a1, double a2, double b1, double b2, long*** ResNearFieldMicroIndexSum, double ResolutionL2, double*** ResFarFieldAlphaAndCenterPotential, double* X1SetL2, double* X2SetL2, double* MassSample, double AveMassTotal)
{
    double* PartialRes = MacroMicroPartialTau(SkyLimitX, SkyLimitY, kappa, gamma, kappaStar, X1, X2, MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
    double Grads = sqrt(PartialRes[0] * PartialRes[0] + PartialRes[1] * PartialRes[1]);
    double Alpha = acos((PartialRes[0] + PartialRes[1])/(sqrt(2)*Grads));
    if((Alpha >= 0)&(Alpha < PI/4))
    {
        Alpha = Alpha;
    }
    else if((Alpha < PI/2)&&(Alpha >= PI/4))
    {
        Alpha = PI/2 - Alpha;
    }
    else if((Alpha >= PI/2)&&(Alpha < 3*PI/4))
    {
        Alpha = Alpha - PI/2;
    }
    else if((Alpha >= 3*PI/4)&&(Alpha < PI))
    {
        Alpha = PI - Alpha;
    }

    // cout << "alpha = " << Alpha << endl;
    
    double t0 = MacroMicroAndMinusSheet(SkyLimitX, SkyLimitY, kappa, gamma, X1, X2, kappaStar, MicroLensCoorXY, NStar, a1, a2, b1, b2, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal) * coeffi;
    double L1 = sqrt(2)/2 * Resolution * cos(Alpha);
    double L2 = sqrt(2)/2 * Resolution * sin(Alpha);
    double t1 = t0 - L2 * Grads * coeffi;
    double t2 = t0 + L2 * Grads * coeffi;
    double t3 = t0 - L1 * Grads * coeffi;
    double t4 = t0 + L1 * Grads * coeffi;
    double* TSUM = new double[5];
    TSUM[0] = t1;
    TSUM[1] = t2;
    TSUM[2] = t3;
    TSUM[3] = t4;
    // cout << "t0 = " << t0 << endl;
    // cout << "t1 = " << t1 << endl;
    // cout << "t2 = " << t2 << endl;
    // cout << "t3 = " << t3 << endl;
    // cout << "t4 = " << t4 << endl;
    double Platform = Resolution / cos(PI/4 - Alpha) / Grads / coeffi;
    TSUM[4] = Platform;
    delete[] PartialRes;
    // cout << "Platform = " << Platform << endl;
    return TSUM;
} 



// 储存微透镜场验证一下。
// int main()
// {
//     double* ReadMap = new double [20000*20000];
//     int SizeLen = sizeof(double);
//     ifstream PsiFile1;
//     PsiFile1.open("/data/pubdata/Data_cxc_to_sxk/saddle_large/phi_stars_and_sheet_500.bin", std::ifstream::binary);
//     PsiFile1.read(reinterpret_cast<char *>(ReadMap), SizeLen*20000*20000); //读整个Map
//     PsiFile1.close ();

//     double* ReadLensPos = new double [19099*2];
//     ifstream PosFile1;
//     PosFile1.open("/data/pubdata/Data_cxc_to_sxk/saddle_large/lens_pos_19099.bin", std::ifstream::binary);
//     PosFile1.read(reinterpret_cast<char *>(ReadLensPos), SizeLen*19099*2); //读lensing position
//     PosFile1.close ();
//     cout << "X[0] = " << ReadLensPos[0] << "Y[0] = " << ReadLensPos[1] << endl;
    
//     printf("Micro and Minus%10.8f \n", 0.5 * ((49.95 * 49.95) + (49.95 * 49.95)) - MacroMicroAndMinusSheet(0., 0, 19099 * PI / 1000 / 1000, -49.95, -49.95, ReadLensPos, 19099, -500, 500, -500, 500));
//     vector<double> X1Set;
//     int X1Len = 0;
//     for(double tmp_x1 = -50 + 0.05 / 2; tmp_x1 <= 50; tmp_x1 += 0.05)
//     {
//         X1Set.push_back(tmp_x1);
//         X1Len += 1;
//     }
//     vector<double> X2Set = X1Set;
//     cout << "X1Len = X1Len = " << X1Len << endl;
    
//     double* MicroMinusPhi = new double[2000*2000];
//     int ThreadCount = 50; //使用的线程数量。

//     double** MicroAndMinus2D = new double*[ThreadCount];
//     for(int i = 0; i < ThreadCount; i ++ )
//     {
//         MicroAndMinus2D[i] = new double[2000*40];
//     }

//     std::thread* threads[ThreadCount];
//     int args[ThreadCount];
    

//     //倍乘因子，主要是算一下一个核算多少行。
//     long Order = X1Len/ThreadCount; 

    


//     auto task = [&] (int i) -> void {
//         int tmp = 0;
//         long IndexUpLimit = i + 1 == ThreadCount ? X1Len : (i + 1) * Order; //三目运算符
//         for(long index_x1 = i*Order; index_x1 < IndexUpLimit; index_x1++ )
//         {
//             for(int index_x2 = 0; index_x2 < X1Len; index_x2 ++)
//             {
//                 MicroAndMinus2D[i][tmp * X1Len + index_x2] = 0.5 * ((X1Set[index_x1] * X1Set[index_x1]) + (X2Set[index_x2] * X2Set[index_x2])) - MacroMicroAndMinusSheet(0., 0, 19099 * PI / 1000 / 1000, X1Set[index_x1], X2Set[index_x2], ReadLensPos, 19099, -500, 500, -500, 500);
//             }
//             tmp += 1;
//             // cout << tmp << endl;
//         }
//         if(tmp != Order)
//         {
//             cout << "tmp Error!!! tmp = " << tmp << endl;
//         }
    
//     };


//     auto run_tasks = [&] (int count) -> void {
//         for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
//             threads[CountIndex] = new std::thread (task, args[CountIndex]);
//         }
//         for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
//             threads[CountIndex]->join();
//             delete threads[CountIndex]; 
//         }
//     };

//     long JThread = 0;
//     for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
//     {
//         args[JThread++] = ThreadIndex;
//         if (JThread == ThreadCount) {
//             run_tasks(ThreadCount);
//             JThread = 0;
//         }
//     }
//     run_tasks(JThread);
    
//     for(int test_i = 0; test_i < ThreadCount; test_i ++)
//     {
//         for(int test_j = 0; test_j < 2000 * 40; test_j ++)
//         {
//             MicroMinusPhi[test_i * 2000 * 40 + test_j] = MicroAndMinus2D[test_i][test_j];    
//         }
        
//     }        


//     ofstream Creat_Micro_Minus;
//     Creat_Micro_Minus.open("./Creat_Micro_Minus.bin", std::ifstream::binary);
//     Creat_Micro_Minus.write(reinterpret_cast<char *>(MicroMinusPhi), sizeof(double)*2000*2000);
//     Creat_Micro_Minus.close();
    
//     return 0;


// }

// //验证近似公式
// int main()
// {
   
//     double coeffi = 0.0016737719515190538;
//     char TFour[100];
//     sprintf(TFour,"Test/TFour.bin");
//     ofstream TFour_file;
//     TFour_file.open(TFour, std::ofstream::binary);

//     char DSDT[100];
//     sprintf(DSDT,"Test/DSDT.bin");
//     ofstream DSDT_file;
//     DSDT_file.open(DSDT, std::ofstream::binary);

//     char TimeName[100];
//     sprintf(TimeName,"Test/Time.bin"); 
//     ofstream fp_time;
//     fp_time.open(TimeName, std::ofstream::binary);

//     // double* SolveResult = SolveLensEquSIS(-0.1, 0.1);
//     // for(int tmpi = 0; tmpi < 4; tmpi ++ )
//     // {
//     //     cout << SolveResult[tmpi] << endl;
//     // }
//     // double Xa = sqrt(SolveResult[0] * SolveResult[0] + SolveResult[1] * SolveResult[1]);
//     // double Xb = sqrt(SolveResult[2] * SolveResult[2] + SolveResult[3] * SolveResult[3]);
//     double kappa = 0.7;
//     double gamma = -0.25;
//     double X1 = 10.1;
//     double X2 = 10;
//     double ImageResolution = 0.05;
//     double TimeDelayEsp = pow(10,-6);
//     double ImageResolutionBoost = ImageResolution / 100;
//     // int test = SameSgn(kappa, gamma, X1, X2, ImageResolution);
//     // double Test = FracK(0.6, 0.3, 20, 20, 0.1);
//     // int NumLayers = LayersNumber(0.6, 0.3, 20, 20, 0.1, 0.1);

//     int NStar = 19099;
//     double* ReadLensPos = new double [NStar*2];
//     double kappaStar = NStar * PI / 1000 / 1000;
//     ifstream PosFile1;
//     PosFile1.open("/data/pubdata/Data_cxc_to_sxk/saddle_large/lens_pos_19099.bin", std::ifstream::binary);
//     PosFile1.read(reinterpret_cast<char *>(ReadLensPos), 8*19099*2); //读lensing position
//     PosFile1.close ();
//     cout << "X[0] = " << ReadLensPos[0] << "Y[0] = " << ReadLensPos[1] << endl;
    


//     double DeltaT = MacroMicroDeltaTauInPixel(kappa, gamma, kappaStar, X1, X2, ImageResolution, ReadLensPos, NStar, -500, 500, -500, 500);
//     double* PartialRes = MacroMicroPartialTau(kappa, gamma, kappaStar, X1, X2, ReadLensPos, NStar, -500, 500, -500, 500);
//     double Grads = sqrt(PartialRes[0] * PartialRes[0] + PartialRes[1] * PartialRes[1]);
//     double Alpha = acos((PartialRes[0] + PartialRes[1])/(sqrt(2)*Grads));
//     cout << "Initial Alpha:" << Alpha << endl;
//     if((Alpha >= 0)&(Alpha < PI/4))
//     {
//         Alpha = Alpha;
//     }
//     else if((Alpha < PI/2)&&(Alpha >= PI/4))
//     {
//         Alpha = PI/2 - Alpha;
//     }
//     else if((Alpha >= PI/2)&&(Alpha < 3*PI/4))
//     {
//         Alpha = Alpha - PI/2;
//     }
//     else if((Alpha >= 3*PI/4)&&(Alpha < PI))
//     {
//         Alpha = PI - Alpha;
//     }

//     cout << "alpha = " << Alpha << endl;
//     double t0 = MacroMicroAndMinusSheet(kappa, gamma, kappaStar, X1, X2, ReadLensPos, NStar, -500, 500, -500, 500) * coeffi;
//     double L1 = sqrt(2)/2 * ImageResolution * cos(Alpha);
//     double L2 = sqrt(2)/2 * ImageResolution * sin(Alpha);
//     double t1 = t0 - L2 * Grads * coeffi;
//     double t2 = t0 + L2 * Grads * coeffi;
//     double t3 = t0 - L1 * Grads * coeffi;
//     double t4 = t0 + L1 * Grads * coeffi;
//     double TSUM[6];
//     TSUM[0] = t1;
//     TSUM[1] = t2;
//     TSUM[2] = t3;
//     TSUM[3] = t4;
//     cout << "t0 = " << t0 << endl;
//     cout << "t1 = " << t1 << endl;
//     cout << "t2 = " << t2 << endl;
//     cout << "t3 = " << t3 << endl;
//     cout << "t4 = " << t4 << endl;
//     double Platform = ImageResolution / cos(PI/4 - Alpha) / Grads / coeffi;
//     TSUM[4] = Platform;
//     cout << "Platform = " << Platform << endl;
//     double* Test = MacroMicroTheoryOut(kappa, gamma, kappaStar, X1, X2, ImageResolution, coeffi, ReadLensPos, NStar, -500, 500, -500, 500);
//     cout << Test[0] << " " << Test[1] << " " << Test[2] << " " << Test[3] << " " << Test[4] << endl;

//     // cout << "Number Layers = " << NumLayers << endl;
//     // cout << "Fraction = " << Test << endl;
//     cout << "Delta t = " << DeltaT << endl;





    
    
//     /*天区坐标*/
//     long LengthX = (ImageResolution + ImageResolutionBoost)/ImageResolutionBoost - 1; 
//     long LengthY = LengthX;
//     double* X1Set = new double [LengthX];
//     long X1Index = 0;
    
    
//     for(double TmpX1 = X1 - ImageResolution/2 + ImageResolutionBoost/2; TmpX1 < X1 + ImageResolution/2; TmpX1 += ImageResolutionBoost)
//     {
//         if(X1Index < LengthX)
//         {
//             X1Set[X1Index] = TmpX1;
//             X1Index ++ ;
//         }
//         else
//         {
//             X1Index ++ ;
//         }
//     } 

//     cout << "X1's length is right "<< endl;
//     cout << "Length of X1 set = " << LengthX << " =  Maximum X1 set index + 1 = " << X1Index << endl; 
    
   
//     double* X2Set = new double [LengthX];
//     long X2Index = 0;
//     for(double TmpX2 = X2 - ImageResolution/2 + ImageResolutionBoost/2; TmpX2 < X2 + ImageResolution/2; TmpX2 += ImageResolutionBoost)
//     {
//         if(X2Index < LengthX)
//         {
//             X2Set[X2Index] = TmpX2;
//             X2Index ++ ;
//         }
//         else
//         {
//             X2Index ++ ;
//         }
//     } 
    
//     printf("X1Set[0] = %2.10lf\n", X1Set[0]);
//     printf("X1Set[-1] = %2.10lf\n", X1Set[LengthX - 1]);
//     printf("X2Set[0] = %2.10f\n", X2Set[0]); 
//     printf("X2Set[-1] = %2.10f\n", X2Set[LengthX - 1]);


//     double RudeImageResolution; // SkyLimit, x2: lens plane coordinate. 
//     //y1, y2: source plane coordinate. esp: lens plane resolution

//     long RudeEspBoost = 10;
//     RudeImageResolution = ImageResolution * RudeEspBoost;

    

//     /*下面是粗略计算最大值点和最小值点的子程序*/
//     double* RudeTimeDelayMinAndMax = new double[2]; //这个是用来存后面得到的最大值和最小值的。
//     RudeTimeDelayMinAndMax[0] = 10000;
//     RudeTimeDelayMinAndMax[1] = -10000;
//     double TestX1, TestX2;
//     long TestX1Index, TestX2Index;
//     for(long RudeX1Index = 0; RudeX1Index < LengthX; RudeX1Index += RudeEspBoost )
//     {
//         double RudeX1Tmp = X1Set[RudeX1Index];

//         // cout << RudeX1Tmp << endl;
//         for(long RudeY1Index = 0; RudeY1Index < LengthX; RudeY1Index += RudeEspBoost )
//         {
//             double RudeX2Tmp = X2Set[RudeY1Index];
//             double RudePsi = MacroMicroAndMinusSheet(kappa, gamma, kappaStar, RudeX1Tmp, RudeX2Tmp, ReadLensPos, NStar, -500, 500, -500, 500);
            
//             if(RudePsi > RudeTimeDelayMinAndMax[1])
//             {
//                 RudeTimeDelayMinAndMax[1] = RudePsi;
//                 TestX1 = X1Set[RudeX1Index];
//                 TestX2 = X2Set[RudeY1Index];
                
//             }
//             if(RudePsi < RudeTimeDelayMinAndMax[0])
//             {
//                 RudeTimeDelayMinAndMax[0] = RudePsi;
//                 TestX1Index = RudeX1Index;
//                 TestX2Index = RudeY1Index;
//                 // cout << "Test = " << X1Set[TestX1Index] << endl;
                
//             }
        
//         }
        
//     }
//     double RudeTimeDelayMin = RudeTimeDelayMinAndMax[0] * coeffi - 0.1;
//     double RudeTimeDelayMax = RudeTimeDelayMinAndMax[1] * coeffi + 0.1;
//     cout << "Rude Minimum Time Delay is " << RudeTimeDelayMin << endl;
//     cout << "Rude Maximum X1 = " << TestX1 << " X2 = " << TestX2 << endl;
//     cout << "Rude Minimum X1 Index = " << TestX1Index << " X2 Index = " << TestX2Index << endl; 
//     cout << "Rude Maximun Time Delay is " << RudeTimeDelayMax << endl;
//     delete[] RudeTimeDelayMinAndMax;



//     //the next is for time range
//     double TimeDelayStart = RudeTimeDelayMin; 
//     double TimeDelayEnd = RudeTimeDelayMax;
//     //
//     cout << "Time Maximum Boundary For Time Delay Curve " << TimeDelayEnd << endl;
//     long LenTimeDelay = (TimeDelayEnd - TimeDelayStart + TimeDelayEsp)/TimeDelayEsp; 
//     // the length of time delay area
//     cout << "Time Delay Curve's Length = " << LenTimeDelay << endl; //print len_time
//     vector <long> TimeLength;
//     TimeLength.push_back(LenTimeDelay);
    
    
//     double* TimeDelayRange = new double [LenTimeDelay];
//     double TimeDelayTmp = TimeDelayStart; //the start of the time delay calculation.
//     //generate time array.
//     for(long TimeDelayIndex = 0; TimeDelayIndex < LenTimeDelay; TimeDelayIndex ++ )
//     {
    
//         TimeDelayRange[TimeDelayIndex] = TimeDelayTmp;
//         TimeDelayTmp += TimeDelayEsp;
        
        

//     }
//     cout << TimeDelayRange[LenTimeDelay - 1] << endl;
//     //cout << int (sizeof(time_range)/sizeof(time_range[0])) << endl;


//     double* TimeDelayArea = new double [LenTimeDelay](); //存放时间延迟面积的数组。

//     for(long Index = 0; Index < LengthX; Index++ )
//     { 
//         double TmpX1 = X1Set[Index]; //算一下这一列的x值。
        
        
//         // double* TimeDelay = CreatPsi(kappa, gamma, TmpX1, MicroLensCoorXY, NStar, SkyLimit, LengthX, X1Set, kappaStar, y1, y2);
//         for(long JIndex = 0; JIndex < LengthY; JIndex += 1)
//         {
//             double TmpX2 = X2Set[JIndex]; //算一下坐标的y值。
//             //判断是不是坐标是不是在双曲线区域内。
//             // if(1)
//             double TimeDelay = MacroMicroAndMinusSheet(kappa, gamma, kappaStar, TmpX1, TmpX2, ReadLensPos, NStar, -500, 500, -500, 500) * coeffi;

//             long TimeDelayIndex = long ((TimeDelay - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) - 1;
//             // cout << time_index << endl;
            
//             TimeDelayArea[TimeDelayIndex] += double(ImageResolutionBoost * ImageResolutionBoost) / TimeDelayEsp; //FIXME 这里可能有错误，比如说多个线程同时相加。
        
//         }
//     }
//     // for(int i = 0; i < LenTimeDelay; i ++ )
//     // {
//     //     cout << TimeDelayArea[i] << ",";
//     // }
//     TSUM[5] = LenTimeDelay;
//     TFour_file.write(reinterpret_cast<char *>(TSUM), sizeof(double)*6);

//     DSDT_file.write(reinterpret_cast<char *>(TimeDelayArea), sizeof(double)*LenTimeDelay);
//     fp_time.write(reinterpret_cast<char *>(TimeDelayRange), sizeof(double)*LenTimeDelay);
    
//     TFour_file.close();
//     DSDT_file.close();
//     fp_time.close();

//     return 0;
// }