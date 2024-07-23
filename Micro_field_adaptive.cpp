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
#include "./spline.h"
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

double MainDiffraction(double kappa, double gamma, double kappaStar_in, double LensRedshift, double SourceRedshift, int ThreadCount)
// int main()
{
    // for(int TotalIndex = 0; TotalIndex < 100; TotalIndex ++ )
    // {
    //总线程数
    // int ThreadCount = 300; 
    time_t TimeStart;
    TimeStart = time(NULL); //start time.
    //时间分辨率
    double TimeDelayEsp = pow(10,-6);

    //Macro配置
    // double kappa = 0.8; //4.250996065574553; //1.2;
    // double gamma = 0.25; //0.25; //1.2136718743511312; //-0.15;
    //Trapezoid approximation precision
    double EpsilonRes = 0.1;

    //Micro配置
    //物理常数以及单位变换
    double M_sun = 1.9884099 * pow(10,30);
    double G = 6.6743 * pow(10,-11);
    double c = 2.9979246 * pow(10,8);
    //微透镜场的配置
    //time length/sky boundary SNR
    int PrecisionFactor = 100;
    int Max_mass = 15; //只是管着文件名的，88，100，28确实代表最大质量，但是50代表的是Minimum和Saddle分别做了50个随机实现。

    // double LensRedshift = 0.5;
    // double SourceRedshift = 1;
    /*粗略计算平均质量，这里的平均质量是根据PDF手算的，跟恒星IMF模型有关系，跟最小最大质量也有关系*/
    /*下面用的是Chariber IMF (0.08~1.5) + IFM*/
    string IMFType = "Chabrier";
    cout << "IMF = " << IMFType << endl;
    double AveMassStellar = 0.3;
    double AveMassRemnant = 1.3974668046341392;
    double Remnant2Stellar = 0.2; //Remnant占Stellar的总量。一般是10%。
    double AveMassTotal = (AveMassRemnant * AveMassStellar + AveMassRemnant * AveMassStellar * Remnant2Stellar) / (AveMassRemnant + AveMassStellar * Remnant2Stellar);
    // double MassPerSolarMass = 1;
    double coeffi = 4 * G * AveMassTotal * M_sun * (1 + LensRedshift) / pow(c, 3);
    printf("Average LensMass from PDF:%2.20f\n", AveMassTotal);
    
    // double kappaStar_in = 0.06;
    //calculate boundary parameters, star total number, maximum time delay.
    double* PreOutPut = Preparation4CreatPhiKappaStar(kappa, gamma, kappaStar_in, coeffi, PrecisionFactor);
    for(int tmp_print_i = 0; tmp_print_i < 6; tmp_print_i ++ )
    {
        cout << "PreOutPut" << " " << PreOutPut[tmp_print_i] << endl;
    }

    double SkyLimit_Micro_file = PreOutPut[0]; //微透镜坐标的边界
    int NStar = PreOutPut[1]; //(2 * SkyLimit_Micro_file) * (2 * SkyLimit_Micro_file) * kappaStar_in / PI;
    /*根据总数目以及平均质量，能算出总质量，再根据质量守恒，计算出各组分的质量*/
    int NStarStellar = (double) NStar * AveMassTotal / (1 + Remnant2Stellar) / AveMassStellar;
    int NStarRemnant = (double) NStarStellar * AveMassStellar * Remnant2Stellar / AveMassRemnant;
    cout << "NStarStellar = " << NStarStellar << endl;
    cout << "NStarRemnant = " << NStarRemnant << endl;
    double TimeNeed2CalMax = PreOutPut[2];
    NStar = NStarRemnant + NStarStellar;
    double* MassSample = SampleResult(NStarStellar, NStarRemnant, IMFType);
    double SumMass = 0;
    for(int SumMassIndex = 0; SumMassIndex < NStar; SumMassIndex ++ )
    {
        SumMass += MassSample[SumMassIndex];
    }
    AveMassTotal = SumMass / NStar;
    coeffi = 4 * G * AveMassTotal * M_sun * (1 + LensRedshift) / pow(c, 3);
    printf("Numerical Average LensMass from Sample:%2.20f\n", AveMassTotal);
    printf("Numerical Maximum LensMass from Sample:%2.20f\n", *max_element(MassSample, MassSample + NStar));
    double kappaStar = SumMass / AveMassTotal * PI / (2 * SkyLimit_Micro_file) / (2 * SkyLimit_Micro_file);
    printf("Numerical exact kappaStar from Sample = %.10f\n", kappaStar);
    

    double* MicroLensCoorXY = CreatMicroLens(SkyLimit_Micro_file, SkyLimit_Micro_file, NStar);
    cout << "Micro Pos X[0] = " << MicroLensCoorXY[0] << " " << "Micro Pos Y[0] = " << MicroLensCoorXY[1] << endl;
    cout << "Minimum Pos XY = " << (*min_element(MicroLensCoorXY, MicroLensCoorXY+2*NStar)) << " " << "Maximum Pos XY = " << (*max_element(MicroLensCoorXY, MicroLensCoorXY+2*NStar)) << endl;
    
    //储存微透镜的质量和坐标
    //质量
    ofstream IMFSampleTestOP;
    char IMFSampleTestFile[100];
    sprintf(IMFSampleTestFile,"./MicroField_%2d/Lens_Mass_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
    IMFSampleTestOP.open(IMFSampleTestFile, std::ofstream::binary);
    IMFSampleTestOP.write(reinterpret_cast<char *>(MassSample), sizeof(double)*(NStarStellar + NStarRemnant));
    IMFSampleTestOP.close();
    //坐标
    ofstream MicroLensCoorXYOP;
    char MicroLensCoorXYFile[100];
    sprintf(MicroLensCoorXYFile,"./MicroField_%2d/MicroLensCoorXY_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
    MicroLensCoorXYOP.open(MicroLensCoorXYFile, std::ofstream::binary);
    MicroLensCoorXYOP.write(reinterpret_cast<char *>(MicroLensCoorXY), sizeof(double)*2*(NStarStellar + NStarRemnant));
    MicroLensCoorXYOP.close();

    
    double mur = 1 - kappa + gamma;
    double mut = 1 - kappa - gamma;
    
    
    
    /*下面设置1级网格储存微透镜的位置*/
    /*与xuechun的不一样，我这里直接计算每个微透镜所占的正方形盒子大小*/
    double ResolutionL1 = 2 * 2 * SkyLimit_Micro_file / sqrt(NStar) < 2 * SkyLimit_Micro_file / 10 ? 2 * 2 * SkyLimit_Micro_file / sqrt(NStar) : 2 * SkyLimit_Micro_file / 10;
    cout << "ResolutionL1 = " << ResolutionL1 << endl;
    vector <double> X1SetL1 = {};
    for(double X1L1Tmp = - SkyLimit_Micro_file + ResolutionL1 / 2; X1L1Tmp < SkyLimit_Micro_file + ResolutionL1; X1L1Tmp += ResolutionL1)
    {
        X1SetL1.push_back(X1L1Tmp);
    }
    vector <double> X2SetL1 = X1SetL1;
    long X1L1Length = X1SetL1.size();
    long X2L1Length = X2SetL1.size();
    cout << "X1L1Length = " << X1L1Length << endl;
    //创建一个X1L1Length宽，X2L1Length高，至少一个元素-1的三维数组，用来储存每个一级网格里的微透镜索引。
    vector<vector<vector<long>>> X1X2SetL1IncludeMicroIndex(X1L1Length, vector<vector<long>>(X2L1Length));
   



    for(long X1L1TmpIndex = 0; X1L1TmpIndex < X1L1Length; X1L1TmpIndex ++ )
    {
        for(long X2L1TmpIndex = 0; X2L1TmpIndex < X2L1Length; X2L1TmpIndex ++ )
        {
            X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex].insert(X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex].begin(), -1); 
        }
    }
    // long* MicroAtX1SetL1 = new long [NStar]();
    // long* MicroAtX2SetL1 = new long [NStar]();
    // long* MicroIndexFake = new long [NStar](); //这是一个假的Index，为的是下面判断微透镜在哪个盒子里时，一但已经分配了格子，就不要再继续运行了。
    for(long MicroTmpI = 0; MicroTmpI < NStar; MicroTmpI ++ )
    {
        bool MicroFoundLoc = false;  // 添加标志表示微透镜是否找到网格。
        for(long X1L1TmpIndex = 0; X1L1TmpIndex < X1L1Length && !MicroFoundLoc; X1L1TmpIndex ++ )
        {
            for(long X2L1TmpIndex = 0; X2L1TmpIndex < X2L1Length && !MicroFoundLoc; X2L1TmpIndex ++ )
            {
                if((abs(MicroLensCoorXY[2 * MicroTmpI] - X1SetL1[X1L1TmpIndex]) <= ResolutionL1 / 2) && (abs(MicroLensCoorXY[2 * MicroTmpI + 1] - X2SetL1[X2L1TmpIndex]) <= ResolutionL1 / 2))
                {
                    X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex].insert(X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex].begin(), MicroTmpI);
                    MicroFoundLoc = true;
                }
            }
        }
    }
    long SumNumMicro = 0; //检验是否所有的微透镜都分配到网格里了
    for(long X1L1TmpIndex = 0; X1L1TmpIndex < X1L1Length; X1L1TmpIndex ++ )
    {
        for(long X2L1TmpIndex = 0; X2L1TmpIndex < X2L1Length; X2L1TmpIndex ++ )
        { 
            for(long testi = 0; testi < NStar && X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex][testi] != -1; testi ++ )
            {
                // cout << "X1X2SetL1IncludeMicroIndex[0][0] = " << X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex][testi] << endl;
                // cout << "MicroLensCoorX[] = " << MicroLensCoorXY[2 * X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex][testi]]<< endl;
                // cout << "MicroLensCoorY[] = " << MicroLensCoorXY[2 * X1X2SetL1IncludeMicroIndex[X1L1TmpIndex][X2L1TmpIndex][testi] + 1]<< endl;
                SumNumMicro += 1;
            }
        }
    }
    if(SumNumMicro == NStar)
    {
        cout << "All Micro lenses have been assigned to a L1 grid !!!!!!"<< endl;
        cout << "SumNumMicro = " << SumNumMicro << endl;
        cout << "NStar = " << NStar << endl;
    }
    else
    {
        cout << "Wrong!!!!!!!! Not All Micro lenses have been assigned to a L1 grid !!!!!!"<< endl;
        cout << "SumNumMicro = " << SumNumMicro << endl;
        cout << "NStar = " << NStar << endl;
    }
    
    
    // cout << "MicroAtX1SetL1[0] = " << MicroAtX1SetL1[0] << endl;
    // cout << "MicroAtX2SetL1[0] = " << MicroAtX2SetL1[0] << endl;
    /*验证一下在第一个网格内的微透镜*/
    // for(long MicroTmpI = 0; MicroTmpI < NStar; MicroTmpI ++ )
    // {
    //     if(MicroAtX1SetL1[MicroTmpI] == 0 && MicroAtX2SetL1[MicroTmpI] == 0)
    //     {
    //         cout << "MicroTmpI = " << MicroTmpI << endl;
    //         cout << "MicroCoorXY[0] = " << MicroLensCoorXY[2 * MicroTmpI] << endl;
    //         cout << "MicroCoorXY[1] = " << MicroLensCoorXY[2 * MicroTmpI + 1] << endl;
    //     }

    // }
    /*下面设置二级网格的分辨率和二级网格坐标*/
    double ResolutionL2 = ResolutionL1 / 20;
    cout << "ResolutionL2 = " << ResolutionL2 << endl;
    vector <double> X1SetL2Vec = {};
    for(double X1L2Tmp = - SkyLimit_Micro_file + ResolutionL2 / 2; X1L2Tmp < SkyLimit_Micro_file + ResolutionL2; X1L2Tmp += ResolutionL2)
    //注意，这里也在上界处加了一个像素，但是是L2的像素不是L1的像素。这么做的原因，有两点：1、为的是保证下一级中的网格中，-limit-limit之间的所有光线，都有二级网格可寻。2、二级网格不会超过1级网格。
    {
        X1SetL2Vec.push_back(X1L2Tmp);
    }
    vector <double> X2SetL2Vec = X1SetL2Vec;
    long X1L2Length = X1SetL2Vec.size();
    long X2L2Length = X2SetL2Vec.size();
    cout << "X1L2Length = " << X1L2Length << endl;
    double* X1SetL2 = new double [X1L2Length];
    double* X2SetL2 = new double [X2L2Length];
    for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex ++ )
    {
        X1SetL2[X1L2TmpIndex] = X1SetL2Vec[X1L2TmpIndex];
    }
    for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex ++ )
    {
        X2SetL2[X2L2TmpIndex] = X2SetL2Vec[X2L2TmpIndex];
    }
    printf("X2SetL2 boundary = %.10f\n", X2SetL2[X2L2Length - 1] + ResolutionL2 / 2);
    printf("X2SetL1 boundary = %.10f\n", X2SetL1[X2L1Length - 1] + ResolutionL1 / 2);
     
    /*输出结果判断是否正确*/
    if((X1SetL2Vec[X1L2Length - 1] != X1SetL2[X1L2Length - 1]) || (X2SetL2Vec[X2L2Length - 1] != X2SetL2[X2L2Length - 1]) || (X1SetL1[X1L1Length - 1] < X1SetL2[X1L2Length - 1]) || (X2SetL1[X2L1Length - 1] < X2SetL2[X2L2Length - 1]))
    {
        cout << "L1 or L2 grid wrong !!! " << endl;
        cout << "X1SetL1Vec[-1]" << X1SetL1[X1L1Length - 1] << endl;
        cout << "X2SetL1Vec[-1]" << X2SetL1[X2L1Length - 1] << endl;
        cout << "X1SetL2Vec[-1]" << X1SetL2Vec[X1L2Length - 1] << endl;
        cout << "X2SetL2Vec[-1]" << X2SetL2Vec[X2L2Length - 1] << endl;
        cout << "X1SetL2[-1]" << X1SetL2[X1L2Length - 1] << endl;
        cout << "X2SetL2[-1]" << X2SetL2[X2L2Length - 1] << endl;
    }
    else
    {
        cout << "L1 and L2 grid right !!! " << endl;
    }
    /*不能把所有的二级网格所对应的远场和近场都保存，这样的话会消耗好几百G的储存空间*/
    /*所以改变了策略，那就是将判断远场和近场的部分写成子程序*/
    /*只保存每个二级网格所对应远场的展开系数以及中心点的引力势*/
    // for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex ++ )
    // {
    //     for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex ++ )
    //     {
    //         long* ResNearFieldMicroIndex = NearFieldMicroIndex(SkyLimit_Micro_file, SkyLimit_Micro_file, NStar, X1X2SetL1IncludeMicroIndex, ResolutionL1, X1L1Length, X2L1Length, X1SetL2[X1L2TmpIndex], X2SetL2[X2L2TmpIndex]);
    //     }
    //     cout << X1L2TmpIndex << endl;
    // }
    
    

    /*检验二级网格设置的正确与否，其中包含二级网格在一级网格中的位置是否找准，以及二级网格周围的近场星是否找准
    经验证是准确的*/
    /*
    long TestTestX1IndexL1;
    long TestTestX2IndexL1;
    long TestTestX1IndexL2 = 1000;
    long TestTestX2IndexL2 = 2000;
    
    for(long X1L1TmpIndex = 0; X1L1TmpIndex < X1L1Length; X1L1TmpIndex ++ )
    {
        double LowX1 = X1SetL1[X1L1TmpIndex] - ResolutionL1 / 2;
        double HighX1 = X1SetL1[X1L1TmpIndex] + ResolutionL1 / 2;
        if((X1SetL2[TestTestX1IndexL2] < HighX1) && (X1SetL2[TestTestX1IndexL2] > LowX1))
        {
            TestTestX1IndexL1 = X1L1TmpIndex;
            break;
        }
    }
    for(long X2L1TmpIndex = 0; X2L1TmpIndex < X2L1Length; X2L1TmpIndex ++ )
    {
        double LowX2 = X2SetL1[X2L1TmpIndex] - ResolutionL1 / 2;
        double HighX2 = X2SetL1[X2L1TmpIndex] + ResolutionL1 / 2;
        if((X2SetL2[TestTestX2IndexL2] < HighX2) && (X2SetL2[TestTestX2IndexL2] > LowX2))
        {
            TestTestX2IndexL1 = X2L1TmpIndex;
            break;
        }
    }
    cout << "Test X1L2TmpAtX1SetL1Index " << TestTestX1IndexL1 << endl;
    cout << "Test X2L2TmpAtX2SetL1Index " << TestTestX2IndexL1 << endl;
    for(long MicroTmpI = 0; MicroTmpI < NStar; MicroTmpI ++ )
    {
        if((abs(MicroLensCoorXY[2 * MicroTmpI] - X1SetL1[TestTestX1IndexL1]) <= 3/2. * ResolutionL1) && (abs(MicroLensCoorXY[2 * MicroTmpI + 1] - X2SetL1[TestTestX2IndexL1]) <= 3/2. * ResolutionL1))
        {
            cout << "Test MicroTmpI = " << MicroTmpI << endl;
            // cout << 1/2 * ResolutionL1 << endl;
        }
    }
    
    long* ResNearFieldMicroIndex = NearFieldMicroIndex(SkyLimit_Micro_file, SkyLimit_Micro_file, NStar, X1X2SetL1IncludeMicroIndex, ResolutionL1, X1L1Length, X2L1Length, X1SetL2[TestTestX1IndexL2], X2SetL2[TestTestX2IndexL2]); 
    cout << "Near Micro number" <<ResNearFieldMicroIndex[0] << endl;  
    for(long TestTmpI = 1; TestTmpI <= ResNearFieldMicroIndex[0]; TestTmpI ++ )
    {
        cout << "Near Field Micro Index = " << ResNearFieldMicroIndex[TestTmpI] << endl;
    }
    
    double* ResFarFieldAlphaAndCenterPotential =  FarFieldAlphaAndCenterPotential(MicroLensCoorXY, NStar, ResNearFieldMicroIndex, ResolutionL2, X1SetL2[TestTestX1IndexL2], X2SetL2[TestTestX2IndexL2]);
    for(int TestResFarFieldAlphaAndCenterPotentialIndex = 0; TestResFarFieldAlphaAndCenterPotentialIndex < 9; TestResFarFieldAlphaAndCenterPotentialIndex ++ )
    {
        printf("%.10f\n", ResFarFieldAlphaAndCenterPotential[TestResFarFieldAlphaAndCenterPotentialIndex]); 
    }
    */
    
    // long* ResNearFieldMicroIndex = NearFieldMicroIndex(SkyLimit_Micro_file, SkyLimit_Micro_file, NStar, X1X2SetL1IncludeMicroIndex, ResolutionL1, X1L1Length, X2L1Length, X1SetL2[X1L2Length-1], X2SetL2[X2L2Length-1]);
    // cout << "Test = " << ResNearFieldMicroIndex[0] << endl;

    /*并行插值系数*/
    double*** ResFarFieldAlphaAndCenterPotential = new double** [X1L2Length];
    //与ResNearFieldMicroIndex不同，此变量储存所有二级网格处的近场index, 为的是在后面插值的时候，避免再次计算。
    long*** ResNearFieldMicroIndexSum = new long** [X1L2Length]; 
    for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex ++ )
    {
        ResFarFieldAlphaAndCenterPotential[X1L2TmpIndex] = new double* [X2L2Length];
        ResNearFieldMicroIndexSum[X1L2TmpIndex] = new long* [X2L2Length];
        // for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex ++ )
        // {
        //     ResFarFieldAlphaAndCenterPotential[X1L2TmpIndex][X2L2TmpIndex] = new double [9];
        // }
    }
    std::thread* threads_L2Coeffi[ThreadCount];
    int args_L2Coeffi[ThreadCount];
    long Order_L2Coeffi = X1L2Length/(ThreadCount);
    long Reminder_L2Coeffi = X1L2Length%ThreadCount;
    long RivetingPoint_L2Coeffi_UpLimit;
    long* RivetingPoint_L2Coeffi = new long [ThreadCount];
    if(Reminder_L2Coeffi == 0)
    {
        RivetingPoint_L2Coeffi_UpLimit = Order_L2Coeffi;
    }
    else
    {
        RivetingPoint_L2Coeffi_UpLimit = Order_L2Coeffi + 1;
    }
    
    for(long RivetingPoint_L2Coeffi_Index = 0; RivetingPoint_L2Coeffi_Index < ThreadCount; RivetingPoint_L2Coeffi_Index ++ )
    {
        Reminder_L2Coeffi -= 1;
        if(Reminder_L2Coeffi > 0)
        {
            RivetingPoint_L2Coeffi[RivetingPoint_L2Coeffi_Index] = RivetingPoint_L2Coeffi_UpLimit;
            RivetingPoint_L2Coeffi_UpLimit += Order_L2Coeffi + 1;
        }
        else
        {
            RivetingPoint_L2Coeffi[RivetingPoint_L2Coeffi_Index] = RivetingPoint_L2Coeffi_UpLimit;
            RivetingPoint_L2Coeffi_UpLimit += Order_L2Coeffi;
        }

    }
    cout << "L2Coeffi Runsing One Core Has " << Order_L2Coeffi << " Task" << endl;
    cout << "L2Coeffi Riveting Point [-1] = " << RivetingPoint_L2Coeffi[ThreadCount - 1] << " X1L2Length = " << X1L2Length << endl;

    // for(int tmpi = 0; tmpi < ThreadCount; tmpi ++ )
    // {
    //     long IndexUpLimit = tmpi + 1 == ThreadCount ? X1L2Length : (tmpi + 1) * Order_L2Coeffi; //三目运算符
    //     cout << "LowLimit = " << tmpi*Order_L2Coeffi << " " << "UpLimit = " << IndexUpLimit << endl;
    // }
    auto task_L2Coeffi = [&] (int i) -> void {
        long IndexLowLimit = i == 0 ? 0 : RivetingPoint_L2Coeffi[i-1]; //三目运算符
        for(long X1L2TmpIndex = IndexLowLimit; X1L2TmpIndex < RivetingPoint_L2Coeffi[i]; X1L2TmpIndex++ )
        { 
            double X1L2Tmp = X1SetL2[X1L2TmpIndex]; //算一下L2网格下这一列的x值。
            // cout << "i = " << i << "Index = " << Index << endl;
            for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex ++ )
            {
                double X2L2Tmp = X2SetL2[X2L2TmpIndex]; //算一下坐标的y值。
                // long* ResNearFieldMicroIndex = NearFieldMicroIndex(SkyLimit_Micro_file, SkyLimit_Micro_file, NStar, X1X2SetL1IncludeMicroIndex, ResolutionL1, X1L1Length, X2L1Length, X1L2Tmp, X2L2Tmp);
                ResNearFieldMicroIndexSum[X1L2TmpIndex][X2L2TmpIndex] = NearFieldMicroIndex(SkyLimit_Micro_file, SkyLimit_Micro_file, NStar, X1X2SetL1IncludeMicroIndex, ResolutionL1, X1L1Length, X2L1Length, X1L2Tmp, X2L2Tmp);
                ResFarFieldAlphaAndCenterPotential[X1L2TmpIndex][X2L2TmpIndex] =  FarFieldAlphaAndCenterPotential(MicroLensCoorXY, NStar, ResNearFieldMicroIndexSum[X1L2TmpIndex][X2L2TmpIndex], ResolutionL2, X1L2Tmp, X2L2Tmp, MassSample, AveMassTotal);
            }
        
        }
        
        

    };
    /*并行计算插值系数*/
    auto run_tasks_L2Coeffi = [&] (int count) -> void {
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
                threads_L2Coeffi[CountIndex] = new std::thread (task_L2Coeffi, args_L2Coeffi[CountIndex]);
            }
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
                threads_L2Coeffi[CountIndex]->join();
                delete threads_L2Coeffi[CountIndex]; // 这句话我漏掉了... 你运行的时候记得加上...
            }
        };

    long JThread = 0;
    for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
    {
        args_L2Coeffi[JThread++] = ThreadIndex;
        if (JThread == ThreadCount) {
            run_tasks_L2Coeffi(ThreadCount);
            JThread = 0;
        }
    }
    run_tasks_L2Coeffi(JThread);
    
    // for(int TestPrintResFar = 0; TestPrintResFar < 9; TestPrintResFar ++ )
    // {
    //     cout << "ResFarFieldAlphaAndCenterPotential[0][0]" << ResFarFieldAlphaAndCenterPotential[0][0][TestPrintResFar] << endl;
    // }
    double CostTimeL2Coeffi = time(NULL) - TimeStart;
    cout << "CostTimeL2Coeffi time : " << CostTimeL2Coeffi  << " s" << endl;
    /*以上*/

    // //测试引力势和偏转角是否正确。
    // // long* ResNearFieldMicroIndex = NearFieldMicroIndex(SkyLimit_Micro_file, SkyLimit_Micro_file, NStar, X1X2SetL1IncludeMicroIndex, ResolutionL1, X1L1Length, X2L1Length, X1SetL2[X1L2Length - 2], X2SetL2[X2L2Length - 2]);
    // double TestPsi = MicroCreatPsi(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, X1SetL2[X1L2Length - 2] - 0.1, X2SetL2[X2L2Length - 2] - 0.5, MicroLensCoorXY, NStar, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2);
    // double* TestAlpha = MacroMicroPartialTau(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, kappaStar, X1SetL2[X1L2Length - 2] - 0.1, X2SetL2[X2L2Length - 2] - 0.5, MicroLensCoorXY, NStar, - SkyLimit_Micro_file, SkyLimit_Micro_file, - SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2);
    // printf("%.10f\n", TestPsi);
    // printf("%.10f\n", TestAlpha[0]);
    // printf("%.10f\n", TestAlpha[1]);
    // printf("%.10f\n", X1SetL2[X1L2Length - 2] - 0.1);
    // printf("%.10f\n", X2SetL2[X2L2Length - 2] - 0.5);
    // double TestPsiAndMinusSheet = MacroMicroAndMinusSheet(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, X1SetL2[X1L2Length - 2] - 0.1, X2SetL2[X2L2Length - 2] - 0.5, kappaStar, MicroLensCoorXY, NStar, - SkyLimit_Micro_file, SkyLimit_Micro_file, - SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2);
    // printf("%.10f\n", TestPsiAndMinusSheet);
    // int LayersNumber_tmp = MacroMicroLayersNumber(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, kappaStar, X1SetL2[0], X2SetL2[X2L2Length - 1], ResolutionL2, EpsilonRes, TimeDelayEsp, coeffi, MicroLensCoorXY, NStar, -SkyLimit_Micro_file, SkyLimit_Micro_file, -SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2);
    /*保存微透镜场，画图看一下结果*/
    cout << ">>>>>>>>>>>>>>>Saving Microlensing Field<<<<<<<<<<<<<<"<< endl;
    char MicroFieldName[100];
    sprintf(MicroFieldName,"./MicroField_%d/MicroField_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift); 
    ofstream MicroFieldSave;
    MicroFieldSave.open(MicroFieldName, std::ofstream::binary);

    double** Level2GridMicroPotential = new double* [X1L2Length];
    for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex ++ )
    {
        Level2GridMicroPotential[X1L2TmpIndex] = new double [X2L2Length];
    }

    for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex++ )
    { 
        double X1L2Tmp = X1SetL2[X1L2TmpIndex]; //算一下L2网格下这一列的x值。
        // cout << "i = " << i << "Index = " << Index << endl;
        for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex ++ )
        {
            double X2L2Tmp = X2SetL2[X2L2TmpIndex]; //算一下坐标的y值。
            Level2GridMicroPotential[X1L2TmpIndex][X2L2TmpIndex] =  MacroMicroAndMinusSheet(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, X1L2Tmp, X2L2Tmp, kappaStar, MicroLensCoorXY, NStar, - SkyLimit_Micro_file, SkyLimit_Micro_file, - SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal) - 0.5 * ((1 - kappa + gamma) * X1L2Tmp * X1L2Tmp + (1 - kappa - gamma) * X2L2Tmp * X2L2Tmp);
        }
    }

    for(long i = 0; i < X1L2Length; i ++)
    {
        MicroFieldSave.write(reinterpret_cast<char *>(Level2GridMicroPotential[i]), 8*X2L2Length);
    }
    cout << ">>>>>>>>>>>Saved done<<<<<<<<<<<<<" << endl;
    MicroFieldSave.close();

    //所有文件
    //储存时间序列长度
    vector <long> TimeLength;
    char TimeLengthFile[100];
    //保存加密网格验证一下。
    vector <long> LayersNumber_Res;
    char LayersFile[100];
    //Area和Time的文件名
    char AreaName[100];
    char TimeName[100];


    ofstream TLFile;
    ofstream LF;
    ofstream fp_area;
    ofstream fp_time;

    //保存平均质量，总微透镜数目，L2网格长度（用来画天图）。
    vector <double> AveMassTotalAndTotalNum;
    AveMassTotalAndTotalNum.push_back(AveMassTotal);
    AveMassTotalAndTotalNum.push_back(NStar);
    AveMassTotalAndTotalNum.push_back(X1L2Length);
    char AveMassAndNumName[100];
    sprintf(AveMassAndNumName,"./MicroField_%2d/AveMassAndNum_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift); 
    ofstream AveMassAndNumFile;
    AveMassAndNumFile.open(AveMassAndNumName, std::ofstream::binary);
    AveMassAndNumFile.write(reinterpret_cast<char *>(&AveMassTotalAndTotalNum[0]), sizeof(double)*3);
    AveMassAndNumFile.close();



    
    std::thread* threads[ThreadCount];
    int args[ThreadCount];

    //创建二维的数组用来存放LayerNum，并初始化为0；
    vector <long> LayersNumber_Res2D[ThreadCount];

    /*并行计算加密层数*/
    //注意，下面的计算没有包括最右边和最上面的那一行，因为在后面自适应部分，会碰到这个边界，然后会被认成在下一级的二级网格中，所以会报错。
    //但是在粗略计算时间延迟的时候，依然包含了这一行一列，因为这个计算只会在二级网格中心计算，不会细分，所以不会碰到边界。
    //这个处理与主程序中处理方法一样。
    long LayerNumSum = 0; 
    std::thread* threads_LayerNum[ThreadCount];
    int args_LayerNum[ThreadCount];
    long Order_LayerNum = (X1L2Length - 1)/(ThreadCount); 

    long Reminder_LayerNum = (X1L2Length - 1)%ThreadCount;
    long RivetingPoint_LayerNum_UpLimit;
    long* RivetingPoint_LayerNum = new long [ThreadCount];
    if(Reminder_LayerNum == 0)
    {
        RivetingPoint_LayerNum_UpLimit = Order_LayerNum;
    }
    else
    {
        RivetingPoint_LayerNum_UpLimit = Order_LayerNum + 1;
    }
    
    for(long RivetingPoint_LayerNum_Index = 0; RivetingPoint_LayerNum_Index < ThreadCount; RivetingPoint_LayerNum_Index ++ )
    {
        Reminder_LayerNum -= 1;
        if(Reminder_LayerNum > 0)
        {
            RivetingPoint_LayerNum[RivetingPoint_LayerNum_Index] = RivetingPoint_LayerNum_UpLimit;
            RivetingPoint_LayerNum_UpLimit += Order_LayerNum + 1;
        }
        else
        {
            RivetingPoint_LayerNum[RivetingPoint_LayerNum_Index] = RivetingPoint_LayerNum_UpLimit;
            RivetingPoint_LayerNum_UpLimit += Order_LayerNum;
        }

    }
    
    cout << "LayerNum Runsing One Core Has " << Order_LayerNum << " Task" << endl;
    cout << "LayerNum Riveting Point [-1] = " << RivetingPoint_LayerNum[ThreadCount - 1] << " X1L2Length - 1 = " << X1L2Length - 1 << endl;

    auto task_LayerNum = [&] (int i) -> void {
        long IndexLowLimit = i == 0 ? 0 : RivetingPoint_LayerNum[i - 1]; //三目运算符
        for(long X1L2TmpIndex = IndexLowLimit; X1L2TmpIndex < RivetingPoint_LayerNum[i]; X1L2TmpIndex++ )
        { 
            double X1L2Tmp = X1SetL2[X1L2TmpIndex]; //算一下L2网格下这一列的x值。
            // cout << "i = " << i << "Index = " << Index << endl;
            for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length - 1; X2L2TmpIndex ++ )
            {
                double X2L2Tmp = X2SetL2[X2L2TmpIndex]; //算一下坐标的y值。
                int LayersNumber_tmp = MacroMicroLayersNumber(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, kappaStar, X1L2Tmp, X2L2Tmp, ResolutionL2, EpsilonRes, TimeDelayEsp, coeffi, MicroLensCoorXY, NStar, -SkyLimit_Micro_file, SkyLimit_Micro_file, -SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
                // cout << "LayersNumber_tmp = " << LayersNumber_tmp << endl;
                LayersNumber_Res2D[i].push_back(LayersNumber_tmp);
            }
        
        }
        
        

    };
    /*并行计算网格加密层数*/
    auto run_tasks_LayerNum = [&] (int count) -> void {
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
                threads_LayerNum[CountIndex] = new std::thread (task_LayerNum, args_LayerNum[CountIndex]);
            }
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
                threads_LayerNum[CountIndex]->join();
                delete threads_LayerNum[CountIndex]; // 这句话我漏掉了... 你运行的时候记得加上...
            }
        };

    long JThread_LayerNum = 0;
    for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
    {
        args_LayerNum[JThread_LayerNum++] = ThreadIndex;
        if (JThread_LayerNum == ThreadCount) {
            run_tasks_LayerNum(ThreadCount);
            JThread_LayerNum = 0;
        }
    }
    run_tasks_LayerNum(JThread_LayerNum);
    
    double CostTimeLayerNum = time(NULL) - TimeStart;
    cout << "CostTimeLayerNum time : " << CostTimeLayerNum - CostTimeL2Coeffi  << " s" << endl;
    //计算总加密格点数
    long TotalLengthLayerNum_Res = 0;
    for(long LayerNumIndex = 0; LayerNumIndex < ThreadCount; LayerNumIndex ++)
    {
        long LengthTmp2D = LayersNumber_Res2D[LayerNumIndex].size();
        TotalLengthLayerNum_Res += LengthTmp2D;
        for(long IndexTmp = 0; IndexTmp < LengthTmp2D; IndexTmp ++)
        {
            LayerNumSum += pow(4, LayersNumber_Res2D[LayerNumIndex][IndexTmp]);
            LayersNumber_Res.push_back(LayersNumber_Res2D[LayerNumIndex][IndexTmp]);
        }
        
    }
    cout << "TotalLengthLayerNum_Res = " << TotalLengthLayerNum_Res << " (LengthX - 1)^2 = " << (X1L2Length - 1) * (X2L2Length - 1) << endl; 
    cout << "Total Layer Num = " << LayerNumSum << endl;
    //重新定义每个核心计算的加密以后的网格数目
    long Order = LayerNumSum/(ThreadCount); //这里选取线程数的时候，最好可以整除，要不然会变慢很多。
    cout << "One Core Has LayerNum " << Order << " Task" << endl;
    //找到加密网格对应的二级网格位置
    vector <double> X1SetL2OneDim; //完整的一维坐标，不包括右边和上边边界
    vector <double> X2SetL2OneDim;
    for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length - 1; X1L2TmpIndex ++ )
    {
        for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length - 1; X2L2TmpIndex ++ )
        { 
            X1SetL2OneDim.push_back(X1SetL2[X1L2TmpIndex]);
            X2SetL2OneDim.push_back(X2SetL2[X2L2TmpIndex]);
        }
    }
    cout << "X1SetL2OneDim size = " << X1SetL2OneDim.size() << " (X1L2Length - 1) * (X2L2Length - 1)" << (X1L2Length - 1) * (X2L2Length - 1) << endl;
    long* RivetingPoint = new long [ThreadCount];
    long LayerNumSumTmp = 0;
    int RivetingPoint_i = 0;
    for(int Layer_i = 0; Layer_i < TotalLengthLayerNum_Res; Layer_i ++ )
    {
        if(LayerNumSumTmp > (RivetingPoint_i + 1) * Order)
        {
            RivetingPoint[RivetingPoint_i] = Layer_i;
            RivetingPoint_i += 1;
        }
        LayerNumSumTmp += pow(4, LayersNumber_Res[Layer_i]);

    }
    if(RivetingPoint[ThreadCount - 2] == 0)
    {
        RivetingPoint[ThreadCount - 2] = TotalLengthLayerNum_Res - 1; 
    }
    RivetingPoint[ThreadCount - 1] = TotalLengthLayerNum_Res - 1;
    // for(int Tmp_i = 0; Tmp_i < ThreadCount; Tmp_i ++ )
    // {
    //     cout << Tmp_i << " : " << RivetingPoint[Tmp_i] << endl;
    // }

    
    
    
    /*下面是粗略计算最大值点和最小值点的子程序*/
    double* RudeTimeDelayMinAndMax = new double[2]; //这个是用来存后面得到的最大值和最小值的。
    RudeTimeDelayMinAndMax[0] = 100000000000;
    RudeTimeDelayMinAndMax[1] = -10000000000;
    double TestX1, TestX2;
    long TestX1Index, TestX2Index;

    //Minimum
    if((mur > 0)&&(mut > 0))
    {
        sprintf(LayersFile,"ResultMinimum_%2d/LayersFile_min_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        sprintf(TimeLengthFile,"ResultMinimum_%2d/TimeLength_min_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        LF.open(LayersFile, std::ofstream::binary);
        TLFile.open(TimeLengthFile, std::ofstream::binary);

        double Axisa = sqrt(1/coeffi/mur);
        double Axisb = sqrt(1/coeffi/mut);
        
        sprintf(AreaName,"ResultMinimum_%2d/adptive_Area_min_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        sprintf(TimeName,"ResultMinimum_%2d/adptive_Time_min_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift); 
        
        fp_area.open(AreaName, std::ofstream::binary);
        fp_time.open(TimeName, std::ofstream::binary);
        


        for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex += 1 )
        {
            double X1L2Tmp = X1SetL2[X1L2TmpIndex];

            // cout << RudeX1Tmp << endl;
            for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex += 1 )
            {
                double X2L2Tmp = X2SetL2[X2L2TmpIndex];
                double RudePsi = MacroMicroAndMinusSheet(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, X1L2Tmp, X2L2Tmp, kappaStar, MicroLensCoorXY, NStar, - SkyLimit_Micro_file, SkyLimit_Micro_file, - SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
                // cout << "RudeX1Tmp = " << RudeX1Tmp << endl;
                // cout << "RudeX2Tmp = " << RudeX2Tmp << endl;
                // cout << "coeffi = " << coeffi << endl;
                // cout << "RudePsi = " << RudePsi << endl;
                double MacroT_tmp = pow(X1L2Tmp,2)/(2*Axisa*Axisa) + pow(X2L2Tmp,2)/(2*Axisb*Axisb);
                if((MacroT_tmp > -TimeNeed2CalMax)&&(MacroT_tmp < TimeNeed2CalMax))
                {
                    if(RudePsi > RudeTimeDelayMinAndMax[1])
                    {
                        RudeTimeDelayMinAndMax[1] = RudePsi;
                        
                    }
                    if(RudePsi < RudeTimeDelayMinAndMax[0])
                    {
                        RudeTimeDelayMinAndMax[0] = RudePsi;
                        TestX1 = X1L2Tmp;
                        TestX2 = X2L2Tmp;
                        TestX1Index = X1L2TmpIndex;
                        TestX2Index = X2L2TmpIndex;
                        // cout << "Test = " << X1Set[TestX1Index] << endl;
                        
                    }
                }
            }
            
        }
        double RudeTimeDelayMin = RudeTimeDelayMinAndMax[0] * coeffi - 0.1;
        double RudeTimeDelayMax = RudeTimeDelayMinAndMax[1] * coeffi;
        cout << "Rude Minimum Time Delay is " << RudeTimeDelayMin << endl;
        cout << "Rude Minimum X1 = " << TestX1 << " X2 = " << TestX2 << endl;
        cout << "Rude Minimum X1 Index = " << TestX1Index << " X2 Index = " << TestX2Index << endl; 
        cout << "Rude Maximun Time Delay is " << RudeTimeDelayMax << endl;
        delete[] RudeTimeDelayMinAndMax;
        /*以上*/


        /*下面并行计算时间延迟曲线*/

        //the next is for time range
        double TimeDelayStart = RudeTimeDelayMin; 
        
        double TimeDelayEnd = RudeTimeDelayMax;
        //
        cout << "Time Maximum Boundary For Time Delay Curve " << TimeDelayEnd << endl;
        long LenTimeDelay = (TimeDelayEnd - TimeDelayStart + TimeDelayEsp)/TimeDelayEsp; 
        // the length of time delay area
        cout << "Time Delay Curve's Length = " << LenTimeDelay << endl; //print len_time
        TimeLength.push_back(LenTimeDelay);
        
        
        double* TimeDelayRange = new double [LenTimeDelay];
        double TimeDelayTmp = TimeDelayStart; //the start of the time delay calculation.
        //generate time array.
        for(long TimeDelayIndex = 0; TimeDelayIndex < LenTimeDelay; TimeDelayIndex ++ )
        {
        
            TimeDelayRange[TimeDelayIndex] = TimeDelayTmp;
            TimeDelayTmp += TimeDelayEsp;
            
            

        }
        cout << TimeDelayRange[LenTimeDelay - 1] << endl;
        //cout << int (sizeof(time_range)/sizeof(time_range[0])) << endl;


        double* TimeDelayArea = new double [LenTimeDelay](); //存放时间延迟面积的数组。
        //创建二维的数组用来存放时间延迟面积的数据，并初始化为0；
        double** TimeDelayArea2D = new double*[ThreadCount];
        for(int i = 0; i < ThreadCount; i++)
        {
            TimeDelayArea2D[i] = new double[LenTimeDelay]();
        }
        

        
        auto task = [&] (int i) -> void {

        
            long IndexLowLimit = i == 0 ? 0 : RivetingPoint[i - 1]; //三目运算符
            for(long X1X2L2TmpIndexOneDim = IndexLowLimit; X1X2L2TmpIndexOneDim <= RivetingPoint[i]; X1X2L2TmpIndexOneDim++ )
            { 
                double X1L2Tmp = X1SetL2OneDim[X1X2L2TmpIndexOneDim]; //算一下这一列的x值。
                double X2L2Tmp = X2SetL2OneDim[X1X2L2TmpIndexOneDim]; //算一下坐标的y值。
                // cout << "X2L2Tmp = " << X2L2Tmp << endl;
                double MacroT_tmp = pow(X1L2Tmp,2)/(2*Axisa*Axisa) + pow(X2L2Tmp,2)/(2*Axisb*Axisb);
                int LayersNumber_tmp = LayersNumber_Res[X1X2L2TmpIndexOneDim];
                
                if((MacroT_tmp > -TimeNeed2CalMax - TimeNeed2CalMax / 10)&&(MacroT_tmp < TimeNeed2CalMax + TimeNeed2CalMax / 10)) //加上了一个保护区
                //判断坐标是不是在双曲线区域内。
                {
                    double Resolution = ResolutionL2 / pow(2, LayersNumber_tmp);
                    // cout << "Resolution = " << Resolution << endl;
                
                    double BoostImageRes = Resolution;
                        
                    int Boost = int(ResolutionL2 / Resolution);
                    
                    /*加密以后的天区坐标*/
                    double* X1SetBoost = new double [Boost];
                    double* X2SetBoost = new double [Boost];
                    int X1BoostIndex = 0;
                    int X2BoostIndex = 0;
                    for(double X1SetBoostTmp = X1L2Tmp - ResolutionL2/2 + BoostImageRes/2; X1SetBoostTmp < X1L2Tmp + ResolutionL2/2; X1SetBoostTmp += BoostImageRes)
                    {
                        if(X1BoostIndex >= Boost)
                        {
                            
                            break;
                        }
                        else
                        {
                            X1SetBoost[X1BoostIndex] = X1SetBoostTmp;
                            X1BoostIndex += 1;
                        }
                    }
                    
                    for(double X2SetBoostTmp = X2L2Tmp - ResolutionL2/2 + BoostImageRes/2; X2SetBoostTmp < X2L2Tmp + ResolutionL2/2; X2SetBoostTmp += BoostImageRes)
                    {
                        if(X2BoostIndex >= Boost)
                        {
                            
                            break;
                        }
                        else
                        {
                            X2SetBoost[X2BoostIndex] = X2SetBoostTmp;
                            X2BoostIndex += 1;
                        }

                    }
                    

                    for(long X1Indextmp = 0; X1Indextmp < Boost; X1Indextmp++ )
                    { 
                        double TmpTmpX1 = X1SetBoost[X1Indextmp]; //算一下这一列的x值。
                        
                        for(long X2Indextmp = 0; X2Indextmp < Boost; X2Indextmp += 1)
                        {
                            double TmpTmpX2 = X2SetBoost[X2Indextmp]; //算一下坐标的y值。
                            double MacroT_tmp_tmp = pow(TmpTmpX1,2)/(2*Axisa*Axisa) + pow(TmpTmpX2,2)/(2*Axisb*Axisb);
                            if((MacroT_tmp_tmp > -TimeNeed2CalMax)&&(MacroT_tmp_tmp < TimeNeed2CalMax))
                            {   
                                double* TSUMTmp = MacroMicroTheoryOut(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, kappaStar, TmpTmpX1, TmpTmpX2, BoostImageRes, coeffi, MicroLensCoorXY, NStar, -SkyLimit_Micro_file, SkyLimit_Micro_file, -SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
                                long TimeDelayIndexStart = long ((TSUMTmp[2] - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) - 10;
                                long TimeDelayIndexEnd = long ((TSUMTmp[3] - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) + 10;
                                for(long TimeDelayIndexTmp = TimeDelayIndexStart; TimeDelayIndexTmp <= TimeDelayIndexEnd; TimeDelayIndexTmp ++ )
                                {
                                    if((TimeDelayIndexTmp >=0)&&(TimeDelayIndexTmp < LenTimeDelay))
                                    {
                                        double TimeDelayValue = TimeDelayRange[TimeDelayIndexTmp];
                                        double TimeDelayLowerLimit = TimeDelayValue - TimeDelayEsp;
                                        double AreaTmp = 0;

                                        if((TimeDelayValue >= TSUMTmp[2])&&(TimeDelayValue <= TSUMTmp[0]))
                                        {
                                            if(TimeDelayLowerLimit < TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * pow(TimeDelayValue - TSUMTmp[2], 2) / (TSUMTmp[0] - TSUMTmp[2]) / 2;
                                            }
                                            else
                                            {
                                                AreaTmp = TSUMTmp[4] * (2 * TimeDelayValue - 2 * TSUMTmp[2] - TimeDelayEsp) / (TSUMTmp[0] - TSUMTmp[2]) * TimeDelayEsp / 2;
                                            }

                                        } 
                                        else if((TimeDelayValue > TSUMTmp[0])&&(TimeDelayValue <= TSUMTmp[1]))
                                        {
                                            if(TimeDelayLowerLimit <= TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                + TSUMTmp[4] * (TimeDelayValue - TSUMTmp[0]);
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2 
                                                + TSUMTmp[4] * (TimeDelayValue - TSUMTmp[0]); 
                                            }
                                            else
                                            {
                                                AreaTmp = TSUMTmp[4] * TimeDelayEsp;
                                            }
                                        }
                                        else if((TimeDelayValue > TSUMTmp[1])&&(TimeDelayValue <= TSUMTmp[3]))
                                        {
                                            if(TimeDelayLowerLimit <= TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) 
                                                + (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2;
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2
                                                    + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) + 
                                                    (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2; 
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[1])&&(TimeDelayLowerLimit > TSUMTmp[0]))
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[1] - TimeDelayLowerLimit) + 
                                                    (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2;
                                            } 
                                            else
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TSUMTmp[3] - TimeDelayLowerLimit) / (TSUMTmp[3] - TSUMTmp[1]) + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * TimeDelayEsp / 2; 
                                            }
                                        }
                                        else if(TimeDelayValue > TSUMTmp[3])
                                        {
                                            if(TimeDelayLowerLimit <= TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) 
                                                + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2;
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2
                                                    + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0])
                                                    + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2; 
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[1])&&(TimeDelayLowerLimit > TSUMTmp[0]))
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[1] - TimeDelayLowerLimit)
                                                    + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2;
                                            } 
                                            else if((TimeDelayLowerLimit <= TSUMTmp[3])&&(TimeDelayLowerLimit > TSUMTmp[1]))
                                            {
                                                AreaTmp = TSUMTmp[4] * pow(TSUMTmp[3] - TimeDelayLowerLimit, 2) / (TSUMTmp[3] - TSUMTmp[1]) / 2; 
                                            }
                                            else
                                            {
                                                AreaTmp = 0;
                                            }
                                        }
                                        
                                        TimeDelayArea2D[i][TimeDelayIndexTmp] += AreaTmp; //FIXME 这里可能有错误，比如说多个线程同时相加。
                                    }
                                }       
                                delete[] TSUMTmp;   
                            }
                        }
                    }
                    delete[] X1SetBoost;
                    delete[] X2SetBoost;


                
                }
            
            }
            
            

        };

        auto run_tasks = [&] (int count) -> void {
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
                threads[CountIndex] = new std::thread (task, args[CountIndex]);
            }
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
                threads[CountIndex]->join();
                delete threads[CountIndex]; // 这句话我漏掉了... 你运行的时候记得加上...
            }
        };

        long JThread = 0;
        for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
        {
            args[JThread++] = ThreadIndex;
            if (JThread == ThreadCount) {
                run_tasks(ThreadCount);
                JThread = 0;
            }
        }
        run_tasks(JThread);
        
        for(long AreaIndex = 0; AreaIndex < LenTimeDelay; AreaIndex ++)
        {
            for(long AreaSumIndex = 0; AreaSumIndex < ThreadCount; AreaSumIndex ++)
            {
                TimeDelayArea[AreaIndex] += TimeDelayArea2D[AreaSumIndex][AreaIndex];
            }
        }

        fp_area.write(reinterpret_cast<char *>(TimeDelayArea), sizeof(double)*LenTimeDelay);
        fp_time.write(reinterpret_cast<char *>(TimeDelayRange), sizeof(double)*LenTimeDelay);
        
        

        delete[] TimeDelayRange;
        delete[] TimeDelayArea;
        delete[] TimeDelayArea2D;
    }

    //Saddle
    if(mur * mut < 0)
    {
        sprintf(LayersFile,"ResultSaddle_%2d/LayersFile_sad_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift); 
        sprintf(TimeLengthFile,"ResultSaddle_%2d/TimeLength_sad_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        LF.open(LayersFile, std::ofstream::binary); 
        TLFile.open(TimeLengthFile, std::ofstream::binary);

        mur = 1 - kappa + gamma;
        mut = kappa + gamma - 1; //FIXME注意，这里有可能反号，好处是在SIE模型中kappa = gamma，所以mur = 1, mut<0
        double Axisa = sqrt(1/coeffi/mur);
        double Axisb = sqrt(1/coeffi/mut);
        cout << "Axisa = " << Axisa << endl;
        cout << "Axisb = " << Axisb << endl;
        double Epsilon1 = PreOutPut[8];// 38.355;//30;
        double Epsilon2 = PreOutPut[9];// Axisb/Axisa*Epsilon1;;
        double X10New = PreOutPut[6];// Axisb*(-pow(Epsilon1,2) + 2 * pow(Axisa,2) * TimeNeed2CalMax)/(2 * Axisa * Epsilon1);
        double X20New = PreOutPut[7];// Axisa * Epsilon2 / 2 / Axisb + Axisa * Axisb * TimeNeed2CalMax / Epsilon2;
        cout << "X10New = " << X10New << " " << "X20New = " << X20New << endl;

        //保存X10New和X20New。
        vector <double> X1020New;
        X1020New.push_back(X10New);
        X1020New.push_back(X20New);
        char X1020NewName[100];
        sprintf(X1020NewName,"./ResultSaddle_%2d/X1020New_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift); 
        ofstream X1020NewFile;
        X1020NewFile.open(X1020NewName, std::ofstream::binary);
        X1020NewFile.write(reinterpret_cast<char *>(&X1020New[0]), sizeof(double)*2);
        X1020NewFile.close();
        //以上
        
        sprintf(AreaName,"ResultSaddle_%2d/adptive_Area_sad_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        sprintf(TimeName,"ResultSaddle_%2d/adptive_Time_sad_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift); 
        
        fp_area.open(AreaName, std::ofstream::binary);
        fp_time.open(TimeName, std::ofstream::binary);
        

        for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex += 1 )
        {
            double X1L2Tmp = X1SetL2[X1L2TmpIndex];
            // cout << RudeX1Tmp << endl;
            for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex += 1 )
            {
                double X2L2Tmp = X2SetL2[X2L2TmpIndex];
                double RudePsi = MacroMicroAndMinusSheet(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, X1L2Tmp, X2L2Tmp, kappaStar, MicroLensCoorXY, NStar, - SkyLimit_Micro_file, SkyLimit_Micro_file, - SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
                double MacroT_tmp = pow(X1L2Tmp,2)/(2*Axisa*Axisa) - pow(X2L2Tmp,2)/(2*Axisb*Axisb);
                if((MacroT_tmp > -TimeNeed2CalMax)&&(MacroT_tmp < TimeNeed2CalMax))
                {
                    if(((MacroT_tmp >= 0)&(X2L2Tmp > - X20New)&&(X2L2Tmp < X20New)) || ((MacroT_tmp < 0)&&(X1L2Tmp > - X10New)&&(X1L2Tmp < X10New)))
                    {
                        double RudeTimeDelayTmp = RudePsi;
                        if(RudeTimeDelayTmp > RudeTimeDelayMinAndMax[1])
                        {
                            RudeTimeDelayMinAndMax[1] = RudeTimeDelayTmp;
                        }
                        if(RudePsi < RudeTimeDelayMinAndMax[0])
                        {
                            RudeTimeDelayMinAndMax[0] = RudePsi;
                            TestX1 = X1SetL2[X1L2TmpIndex];
                            TestX2 = X2SetL2[X2L2TmpIndex];
                            TestX1Index = X1L2TmpIndex;
                            TestX2Index = X2L2TmpIndex;
                            // cout << "Test = " << X1Set[TestX1Index] << endl;
                            
                        }
                    }
                    
                }
            }
            
        }
        double RudeTimeDelayMin = RudeTimeDelayMinAndMax[0] * coeffi;
        double RudeTimeDelayMax = RudeTimeDelayMinAndMax[1] * coeffi;
        cout << "Rude Minimum Time Delay is " << RudeTimeDelayMin << endl;
        cout << "Rude Minimum X1 = " << TestX1 << " X2 = " << TestX2 << endl;
        cout << "Rude Minimum X1 Index = " << TestX1Index << " X2 Index = " << TestX2Index << endl; 
        cout << "Rude Maximun Time Delay is " << RudeTimeDelayMax << endl;
        delete[] RudeTimeDelayMinAndMax;
        /*以上*/


        /*下面并行计算时间延迟曲线*/

        //the next is for time range
        double TimeDelayStart = RudeTimeDelayMin; 
        
        double TimeDelayEnd = RudeTimeDelayMax;
        //
        cout << "Time Maximum Boundary For Time Delay Curve " << TimeDelayEnd << endl;
        long LenTimeDelay = (TimeDelayEnd - TimeDelayStart + TimeDelayEsp)/TimeDelayEsp; 
        // the length of time delay area
        cout << "Time Delay Curve's Length = " << LenTimeDelay << endl; //print len_time
        TimeLength.push_back(LenTimeDelay);
        
        
        double* TimeDelayRange = new double [LenTimeDelay];
        double TimeDelayTmp = TimeDelayStart; //the start of the time delay calculation.
        //generate time array.
        for(long TimeDelayIndex = 0; TimeDelayIndex < LenTimeDelay; TimeDelayIndex ++ )
        {
        
            TimeDelayRange[TimeDelayIndex] = TimeDelayTmp;
            TimeDelayTmp += TimeDelayEsp;
            
            

        }
        cout << TimeDelayRange[LenTimeDelay - 1] << endl;
        //cout << int (sizeof(time_range)/sizeof(time_range[0])) << endl;


        double* TimeDelayArea = new double [LenTimeDelay](); //存放时间延迟面积的数组。



        
        //创建二维的数组用来存放时间延迟面积的数据，并初始化为0；
        double** TimeDelayArea2D = new double*[ThreadCount];
        for(int i = 0; i < ThreadCount; i++)
        {
            TimeDelayArea2D[i] = new double[LenTimeDelay]();
        }
        

        auto task = [&] (int i) -> void {

          
        
            long IndexLowLimit = i == 0 ? 0 : RivetingPoint[i - 1]; //三目运算符
            for(long X1X2L2TmpIndexOneDim = IndexLowLimit; X1X2L2TmpIndexOneDim <= RivetingPoint[i]; X1X2L2TmpIndexOneDim++ )
            { 
                double X1L2Tmp = X1SetL2OneDim[X1X2L2TmpIndexOneDim]; //算一下这一列的x值。
                double X2L2Tmp = X2SetL2OneDim[X1X2L2TmpIndexOneDim]; //算一下坐标的y值。
                double MacroT_tmp = pow(X1L2Tmp,2)/(2*Axisa*Axisa) - pow(X2L2Tmp,2)/(2*Axisb*Axisb);
                int LayersNumber_tmp = LayersNumber_Res[X1X2L2TmpIndexOneDim];
                if((MacroT_tmp > -TimeNeed2CalMax - TimeNeed2CalMax / 10)&&(MacroT_tmp < TimeNeed2CalMax + TimeNeed2CalMax / 10)) //加上了一个保护区
                //判断是不是坐标是不是在双曲线区域内。
                {
                    // double GradsMacro = sqrt(pow(X1L2Tmp, 2) / pow(Axisa, 4) + pow(X2L2Tmp, 2) / pow(Axisb, 4));
                    // double AlphaMacro = acos((X1L2Tmp / pow(Axisa, 2) - X2L2Tmp / pow(Axisb, 2))/(sqrt(2)*GradsMacro));
                    // if((AlphaMacro >= 0)&(AlphaMacro < PI/4))
                    // {
                    //     AlphaMacro = AlphaMacro;
                    // }
                    // else if((AlphaMacro < PI/2)&&(AlphaMacro >= PI/4))
                    // {
                    //     AlphaMacro = PI/2 - AlphaMacro;
                    // }
                    // else if((AlphaMacro >= PI/2)&&(AlphaMacro < 3*PI/4))
                    // {
                    //     AlphaMacro = AlphaMacro - PI/2;
                    // }
                    // else if((AlphaMacro >= 3*PI/4)&&(AlphaMacro < PI))
                    // {
                    //     AlphaMacro = PI - AlphaMacro;
                    // }
                    // double L1Macro = sqrt(2) * ResolutionL2 * cos(AlphaMacro);
                    // double MaxDeltaTInPixel = L1Macro * GradsMacro;
                    // // cout << "MaxDeltaTInPixel = " << MaxDeltaTInPixel << endl;
                    // if((MacroT_tmp >=0)&&(pow(MacroT_tmp - TimeNeed2CalMax,2) < pow(MaxDeltaTInPixel,2)) || (MacroT_tmp <=0)&&(pow(MacroT_tmp + TimeNeed2CalMax, 2) < pow(MaxDeltaTInPixel, 2)))
                    // {
                    //     LayersNumber_tmp = 7;
                    //     LayersNumber_Res[X1X2L2TmpIndexOneDim] = LayersNumber_tmp;
                    //     // cout << "MaxDeltaTInPixel = " << MaxDeltaTInPixel << endl;
                    //     // cout << "X1L2Tmp = " << X1L2Tmp << endl;
                    //     // cout << "X2L2Tmp = " << X2L2Tmp << endl;
                    // }

                    if((MacroT_tmp >=0)&&(X2L2Tmp > - X20New - X20New / 10)&&(X2L2Tmp < X20New + X20New / 10) || (MacroT_tmp <0)&&(X1L2Tmp > - X10New - X10New / 10)&&(X1L2Tmp < X10New + X10New / 10)) //加了一个保护区
                    {
                        double Resolution = ResolutionL2 / pow(2, LayersNumber_tmp);
                        // cout << "Resolution = " << Resolution << endl;
                        double BoostImageRes = Resolution;
                        
                        int Boost = int(ResolutionL2 / Resolution);
                        

                        /*加密以后的天区坐标*/
                        double* X1SetBoost = new double [Boost];
                        double* X2SetBoost = new double [Boost];
                        int X1BoostIndex = 0;
                        int X2BoostIndex = 0;
                        for(double X1SetBoostTmp = X1L2Tmp - ResolutionL2/2 + BoostImageRes/2; X1SetBoostTmp < X1L2Tmp + ResolutionL2/2; X1SetBoostTmp += BoostImageRes)
                        {
                            if(X1BoostIndex >= Boost)
                            {
                                
                                break;
                            }
                            else
                            {
                                X1SetBoost[X1BoostIndex] = X1SetBoostTmp;
                                X1BoostIndex += 1;
                            }
                        }

                        
                        
                        
                        for(double X2SetBoostTmp = X2L2Tmp - ResolutionL2/2 + BoostImageRes/2; X2SetBoostTmp < X2L2Tmp + ResolutionL2/2; X2SetBoostTmp += BoostImageRes)
                        {
                            if(X2BoostIndex >= Boost)
                            {
                                
                                break;
                            }
                            else
                            {
                                X2SetBoost[X2BoostIndex] = X2SetBoostTmp;
                                X2BoostIndex += 1;
                            }

                        }
        
                        for(long X1Indextmp = 0; X1Indextmp < Boost; X1Indextmp++ )
                        { 
                            double TmpTmpX1 = X1SetBoost[X1Indextmp]; //算一下这一列的x值。
                            
                            for(long X2Indextmp = 0; X2Indextmp < Boost; X2Indextmp += 1)
                            {
                                double TmpTmpX2 = X2SetBoost[X2Indextmp]; //算一下坐标的y值。
                                double MacroT_tmp_tmp = pow(TmpTmpX1,2)/(2*Axisa*Axisa) - pow(TmpTmpX2,2)/(2*Axisb*Axisb);
                                if((MacroT_tmp_tmp > -TimeNeed2CalMax)&&(MacroT_tmp_tmp < TimeNeed2CalMax))
                                {   

                                    if((MacroT_tmp_tmp >=0)&&(X2L2Tmp > - X20New)&&(X2L2Tmp < X20New) || (MacroT_tmp_tmp <0)&&(X1L2Tmp > - X10New)&&(X1L2Tmp < X10New))
                                    {
                                        double* TSUMTmp = MacroMicroTheoryOut(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, kappaStar, TmpTmpX1, TmpTmpX2, BoostImageRes, coeffi, MicroLensCoorXY, NStar, -SkyLimit_Micro_file, SkyLimit_Micro_file, -SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
                        
                                        long TimeDelayIndexStart = long ((TSUMTmp[2] - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) - 10;
                                        long TimeDelayIndexEnd = long ((TSUMTmp[3] - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) + 10;
                                        for(long TimeDelayIndexTmp = TimeDelayIndexStart; TimeDelayIndexTmp <= TimeDelayIndexEnd; TimeDelayIndexTmp ++ )
                                        {
                                            if((TimeDelayIndexStart >=0)&&(TimeDelayIndexEnd < LenTimeDelay))
                                            {
                                                double TimeDelayValue = TimeDelayRange[TimeDelayIndexTmp];
                                                double TimeDelayLowerLimit = TimeDelayValue - TimeDelayEsp;
                                                double AreaTmp = 0;

                                                if((TimeDelayValue >= TSUMTmp[2])&&(TimeDelayValue <= TSUMTmp[0]))
                                                {
                                                    if(TimeDelayLowerLimit < TSUMTmp[2])
                                                    {
                                                        AreaTmp = TSUMTmp[4] * pow(TimeDelayValue - TSUMTmp[2], 2) / (TSUMTmp[0] - TSUMTmp[2]) / 2;
                                                    }
                                                    else
                                                    {
                                                        AreaTmp = TSUMTmp[4] * (2 * TimeDelayValue - 2 * TSUMTmp[2] - TimeDelayEsp) / (TSUMTmp[0] - TSUMTmp[2]) * TimeDelayEsp / 2;
                                                    }

                                                } 
                                                else if((TimeDelayValue > TSUMTmp[0])&&(TimeDelayValue <= TSUMTmp[1]))
                                                {
                                                    if(TimeDelayLowerLimit <= TSUMTmp[2])
                                                    {
                                                        AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                        + TSUMTmp[4] * (TimeDelayValue - TSUMTmp[0]);
                                                    }
                                                    else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                                    {
                                                        AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2 
                                                        + TSUMTmp[4] * (TimeDelayValue - TSUMTmp[0]); 
                                                    }
                                                    else
                                                    {
                                                        AreaTmp = TSUMTmp[4] * TimeDelayEsp;
                                                    }
                                                }
                                                else if((TimeDelayValue > TSUMTmp[1])&&(TimeDelayValue <= TSUMTmp[3]))
                                                {
                                                    if(TimeDelayLowerLimit <= TSUMTmp[2])
                                                    {
                                                        AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                        + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) 
                                                        + (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2;
                                                    }
                                                    else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                                    {
                                                        AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2
                                                            + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) + 
                                                            (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2; 
                                                    }
                                                    else if((TimeDelayLowerLimit <= TSUMTmp[1])&&(TimeDelayLowerLimit > TSUMTmp[0]))
                                                    {
                                                        AreaTmp = TSUMTmp[4] * (TSUMTmp[1] - TimeDelayLowerLimit) + 
                                                            (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2;
                                                    } 
                                                    else
                                                    {
                                                        AreaTmp = (TSUMTmp[4] * (TSUMTmp[3] - TimeDelayLowerLimit) / (TSUMTmp[3] - TSUMTmp[1]) + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * TimeDelayEsp / 2; 
                                                    }
                                                }
                                                else if(TimeDelayValue > TSUMTmp[3])
                                                {
                                                    if(TimeDelayLowerLimit <= TSUMTmp[2])
                                                    {
                                                        AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                        + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) 
                                                        + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2;
                                                    }
                                                    else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                                    {
                                                        AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2
                                                            + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0])
                                                            + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2; 
                                                    }
                                                    else if((TimeDelayLowerLimit <= TSUMTmp[1])&&(TimeDelayLowerLimit > TSUMTmp[0]))
                                                    {
                                                        AreaTmp = TSUMTmp[4] * (TSUMTmp[1] - TimeDelayLowerLimit)
                                                            + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2;
                                                    } 
                                                    else if((TimeDelayLowerLimit <= TSUMTmp[3])&&(TimeDelayLowerLimit > TSUMTmp[1]))
                                                    {
                                                        AreaTmp = TSUMTmp[4] * pow(TSUMTmp[3] - TimeDelayLowerLimit, 2) / (TSUMTmp[3] - TSUMTmp[1]) / 2; 
                                                    }
                                                    else
                                                    {
                                                        AreaTmp = 0;
                                                    }
                                                }
                                                
                                                TimeDelayArea2D[i][TimeDelayIndexTmp] += AreaTmp; //FIXME 这里可能有错误，比如说多个线程同时相加。
                                            }
                                        }       
                                        delete[] TSUMTmp;   
                                    }
                                }
                            }
                        }
                        delete[] X1SetBoost;
                        delete[] X2SetBoost;
                    }
                }
            
            }
            
            

        };

        auto run_tasks = [&] (int count) -> void {
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
                threads[CountIndex] = new std::thread (task, args[CountIndex]);
            }
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
                threads[CountIndex]->join();
                delete threads[CountIndex]; // 这句话我漏掉了... 你运行的时候记得加上...
            }
        };

        long JThread = 0;
        for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
        {
            args[JThread++] = ThreadIndex;
            if (JThread == ThreadCount) {
                run_tasks(ThreadCount);
                JThread = 0;
            }
        }
        run_tasks(JThread);
        
        for(long AreaIndex = 0; AreaIndex < LenTimeDelay; AreaIndex ++)
        {
            for(long AreaSumIndex = 0; AreaSumIndex < ThreadCount; AreaSumIndex ++)
            {
                TimeDelayArea[AreaIndex] += TimeDelayArea2D[AreaSumIndex][AreaIndex];
            }
        }


        fp_area.write(reinterpret_cast<char *>(TimeDelayArea), sizeof(double)*LenTimeDelay);
        fp_time.write(reinterpret_cast<char *>(TimeDelayRange), sizeof(double)*LenTimeDelay);
        

        
        delete[] TimeDelayRange;
        delete[] TimeDelayArea;
        delete[] TimeDelayArea2D;
        
    }
    //Maximum
    if((mur < 0)&&(mut < 0))
    {
        sprintf(LayersFile,"ResultMaximum_%2d/LayersFile_max_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        sprintf(TimeLengthFile,"ResultMaximum_%2d/TimeLength_max_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        LF.open(LayersFile, std::ofstream::binary);
        TLFile.open(TimeLengthFile, std::ofstream::binary);
        mur = - mur;
        mut = - mut;
        double Axisa = sqrt(1/coeffi/mur);
        double Axisb = sqrt(1/coeffi/mut);
        
        sprintf(AreaName,"ResultMaximum_%2d/adptive_Area_max_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift);
        sprintf(TimeName,"ResultMaximum_%2d/adptive_Time_max_%2.2f_%2.2f_%2.2f_%2.2f_%2.2f.bin", Max_mass, kappa, gamma, kappaStar_in, LensRedshift, SourceRedshift); 
        
        fp_area.open(AreaName, std::ofstream::binary);
        fp_time.open(TimeName, std::ofstream::binary);
        


        for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex += 1 )
        {
            double X1L2Tmp = X1SetL2[X1L2TmpIndex];

            // cout << RudeX1Tmp << endl;
            for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex += 1 )
            {
                double X2L2Tmp = X2SetL2[X2L2TmpIndex];
                double RudePsi = MacroMicroAndMinusSheet(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, X1L2Tmp, X2L2Tmp, kappaStar, MicroLensCoorXY, NStar, - SkyLimit_Micro_file, SkyLimit_Micro_file, - SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
                // cout << "RudeX1Tmp = " << RudeX1Tmp << endl;
                // cout << "RudeX2Tmp = " << RudeX2Tmp << endl;
                // cout << "coeffi = " << coeffi << endl;
                // cout << "RudePsi = " << RudePsi << endl;
                double MacroT_tmp = pow(X1L2Tmp,2)/(2*Axisa*Axisa) + pow(X2L2Tmp,2)/(2*Axisb*Axisb);
                if((MacroT_tmp > -TimeNeed2CalMax)&&(MacroT_tmp < TimeNeed2CalMax))
                {
                    if(RudePsi > RudeTimeDelayMinAndMax[1])
                    {
                        RudeTimeDelayMinAndMax[1] = RudePsi;
                        
                    }
                    if(RudePsi < RudeTimeDelayMinAndMax[0])
                    {
                        RudeTimeDelayMinAndMax[0] = RudePsi;
                        TestX1 = X1L2Tmp;
                        TestX2 = X2L2Tmp;
                        TestX1Index = X1L2TmpIndex;
                        TestX2Index = X2L2TmpIndex;
                        // cout << "Test = " << X1Set[TestX1Index] << endl;
                        
                    }
                }
            }
            
        }
        double RudeTimeDelayMin = RudeTimeDelayMinAndMax[0] * coeffi;
        double RudeTimeDelayMax = RudeTimeDelayMinAndMax[1] * coeffi + 0.1;
        cout << "Rude Minimum Time Delay is " << RudeTimeDelayMin << endl;
        cout << "Rude Minimum X1 = " << TestX1 << " X2 = " << TestX2 << endl;
        cout << "Rude Minimum X1 Index = " << TestX1Index << " X2 Index = " << TestX2Index << endl; 
        cout << "Rude Maximun Time Delay is " << RudeTimeDelayMax << endl;
        delete[] RudeTimeDelayMinAndMax;
        /*以上*/


        /*下面并行计算时间延迟曲线*/

        //the next is for time range
        double TimeDelayStart = RudeTimeDelayMin; 
        
        double TimeDelayEnd = RudeTimeDelayMax;
        //
        cout << "Time Maximum Boundary For Time Delay Curve " << TimeDelayEnd << endl;
        long LenTimeDelay = (TimeDelayEnd - TimeDelayStart + TimeDelayEsp)/TimeDelayEsp; 
        // the length of time delay area
        cout << "Time Delay Curve's Length = " << LenTimeDelay << endl; //print len_time
        TimeLength.push_back(LenTimeDelay);
        
        
        double* TimeDelayRange = new double [LenTimeDelay];
        double TimeDelayTmp = TimeDelayStart; //the start of the time delay calculation.
        //generate time array.
        for(long TimeDelayIndex = 0; TimeDelayIndex < LenTimeDelay; TimeDelayIndex ++ )
        {
        
            TimeDelayRange[TimeDelayIndex] = TimeDelayTmp;
            TimeDelayTmp += TimeDelayEsp;
            
            

        }
        cout << TimeDelayRange[LenTimeDelay - 1] << endl;
        //cout << int (sizeof(time_range)/sizeof(time_range[0])) << endl;


        double* TimeDelayArea = new double [LenTimeDelay](); //存放时间延迟面积的数组。
        //创建二维的数组用来存放时间延迟面积的数据，并初始化为0；
        double** TimeDelayArea2D = new double*[ThreadCount];
        for(int i = 0; i < ThreadCount; i++)
        {
            TimeDelayArea2D[i] = new double[LenTimeDelay]();
        }
        

        
        auto task = [&] (int i) -> void {

        
            long IndexLowLimit = i == 0 ? 0 : RivetingPoint[i - 1]; //三目运算符
            for(long X1X2L2TmpIndexOneDim = IndexLowLimit; X1X2L2TmpIndexOneDim <= RivetingPoint[i]; X1X2L2TmpIndexOneDim++ )
            { 
                double X1L2Tmp = X1SetL2OneDim[X1X2L2TmpIndexOneDim]; //算一下这一列的x值。
                double X2L2Tmp = X2SetL2OneDim[X1X2L2TmpIndexOneDim]; //算一下坐标的y值。
                // cout << "X2L2Tmp = " << X2L2Tmp << endl;
                double MacroT_tmp = pow(X1L2Tmp,2)/(2*Axisa*Axisa) + pow(X2L2Tmp,2)/(2*Axisb*Axisb);
                int LayersNumber_tmp = LayersNumber_Res[X1X2L2TmpIndexOneDim];
                
                if((MacroT_tmp > -TimeNeed2CalMax - TimeNeed2CalMax / 10)&&(MacroT_tmp < TimeNeed2CalMax + TimeNeed2CalMax / 10)) //加上了一个保护区
                //判断坐标是不是在双曲线区域内。
                {
                    double Resolution = ResolutionL2 / pow(2, LayersNumber_tmp);
                    // cout << "Resolution = " << Resolution << endl;
                
                    double BoostImageRes = Resolution;
                        
                    int Boost = int(ResolutionL2 / Resolution);
                    
                    /*加密以后的天区坐标*/
                    double* X1SetBoost = new double [Boost];
                    double* X2SetBoost = new double [Boost];
                    int X1BoostIndex = 0;
                    int X2BoostIndex = 0;
                    for(double X1SetBoostTmp = X1L2Tmp - ResolutionL2/2 + BoostImageRes/2; X1SetBoostTmp < X1L2Tmp + ResolutionL2/2; X1SetBoostTmp += BoostImageRes)
                    {
                        if(X1BoostIndex >= Boost)
                        {
                            
                            break;
                        }
                        else
                        {
                            X1SetBoost[X1BoostIndex] = X1SetBoostTmp;
                            X1BoostIndex += 1;
                        }
                    }
                    
                    for(double X2SetBoostTmp = X2L2Tmp - ResolutionL2/2 + BoostImageRes/2; X2SetBoostTmp < X2L2Tmp + ResolutionL2/2; X2SetBoostTmp += BoostImageRes)
                    {
                        if(X2BoostIndex >= Boost)
                        {
                            
                            break;
                        }
                        else
                        {
                            X2SetBoost[X2BoostIndex] = X2SetBoostTmp;
                            X2BoostIndex += 1;
                        }

                    }
                    

                    for(long X1Indextmp = 0; X1Indextmp < Boost; X1Indextmp++ )
                    { 
                        double TmpTmpX1 = X1SetBoost[X1Indextmp]; //算一下这一列的x值。
                        
                        for(long X2Indextmp = 0; X2Indextmp < Boost; X2Indextmp += 1)
                        {
                            double TmpTmpX2 = X2SetBoost[X2Indextmp]; //算一下坐标的y值。
                            double MacroT_tmp_tmp = pow(TmpTmpX1,2)/(2*Axisa*Axisa) + pow(TmpTmpX2,2)/(2*Axisb*Axisb);
                            if((MacroT_tmp_tmp > -TimeNeed2CalMax)&&(MacroT_tmp_tmp < TimeNeed2CalMax))
                            {   
                                double* TSUMTmp = MacroMicroTheoryOut(SkyLimit_Micro_file, SkyLimit_Micro_file, kappa, gamma, kappaStar, TmpTmpX1, TmpTmpX2, BoostImageRes, coeffi, MicroLensCoorXY, NStar, -SkyLimit_Micro_file, SkyLimit_Micro_file, -SkyLimit_Micro_file, SkyLimit_Micro_file, ResNearFieldMicroIndexSum, ResolutionL2, ResFarFieldAlphaAndCenterPotential, X1SetL2, X2SetL2, MassSample, AveMassTotal);
                                long TimeDelayIndexStart = long ((TSUMTmp[2] - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) - 10;
                                long TimeDelayIndexEnd = long ((TSUMTmp[3] - RudeTimeDelayMin + TimeDelayEsp)/TimeDelayEsp) + 10;
                                for(long TimeDelayIndexTmp = TimeDelayIndexStart; TimeDelayIndexTmp <= TimeDelayIndexEnd; TimeDelayIndexTmp ++ )
                                {
                                    if((TimeDelayIndexTmp >=0)&&(TimeDelayIndexTmp < LenTimeDelay))
                                    {
                                        double TimeDelayValue = TimeDelayRange[TimeDelayIndexTmp];
                                        double TimeDelayLowerLimit = TimeDelayValue - TimeDelayEsp;
                                        double AreaTmp = 0;

                                        if((TimeDelayValue >= TSUMTmp[2])&&(TimeDelayValue <= TSUMTmp[0]))
                                        {
                                            if(TimeDelayLowerLimit < TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * pow(TimeDelayValue - TSUMTmp[2], 2) / (TSUMTmp[0] - TSUMTmp[2]) / 2;
                                            }
                                            else
                                            {
                                                AreaTmp = TSUMTmp[4] * (2 * TimeDelayValue - 2 * TSUMTmp[2] - TimeDelayEsp) / (TSUMTmp[0] - TSUMTmp[2]) * TimeDelayEsp / 2;
                                            }

                                        } 
                                        else if((TimeDelayValue > TSUMTmp[0])&&(TimeDelayValue <= TSUMTmp[1]))
                                        {
                                            if(TimeDelayLowerLimit <= TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                + TSUMTmp[4] * (TimeDelayValue - TSUMTmp[0]);
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2 
                                                + TSUMTmp[4] * (TimeDelayValue - TSUMTmp[0]); 
                                            }
                                            else
                                            {
                                                AreaTmp = TSUMTmp[4] * TimeDelayEsp;
                                            }
                                        }
                                        else if((TimeDelayValue > TSUMTmp[1])&&(TimeDelayValue <= TSUMTmp[3]))
                                        {
                                            if(TimeDelayLowerLimit <= TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) 
                                                + (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2;
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2
                                                    + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) + 
                                                    (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2; 
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[1])&&(TimeDelayLowerLimit > TSUMTmp[0]))
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[1] - TimeDelayLowerLimit) + 
                                                    (TSUMTmp[4] + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * (TimeDelayValue - TSUMTmp[1]) / 2;
                                            } 
                                            else
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TSUMTmp[3] - TimeDelayLowerLimit) / (TSUMTmp[3] - TSUMTmp[1]) + TSUMTmp[4] * (TSUMTmp[3] - TimeDelayValue) / (TSUMTmp[3] - TSUMTmp[1])) * TimeDelayEsp / 2; 
                                            }
                                        }
                                        else if(TimeDelayValue > TSUMTmp[3])
                                        {
                                            if(TimeDelayLowerLimit <= TSUMTmp[2])
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[0] - TSUMTmp[2]) / 2 
                                                + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0]) 
                                                + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2;
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[0])&&(TimeDelayLowerLimit > TSUMTmp[2]))
                                            {
                                                AreaTmp = (TSUMTmp[4] * (TimeDelayLowerLimit - TSUMTmp[2]) / (TSUMTmp[0] - TSUMTmp[2]) + TSUMTmp[4]) * (TSUMTmp[0] - TimeDelayLowerLimit) / 2
                                                    + TSUMTmp[4] * (TSUMTmp[1] - TSUMTmp[0])
                                                    + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2; 
                                            }
                                            else if((TimeDelayLowerLimit <= TSUMTmp[1])&&(TimeDelayLowerLimit > TSUMTmp[0]))
                                            {
                                                AreaTmp = TSUMTmp[4] * (TSUMTmp[1] - TimeDelayLowerLimit)
                                                    + TSUMTmp[4] * (TSUMTmp[3] - TSUMTmp[1]) / 2;
                                            } 
                                            else if((TimeDelayLowerLimit <= TSUMTmp[3])&&(TimeDelayLowerLimit > TSUMTmp[1]))
                                            {
                                                AreaTmp = TSUMTmp[4] * pow(TSUMTmp[3] - TimeDelayLowerLimit, 2) / (TSUMTmp[3] - TSUMTmp[1]) / 2; 
                                            }
                                            else
                                            {
                                                AreaTmp = 0;
                                            }
                                        }
                                        
                                        TimeDelayArea2D[i][TimeDelayIndexTmp] += AreaTmp; //FIXME 这里可能有错误，比如说多个线程同时相加。
                                    }
                                }       
                                delete[] TSUMTmp;   
                            }
                        }
                    }
                    delete[] X1SetBoost;
                    delete[] X2SetBoost;


                
                }
            
            }
            
            

        };

        auto run_tasks = [&] (int count) -> void {
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex ) {
                threads[CountIndex] = new std::thread (task, args[CountIndex]);
            }
            for (int CountIndex = 0; CountIndex < count; ++ CountIndex) {
                threads[CountIndex]->join();
                delete threads[CountIndex]; // 这句话我漏掉了... 你运行的时候记得加上...
            }
        };

        long JThread = 0;
        for(int ThreadIndex = 0; ThreadIndex < ThreadCount ; ThreadIndex += 1)
        {
            args[JThread++] = ThreadIndex;
            if (JThread == ThreadCount) {
                run_tasks(ThreadCount);
                JThread = 0;
            }
        }
        run_tasks(JThread);
        
        for(long AreaIndex = 0; AreaIndex < LenTimeDelay; AreaIndex ++)
        {
            for(long AreaSumIndex = 0; AreaSumIndex < ThreadCount; AreaSumIndex ++)
            {
                TimeDelayArea[AreaIndex] += TimeDelayArea2D[AreaSumIndex][AreaIndex];
            }
        }

        fp_area.write(reinterpret_cast<char *>(TimeDelayArea), sizeof(double)*LenTimeDelay);
        fp_time.write(reinterpret_cast<char *>(TimeDelayRange), sizeof(double)*LenTimeDelay);
        
        

        delete[] TimeDelayRange;
        delete[] TimeDelayArea;
        delete[] TimeDelayArea2D;
    }



    fp_area.close();
    fp_time.close();
    
    //print run time.
    double CostTime = time(NULL) - TimeStart;
    cout << "CostTimeAdaptiveCal =  " << CostTime - CostTimeLayerNum << " s" << endl;
    delete[] X1SetL2; 
    delete[] X2SetL2;
    for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex++ )
    {
        for(long X2L2TmpIndex = 0; X2L2TmpIndex < X2L2Length; X2L2TmpIndex ++ )
        {
            delete[] ResFarFieldAlphaAndCenterPotential[X1L2TmpIndex][X2L2TmpIndex];
            delete[] ResNearFieldMicroIndexSum[X1L2TmpIndex][X2L2TmpIndex];
        }
        delete[] ResFarFieldAlphaAndCenterPotential[X1L2TmpIndex];
        delete[] ResNearFieldMicroIndexSum[X1L2TmpIndex];
    
    }
    delete[] ResFarFieldAlphaAndCenterPotential;
    delete[] ResNearFieldMicroIndexSum;
    delete[] RivetingPoint_L2Coeffi;
    for(long X1L2TmpIndex = 0; X1L2TmpIndex < X1L2Length; X1L2TmpIndex ++ )
    {
        delete[] Level2GridMicroPotential[X1L2TmpIndex];
    }
    delete[] Level2GridMicroPotential; 
    delete[] RivetingPoint_LayerNum;
    delete[] RivetingPoint;
    delete[] PreOutPut;
    delete[] MassSample;
    delete[] MicroLensCoorXY;

    TLFile.write(reinterpret_cast<char *>(&TimeLength[0]), sizeof(long)*1);

    TLFile.close();
    LF.write((char*)&LayersNumber_Res[0], sizeof(long)*TotalLengthLayerNum_Res);
    LF.close();
    
    

    return 0;
}

