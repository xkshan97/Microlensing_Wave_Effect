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
#include <numeric>
#include "../../spline.h"


using std::sin;
using std::cos;
using std::tan;
using std::exp; 
using std::pow;
using std::log;
using std::sqrt;
using std::atan;
using std::acos;
using std::log10;

using namespace std;


#define PI acos(-1)



double SalpeterIMF(double StellarMass, double StellarMassMin, double StellarMassMax)
{
    return pow(StellarMass, -2.35) / (1 / (-1.35) * pow(StellarMassMax, -1.35) - 1 / (-1.35) * pow(StellarMassMin, -1.35));
}


double NormalizationChabrierIMF(double StellarMassMin, double StellarMassMax)
{
    double Step = 0.001;
    int StepNum = (StellarMassMax - StellarMassMin) / Step;
    double IMFPDF = 0;
    double StellarMassTmp = StellarMassMin;
    for(int StepI = 0; StepI < StepNum; StepI ++)
    {
        if(StellarMassTmp < 1)
        {
            IMFPDF += 0.158 / StellarMassTmp * exp(-1/2. * pow((log10(StellarMassTmp) - log10(0.079)) / 0.69, 2));
            StellarMassTmp += Step;
        }
        else
        {
            IMFPDF += 0.0443 * pow(StellarMassTmp, -2.3);
            StellarMassTmp += Step;
        }
    }
    return IMFPDF * Step;
}

double ChabrierIMF(double StellarMass, double StellarMassMin, double StellarMassMax, double NormalizationChabrierIMFValue)
{
    if(StellarMass < 1)
    {
        return 0.158 / StellarMass * exp(-1/2. * pow((log10(StellarMass) - log10(0.079)) / 0.69, 2)) / NormalizationChabrierIMFValue;
    }
    else
    {
        return 0.0443 * pow(StellarMass, -2.3) / NormalizationChabrierIMFValue;
    }
}



double* SampleResult(int NStarStellar, int NStarRemnant, string IMFType)
// int main()
{
    
    double* OutPut = new double [NStarStellar + NStarRemnant]; //前面的储存质量
    
    /*随机种子*/
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> Uniform01(0, 1); //用来生成0～1之间的均匀分布。

    /*差值remnant质量函数，输入一个质量给出一个PDF*/
    ifstream RemnantMFFile("/disk1/home/shanxk/work/Paper3_adaptive_Micro_field_release/SampleMethod/Remnant_MF.csv", ios::in);
    string lineRemnantMF;
   
    vector<double> RemnantMFMass;
    vector<double> RemnantMFPDF;
    int ReadIndex = 0;
    while (getline(RemnantMFFile, lineRemnantMF))
    {
        istringstream sin(lineRemnantMF);
        
        string lineRemnantMF_tmp;
        while (getline(sin, lineRemnantMF_tmp, ','))
        {
            if(ReadIndex % 2 == 0)
            {
                RemnantMFMass.push_back(stod(lineRemnantMF_tmp));
            }
            else
            {
                RemnantMFPDF.push_back(stod(lineRemnantMF_tmp));
            }
        }
        ReadIndex += 1;
    }
    
    tk::spline InterRemnantMF(RemnantMFMass,RemnantMFPDF);

    /*Remnant拒绝接受采样*/
    double ConstCRemnant = *max_element(RemnantMFPDF.begin(), RemnantMFPDF.end());
    double RemnantMassMin = *min_element(RemnantMFMass.begin(), RemnantMFMass.end()); 
    double RemnantMassMax = *max_element(RemnantMFMass.begin(), RemnantMFMass.end()); //88; //28; //
    cout << "ConstCRemnant = " << ConstCRemnant << endl;
    cout << "RemnantMassMin = " << RemnantMassMin << endl;
    cout << "RemnantMassMax = " << RemnantMassMax << endl;
    
    std::uniform_real_distribution<> GRemnant(RemnantMassMin, RemnantMassMax); //用来生成候选样本。
    for(int NStarRemnantIndex = 0; NStarRemnantIndex < NStarRemnant; NStarRemnantIndex ++ )
    {
        double Uniform01Sample = Uniform01(gen);
        double GRemnantUniform = GRemnant(gen);
        // while(GRemnantUniform > RemnantMassMax || GRemnantUniform < RemnantMassMin)
        // {
        //     GRemnantUniform = GRemnant(gen);
        // }
        // cout << "!!!!!!!!!!!!!" << endl;
        while(Uniform01Sample > InterRemnantMF(GRemnantUniform) / ConstCRemnant)
        {
            Uniform01Sample = Uniform01(gen);
            GRemnantUniform = GRemnant(gen);
        }
        OutPut[NStarStellar + NStarRemnantIndex] = GRemnantUniform;

    }

    // /*检验差值误差*/
    // double* RelativeError = new double [100];
    // double* InterPDF = new double [100];
    // for(int TestErrorI = 0; TestErrorI < 100; TestErrorI ++ )
    // {
    //     RelativeError[TestErrorI] = (InterRemnantMF(RemnantMFMass[TestErrorI]) - RemnantMFPDF[TestErrorI]) / RemnantMFPDF[TestErrorI];
    //     InterPDF[TestErrorI] = InterRemnantMF(RemnantMFMass[TestErrorI]); 
    //     printf("RelativeError = %1.10f and InterPDF = %1.10f\n", RelativeError[TestErrorI], InterPDF[TestErrorI]);
    // }
    // /*以上*/


    // /*检验IMF归一化*/
    // cout << NormalizationChabrierIMF(0.01, 1.5) << endl;
    // cout << NormalizationChabrierIMF(0.01, 0.08) << endl;

    /*检验IMF PDF*/
    double StellarMassMin = 0.08;
    double StellarMassMax = 1.5;
    if(IMFType == "Chabrier")
    {
        double NormalizationChabrierIMFValue = NormalizationChabrierIMF(StellarMassMin, StellarMassMax);
        double StepTest = 0.01;
        int StepNum = (StellarMassMax - StellarMassMin) / StepTest;
        double* IMFPDFTest = new double [StepNum];
        for(int StepI = 0; StepI < StepNum; StepI ++ )
        {
            IMFPDFTest[StepI] = ChabrierIMF(0.08 + StepI * StepTest, StellarMassMin, StellarMassMax, NormalizationChabrierIMFValue);
        }
        // ofstream IMFPDFTestOP;
        // char IMFPDFTestFile[100];
        // sprintf(IMFPDFTestFile,"/disk1/home/shanxk/work/Paper3_adaptive_Micro_field_release/SampleMethod/IMFPDFTest.bin");
        // IMFPDFTestOP.open(IMFPDFTestFile, std::ofstream::binary);
        // IMFPDFTestOP.write(reinterpret_cast<char *>(IMFPDFTest), sizeof(double)*StepNum);
        // IMFPDFTestOP.close();
        // /*以上*/
        
        /*Stellar拒绝接受采样*/
        double ConstCStellar = *max_element(IMFPDFTest, IMFPDFTest + StepNum) * (StellarMassMax - StellarMassMin);
        cout << "ConstCStellar = " << ConstCStellar << endl;
        
        std::uniform_real_distribution<> GStellar(StellarMassMin, StellarMassMax); //用来生成候选样本。
        for(int NStarStellarIndex = 0; NStarStellarIndex < NStarStellar; NStarStellarIndex ++ )
        {
            double Uniform01Sample = Uniform01(gen);
            double GStellarUniform = GStellar(gen);
            while(Uniform01Sample > ChabrierIMF(GStellarUniform, StellarMassMin, StellarMassMax, NormalizationChabrierIMFValue) / ConstCStellar * (StellarMassMax - StellarMassMin))
            {
                Uniform01Sample = Uniform01(gen);
                GStellarUniform = GStellar(gen);
            }
            OutPut[NStarStellarIndex] = GStellarUniform;

        }
    }
    if(IMFType == "Salpeter")
    {
        double StepTest = 0.01;
        int StepNum = (StellarMassMax - StellarMassMin) / StepTest;
        double* IMFPDFTest = new double [StepNum];
        for(int StepI = 0; StepI < StepNum; StepI ++ )
        {
            IMFPDFTest[StepI] = SalpeterIMF(0.08 + StepI * StepTest, StellarMassMin, StellarMassMax);
        }
        // ofstream IMFPDFTestOP;
        // char IMFPDFTestFile[100];
        // sprintf(IMFPDFTestFile,"/disk1/home/shanxk/work/Paper3_adaptive_Micro_field_release/SampleMethod/IMFPDFTest.bin");
        // IMFPDFTestOP.open(IMFPDFTestFile, std::ofstream::binary);
        // IMFPDFTestOP.write(reinterpret_cast<char *>(IMFPDFTest), sizeof(double)*StepNum);
        // IMFPDFTestOP.close();
        // /*以上*/
        
        /*Stellar拒绝接受采样*/
        double ConstCStellar = *max_element(IMFPDFTest, IMFPDFTest + StepNum) * (StellarMassMax - StellarMassMin);
        cout << "ConstCStellar = " << ConstCStellar << endl;
        
        std::uniform_real_distribution<> GStellar(StellarMassMin, StellarMassMax); //用来生成候选样本。
        for(int NStarStellarIndex = 0; NStarStellarIndex < NStarStellar; NStarStellarIndex ++ )
        {
            double Uniform01Sample = Uniform01(gen);
            double GStellarUniform = GStellar(gen);
            while(Uniform01Sample > SalpeterIMF(GStellarUniform, StellarMassMin, StellarMassMax) / ConstCStellar * (StellarMassMax - StellarMassMin))
            {
                Uniform01Sample = Uniform01(gen);
                GStellarUniform = GStellar(gen);
            }
            OutPut[NStarStellarIndex] = GStellarUniform;

        }

    }

    /*将采样的数据保存成文件*/
    ofstream IMFSampleTestOP;
    char IMFSampleTestFile[100];
    sprintf(IMFSampleTestFile,"/disk1/home/shanxk/work/Paper3_adaptive_Micro_field_release/SampleMethod/SampleTest.bin");
    IMFSampleTestOP.open(IMFSampleTestFile, std::ofstream::binary);
    IMFSampleTestOP.write(reinterpret_cast<char *>(OutPut), sizeof(double)*(NStarStellar + NStarRemnant));
    IMFSampleTestOP.close();

    return OutPut;

}

// int main()
// {
//     double* MassSample = SampleResult(1997, 46);
//     double SumMass = 0;
//     for(int SumMassIndex = 0; SumMassIndex < 1997 + 46; SumMassIndex ++ )
//     {
//         SumMass += MassSample[SumMassIndex];
//     }
//     double kappaStar = SumMass / 0.35224277772719608803 * PI / (2 * 36.5843) / (2 * 36.5843);
//     printf("Numerical exact kappaStar = %.10f\n", kappaStar);
// }