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
    // //读取透镜的配置
    // ifstream kappafile("/disk1/home/shanxk/work/Paper4_CE_Modify/SampleResult/kappa.csv", ios::in);
    // ifstream gammafile("/disk1/home/shanxk/work/Paper4_CE_Modify/SampleResult/gamma.csv", ios::in);
    // ifstream kappastarfile("/disk1/home/shanxk/work/Paper4_CE_Modify/SampleResult/kappa_s.csv", ios::in);
    // ifstream AcceptLensindexfile("/disk1/home/shanxk/work/Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv", ios::in);
    // ifstream lenszfile("/disk1/home/shanxk/work/Paper4_CE_Modify/SampleResult/lens_z.csv", ios::in);
    // ifstream imagenumfile("/disk1/home/shanxk/work/Paper4_CE_Modify/SampleResult/imagenumber.csv", ios::in); 
    // ifstream SNROnlyMacrofile("/disk1/home/shanxk/work/Paper4_CE_Modify/Lensed_SampleResult/SNR_network_only_macro.csv", ios::in);
    // string linekappa;
    // string linegamma;
    // string linekappastar;
    // string lineindex;
    // string linelensz;
    // string lineimagenum;
    // string linesnr;
    // vector<double> lenskappa;
    // vector<double> lensgamma;
    // vector<double> lenskappastar;
    // vector<double> lensindex;
    // vector<double> lensz;
    // vector<double> lensimagenum;
    // vector<double> lenssnr;
    // int index_num = 0;
    // while (getline(kappafile, linekappa))
    // {
    //     istringstream sin(linekappa);
        
    //     string kappa_tmp;
    //     while (getline(sin, kappa_tmp, ','))
    //     {
    //         lenskappa.push_back(stod(kappa_tmp));
    //     }
    // }
    // while (getline(gammafile, linegamma))
    // {
    //     istringstream sin(linegamma);
        
    //     string gamma_tmp;
    //     while (getline(sin, gamma_tmp, ','))
    //     {
    //         lensgamma.push_back(stod(gamma_tmp));
    //     }
    // }
    // while (getline(kappastarfile, linekappastar))
    // {
    //     istringstream sin(linekappastar);
        
    //     string kappastar_tmp;
    //     while (getline(sin, kappastar_tmp, ','))
    //     {
    //         lenskappastar.push_back(stod(kappastar_tmp));
    //     }
    // }
    // while (getline(AcceptLensindexfile, lineindex))
    // {
    //     istringstream sin(lineindex);
        
    //     string index_tmp;
    //     while (getline(sin, index_tmp, ','))
    //     {
    //         lensindex.push_back(stod(index_tmp));
    //         index_num ++;
    //     }
    // }
    // while (getline(lenszfile, linelensz))
    // {
    //     istringstream sin(linelensz);
        
    //     string lensz_tmp;
    //     while (getline(sin, lensz_tmp, ','))
    //     {
    //         lensz.push_back(stod(lensz_tmp));
    //     }
    // }
    // while (getline(imagenumfile, lineimagenum))
    // {
    //     istringstream sin(lineimagenum);
        
    //     string imagenum_tmp;
    //     while (getline(sin, imagenum_tmp, ','))
    //     {
    //         lensimagenum.push_back(stod(imagenum_tmp));
    //     }
    // }
    // while (getline(SNROnlyMacrofile, linesnr))
    // {
    //     istringstream sin(linesnr);
        
    //     string lenssnr_tmp;
    //     while (getline(sin, lenssnr_tmp, ','))
    //     {
    //         lenssnr.push_back(stod(lenssnr_tmp));
    //     }
    // }
    
    // //以上
    
    
    // int IndexSumimage = 0;
    // for(int ImageIndex = 0; ImageIndex < 300; ImageIndex ++)
    // {
    //     int IndexInAccept = lensindex[ImageIndex];
    //     double LensRedshift = lensz[ImageIndex];
    //     for(int subImageIndex = 0; subImageIndex < lensimagenum[ImageIndex]; subImageIndex ++ )
    //     {
    //         // if(1)
    //         if(lenssnr[IndexSumimage] >= 12)
    //         {
    //             double kappa = lenskappa[IndexSumimage];
    //             double gamma = lensgamma[IndexSumimage];
    //             double kappaStar_stellar = lenskappastar[IndexSumimage];
    //             if(1.2 * kappaStar_stellar > kappa)
    //             {
    //                 kappaStar_stellar = kappa / 1.2;
    //                 cout << "f_* > 1 and reassignment kappa_* = kappa" << endl;
    //             }

    //             printf("Input Kappa star = %2.20f\n", 1.2 * kappaStar_stellar);
    //             cout << "IndexSumimage = " << IndexSumimage << endl;
    //             cout << "ImageIndex = " << ImageIndex << endl;
    //             cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl; 
    //             MainDiffraction(kappa, gamma, 1.2 * kappaStar_stellar, LensRedshift, 1, 250);
                
    //         }
    //         IndexSumimage += 1;
    //     }
    // }
    double kappa = 0.45;
    double gamma = 0.45;
    double kappaStar_stellar = 0.3;
    double LensRedshift = 0.5;
    double SourceRedshift = 1;
    int thread_count = 100;
    
    MainDiffraction(kappa, gamma, kappaStar_stellar, LensRedshift, SourceRedshift, thread_count);
}
