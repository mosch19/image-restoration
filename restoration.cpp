#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include <complex>
#include "CImg.h"

#define cimg_use_jpeg

using namespace cimg_library;
using namespace cv;
using namespace std;

//Compile with:  g++ `pkg-config --cflags opencv` restoration.cpp `pkg-config --libs opencv` -o runt -L/usr/X11R6/lib -lm -lpthread -lX11

void takeDFT(Mat& source, Mat& destination)
{
    Mat dftIN[2] = {source, Mat::zeros(source.size(), CV_32F)};
    Mat dftReady;
    merge(dftIN, 2, dftReady);
    Mat dftOUT;
    
    dft(dftReady, dftOUT, DFT_COMPLEX_OUTPUT);

    destination = dftOUT;
}

Mat showDFT(Mat& source)
{
    Mat splitArray[2] = {Mat::zeros(source.size(), CV_32F), Mat::zeros(source.size(), CV_32F)};
    split(source, splitArray);
    
    Mat dftmagnitude;

    magnitude(splitArray[0], splitArray[1], dftmagnitude);
    dftmagnitude += Scalar::all(1);

    log(dftmagnitude, dftmagnitude);
    normalize(dftmagnitude, dftmagnitude, 0, 1, CV_MINMAX);

    return dftmagnitude;
    //imshow("DFT", dftmagnitude);
    //waitKey();

}

void insertionSort(int frame[])
{
    int temp;
    int i = 0;
    int j = 0;

    for(i = 0; i < 9; i++)
    {
        temp = frame[i];

        for(j = i-1; j >= 0 && temp < frame[j]; j--)
        {
            frame[j+1] = frame[j];
        }
        
        frame[j+1] = temp;
    }
}

double MSE(CImg<float> source, CImg<float> comp)
{
    int area = int(source.width) * int(source.height);
    double sum_squared = 0;
    double mse = 0;

    cimg_forXY(source,x,y)
    {
        double pix1 = source(x,y);
        double pix2 = comp(x,y);

        double error = pix2 - pix1;
        sum_squared += (error * error);
    }

    mse = sum_squared/area;

    return mse;
}

CImg<float> medianFilter(CImg<float> source)
{
    int frame[9];
    int width = source.width;
    int height = source.height;
    CImg<float> destination;
    destination = source;

    for (int x = 1; x < int(width) - 1; x++)
    {
        for (int y = 1; y < int(height) - 1; y++)
        {
            frame[0] = int(source(x - 1 ,y - 1));
            frame[1] = int(source(x, y - 1));
            frame[2] = int(source(x + 1, y - 1));
            frame[3] = int(source(x - 1, y));
            frame[4] = int(source(x, y));
            frame[5] = int(source(x + 1, y));
            frame[6] = int(source(x - 1, y + 1));
            frame[7] = int(source(x, y + 1));
            frame[8] = int(source(x + 1, y + 1));

            insertionSort(frame);
            destination(x,y) = frame[4];
        }
    }

    return destination;
}

CImg<float> meanFilter(CImg<float> source)
{
    int frame[9];
    float mean = 0;
    int width = source.width;
    int height = source.height;
    CImg<float> destination;
    destination = source;

    for (int x = 1; x < int(width) - 1; x++)
    {
        for (int y = 1; y < int(height) - 1; y++)
        {
            mean = 0;

            frame[0] = int(source(x - 1 ,y - 1));
            frame[1] = int(source(x, y - 1));
            frame[2] = int(source(x + 1, y - 1));
            frame[3] = int(source(x - 1, y));
            frame[4] = int(source(x, y));
            frame[5] = int(source(x + 1, y));
            frame[6] = int(source(x - 1, y + 1));
            frame[7] = int(source(x, y + 1));
            frame[8] = int(source(x + 1, y + 1));

            for (int z = 0; z < 9; z++)
            {
                mean += frame[z];
            }
    
            destination(x,y) = (mean / 9);
        }
    }

    return destination;
}

CImg<float> gaussianBlur(CImg<float> source)
{
    double frame[25];
    double mean = 0;
    int width = source.width;
    int height = source.height;
    CImg<float> destination;
    destination = source;

    for (int x = 3; x < int(width) - 3; x++)
    {
        for (int y = 3; y < int(height) - 3; y++)
        {
            mean = 0.0;

            frame[0] = source(x - 2 ,y - 2) * .003765;
            frame[1] = source(x - 1 ,y - 2) * .015019;
            frame[2] = source(x - 0 ,y - 2) * .023792;
            frame[3] = source(x + 1 ,y - 2) * .015019;
            frame[4] = source(x + 2 ,y - 2) * .003765;

            frame[5] = source(x - 2 ,y - 1) * .015019;
            frame[6] = source(x - 1 ,y - 1) * .059912;
            frame[7] = source(x - 0 ,y - 1) * .094907;
            frame[8] = source(x + 1 ,y - 1) * .059912;
            frame[9] = source(x + 2 ,y - 1) * .015019;
            
            frame[10] = source(x - 2 ,y - 0) * .023792;
            frame[11] = source(x - 1 ,y - 0) * .094907;
            frame[12] = source(x - 0 ,y - 0) * .150342;
            frame[13] = source(x + 1 ,y - 0) * .094907;
            frame[14] = source(x + 2 ,y - 0) * .023792;
            
            frame[15] = source(x - 2 ,y + 1) * .015019;
            frame[16] = source(x - 1 ,y + 1) * .059912;
            frame[17] = source(x - 0 ,y + 1) * .094907;
            frame[18] = source(x + 1 ,y + 1) * .059912;
            frame[19] = source(x + 2 ,y + 1) * .015019;
            
            frame[20] = source(x - 2 ,y + 2) * .003765;
            frame[21] = source(x - 1 ,y + 2) * .015019;
            frame[22] = source(x - 0 ,y + 2) * .023792;
            frame[23] = source(x + 1 ,y + 2) * .015019;
            frame[24] = source(x + 2 ,y + 2) * .003765;
            
            for (int z = 0; z < 25; z++)
            {
                mean += frame[z];
            }
            
            //cout << mean << endl;
            destination(x,y) = float(mean);// / 25);
            //cout << destination(x,y) << endl;
        }
    }

    return destination;
}
/*
CImg<float> gaussianKernel(CImg<float> source)
{
    int width = source.width;
    int height = source.height;
    CImg<float> destination;
    destination = source;

    double mean = 0;
    double sigma = 1.0;
    double sum = 0.0;
    double gkern[5][5];
    double top, bottom = 0;

    for (int x = 1; x < int(width) - 1; x++)
    {
        for (int y = 1; y < int(height) - 1; y++)
        {
            //solve for the gaussian kernal
            for (int j = -2; j < 3; j++)
            {
                for (int k = -2; k < 3; k++)
                {
                    top = (k*k) + (j*j);
                    bottom = 2(sigma*sigma);
                    gkern[j + 2][k + 2] = exp(-(top/bottom)) / (bottom * M_PI);
                    sum += gkern[j + 2][k + 2];
                }
            }

            //kernal normalized
            for(int l = 0; l < 5; l++)
            {
                for(int m = 0; m < 5; m++)
                {
                    gkernal[l][m] /= sum;
                }
            }

            //apply kernel
            for(int n = 0; n < 5; n++)
            {
                for(int o = 0; o < 5; o++)
                {
                    mean += gkernal[n][o];
                }
            }
    
            destination(x,y) = (mean / 25);
        }
    }

    return destination;
}
*/
int main()
{
    CImg<float> original;
    CImg<float> image;
    CImg<float> image2;
    CImg<float> image3;
    CImg<float> image4;

    string location = "/home/mosch/Documents/Restoration/lena_full.jpg";
    original.load(location.c_str());
    image = original.get_RGBtoYCbCr().get_channel(0);

    float width = image.width;
    float height = image.height;
    float area = width * height;

    //mess it up
    image.noise(60,0); //30, 2 crazy good for median filter
    //image.blur(2.5);

    image2 = medianFilter(image);
    image3 = meanFilter(image);
    image4 = gaussianBlur(image);

    //calculate MSE between two images
    cout << "MSE: " << MSE(image2, image2) << endl;   

    Mat origin = imread("/home/mosch/Documents/Restoration/lena_full.jpg", CV_LOAD_IMAGE_GRAYSCALE);
    Mat origin_float;
    origin.convertTo(origin_float, CV_32FC1, 1.0 / 255.0);
    Mat dftOUT;

    takeDFT(origin_float, dftOUT);
    showDFT(dftOUT);

    CImgDisplay noise_disp(image, "Original"), median_disp(image2, "Median Filtered"), mean_disp(image3, "Mean Filtered"), gaus_disp(image4, "Gaussian Blur");
    //CImgDisplay org_disp(original, "Original"), main_disp(image, "Blur/Noise"), draw_disp(image2, "Median Filtered");
    //imshow("DFT Magnitude", showDFT(dftOUT));
    //waitKey();

    //display the images until one closed
    while(!mean_disp.is_closed)
    {
        mean_disp.wait();
    }


    return 0;
}