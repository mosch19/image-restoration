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

void showDFT(Mat& source)
{
    Mat splitArray[2] = {Mat::zeros(source.size(), CV_32F), Mat::zeros(source.size(), CV_32F)};
    split(source, splitArray);
    
    Mat dftmagnitude;

    magnitude(splitArray[0], splitArray[1], dftmagnitude);
    dftmagnitude += Scalar::all(1);

    log(dftmagnitude, dftmagnitude);
    normalize(dftmagnitude, dftmagnitude, 0, 1, CV_MINMAX);

    imshow("DFT", dftmagnitude);
    waitKey();

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

/*
float MSE(float image1, float image2)
{
    float sum_squared = 0;
    float mse = 0;
    float area = image1.width * image1.height;

    cimg_forXY(image1,x,y)
    {
        int pix1 = image1(x,y);
        int pix2 = image2(x,y);

        int error = pix2 - pix1;
        sum_squared += (error * error);
    }

    mse = sum_squared/area;
    return mse;
}
*/

int main()
{
    CImg<float> original;
    CImg<float> image;
    CImg<float> image2;

    string location = "/home/mosch/Documents/Restoration/lena_full.jpg";
    original.load(location.c_str());
    image = original.get_RGBtoYCbCr().get_channel(0);

    float width = image.width;
    float height = image.height;
    float area = width * height;
    int frame[9];

    //mess it up
    image.noise(40);
    //image.blur(2.5);
    image2 = image;

    //vector<float> input[int(area)];
    //int z = 0;

    for (int x = 1; x < int(width) - 1; x++)
    {
        for (int y = 1; y < int(height) - 1; y++)
        {
            frame[0] = int(image2(x - 1 ,y - 1));
            frame[1] = int(image2(x, y - 1));
            frame[2] = int(image2(x + 1, y - 1));
            frame[3] = int(image2(x - 1, y));
            frame[4] = int(image2(x, y));
            frame[5] = int(image2(x + 1, y));
            frame[6] = int(image2(x - 1, y + 1));
            frame[7] = int(image2(x, y + 1));
            frame[8] = int(image2(x + 1, y + 1));

            insertionSort(frame);
            image2(x,y) = frame[4];
        }
    }

    //calculate MSE between two images
    float
    
    sum_squared = 0;
    float mse = 0;
    

    cimg_forXY(image,x,y)
    {
        int pix1 = image(x,y);
        int pix2 = image2(x,y);

        int error = pix2 - pix1;
        sum_squared += (error * error);
    }

    mse = sum_squared/area;

    Mat origin = imread("/home/mosch/Documents/Restoration/lena_full.jpg", CV_LOAD_IMAGE_GRAYSCALE);
    Mat origin_float;
    origin.convertTo(origin_float, CV_32FC1, 1.0 / 255.0);
    Mat dftOUT;

    takeDFT(origin_float, dftOUT);
    showDFT(dftOUT);

    CImgDisplay org_disp(original, "Original"), main_disp(image, "Low Quality Lena"), draw_disp(image2, "Median Filtered");
    cout << "MSE: " << mse << endl;

    //display the images

    while(!draw_disp.is_closed)
    {
        draw_disp.wait();
    }


    return 0;
}