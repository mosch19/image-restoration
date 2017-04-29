#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <complex>
#include "CImg.h"

#define cimg_use_jpeg

using namespace cimg_library;
using namespace std;

void computeDft(const vector<double> &inreal, const vector<double> &inimag,
		vector<double> &outreal, vector<double> &outimag) {
	
	unsigned int n = inreal.size();
	for (unsigned int k = 0; k < n; k++) {  /* For each output element */
		double sumreal = 0;
		double sumimag = 0;
		for (unsigned int t = 0; t < n; t++) {  /* For each input element */
			double angle = 2 * M_PI * t * k / n;
			sumreal +=  inreal[t] * cos(angle) + inimag[t] * sin(angle);
			sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
		}
		outreal[k] = sumreal;
		outimag[k] = sumimag;
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

    float area = image.width * image.height;
    

    string location = "/home/mosch/Documents/Restoration/lena_full.jpg";
    original.load(location.c_str());
    image = original.get_RGBtoYCbCr().get_channel(0);
    image2 = image;

    vector<float> input[int(area)];
    int z = 0;

    cimg_forXY(image, x, y)
    {
        input[z] = image(x,y);
        z++;
        
        //if (image(x,y) < 100)
        //{
        //    image(x,y) = 0;
        //    cout << "value: " << image(x,y) << endl;
        //}
    
    }

    computeDft(input);

    //mess it up
    image2.noise(40);
    image2.blur(2.5);

    //calculate MSE between two images
    float sum_squared = 0;
    float mse = 0;
    

    cimg_forXY(image,x,y)
    {
        int pix1 = image(x,y);
        int pix2 = image2(x,y);

        int error = pix2 - pix1;
        sum_squared += (error * error);
    }

    mse = sum_squared/area;

    CImgDisplay org_disp(original, "Original"), main_disp(image, "Greyscale Lena"), draw_disp(image2, "Messed Up");
    cout << "MSE: " << mse << endl;

/*
    //calculate dft 
    cimg_forX(image,x)
    {
        complex<double> sum(0.0,0.0);

        cimg_forY(image,y)
        {
            int integers = -2*x*y;
            complex<float> expo(0.0, M_PI/area*(double)integers);
            sum += image(y) * exp(expo);
        }

        cout << "DFT: " << abs(sum) << endl;    
    }
*/

/*
    complex<float> image3 = CImg<float> image;
    complex<double> output_seq[int(area)];
    double pi2 = 2.0 * M_PI;
    double angleTerm,cosineA,sineA;
    double invs = 1.0 / area;

    cimg_forX(image3,x)
    {
        output_seq[x] = 0;

        cimg_forY(image3,y)
        {
            angleTerm = pi2 * y * x * invs;
            cosineA = cos(angleTerm);
            sineA = sin(angleTerm);
            
            output_seq[y].real += image3(y).real * cosineA - image3(y).imag * sineA;
            output_seq[y].imag += image3(y).real * sineA + image3(y).imag * cosineA;
        }
        output_seq[x] *= invs;
    }

    for (int z = 0; z < area; z++)
    {
        cout << "DFT: " << output_seq[z] << endl;
    }
*/
    //display the images
    while(!draw_disp.is_closed)
    {
        draw_disp.wait();
    }


    return 0;
}