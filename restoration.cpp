#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "CImg.h"
#define cimg_use_jpeg
using namespace cimg_library;
using namespace std;

int main()
{
    CImg<float> original;
    CImg<float> image;
    CImg<float> image2;

    float MSE = 0;

    //int height, width = 0;
    //height = original.height;
    //width = original.width;

    //CImg<int> test(int width, int height);

    string location = "/home/mosch/Documents/Restoration/lena_full.jpg";
    original.load(location.c_str());
    image = original.get_RGBtoYCbCr().get_channel(0);
    image2 = image;

    cimg_forXY(image, x, y)
    {
        if (image(x,y) < 100)
        {
            image(x,y) = 0;
            cout << "value: " << image(x,y) << endl;
        }
    }
    
    //mess it up
    image2.noise(40);
    image2.blur(2.5);
    CImgDisplay org_disp(original, "Original"), main_disp(image, "Greyscale Lena"), draw_disp(image2, "Blur/Noise Lena");

    //display the images
    while(!draw_disp.is_closed)
    {
        draw_disp.wait();
    }


    return 0;
}