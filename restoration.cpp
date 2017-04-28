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
    CImg<float> image;
    CImg<float> image2;
    string location = "/home/mosch/Documents/Restoration/lena_full.jpg";
    image.load(location.c_str());
    image2 = image.get_RGBtoYCbCr().get_channel(0);
    CImgDisplay main_disp(image, "Lena"), draw_disp(image2, "Gray Lena");

    while(!main_disp.is_closed)
    {
        main_disp.wait();
    }


    return 0;
}