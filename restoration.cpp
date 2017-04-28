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
    string location = "/home/mosch/Documents/C++/lena_full.jpg";
    image.load(location.c_str());
    CImgDisplay main_disp(image, "Lena");

    while(!main_disp.is_closed)
    {
        main_disp.wait();
    }


    return 0;
}