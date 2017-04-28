#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#define cimg_use_jpeg
#include "CImg.h"
using namespace cimg_library;

int main()
{
	CImg<unsigned char> image("lenna.jpg"); // , visu(500, 400, 1, 3, 0);

	return 0;
}