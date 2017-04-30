#include <iostream>
#include <math.h>
#include <vector> 


using namespace std;
 
#define PI 3.14159265
int n;
 
int main(int argc, char **argv)
{
    //At this point, instead of asking the user to enter the size, we pass the size as the dimensions of the image we are restoring.
    //It should also be noted that this code considers the width and the height to be the same. We just need to change that to height and width instead of n
    cout << "Enter the size: ";
    cin >> n;
    //double inputData[n][n];
	  vector<vector<double>> inputData;
  	inputData.resize(n,vector<double>(n,n));

    //This is where the user inputs the pixel values of the array. I think instead of this, we just have the for loops that go through each 
    //pixel of the image and extract its amplitude
    cout << "Enter the 2D elements ";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> inputData[i][j];
 
    //double realOut[n][n];
	  vector<vector<double>> realOut;
	  realOut.resize(n,vector<double>(n,n));

	  //double imagOut[n][n];
	  vector<vector<double>> imagOut;
	  imagOut.resize(n,vector<double>(n,n));

    //double amplitudeOut[n][n];
	  vector<vector<double>> amplitudeOut;
	  amplitudeOut.resize(n,vector<double>(n,n));

    //this is another spot that needs to slightly be changed. Instead of n it will be the actual height/width of the image
    int height = n;
    int width = n;
 
    // Two outer loops iterate on output data.
    for (int yWave = 0; yWave < height; yWave++)
    {
        for (int xWave = 0; xWave < width; xWave++)
        {
            // Two inner loops iterate on input data.
            for (int ySpace = 0; ySpace < height; ySpace++)
            {
                for (int xSpace = 0; xSpace < width; xSpace++)
                {
                    // Compute real, imag, and ampltude.
                    realOut[yWave][xWave] += (inputData[ySpace][xSpace] * cos(
                            2 * PI * ((1.0 * xWave * xSpace / width) + (1.0
                                    * yWave * ySpace / height)))) / sqrt(
                            width * height);
                    imagOut[yWave][xWave] -= (inputData[ySpace][xSpace] * sin(
                            2 * PI * ((1.0 * xWave * xSpace / width) + (1.0
                                    * yWave * ySpace / height)))) / sqrt(
                            width * height);
                    amplitudeOut[yWave][xWave] = sqrt(
                            realOut[yWave][xWave] * realOut[yWave][xWave]
                                    + imagOut[yWave][xWave]
                                            * imagOut[yWave][xWave]);
                }
                cout << realOut[yWave][xWave] << " + " << imagOut[yWave][xWave]
                        << " i (" << amplitudeOut[yWave][xWave] << ")\n";
            }
        }
    }
	system("pause");
}
