#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include <complex>

using namespace cv;
using namespace std;

#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include <complex>


using namespace cv;
using namespace std;



Mat padd_image(Mat I){
	Mat padded;
	int m = getOptimalDFTSize( I.rows );
	int n = getOptimalDFTSize( I.cols ); // on the border add zero pixels
	copyMakeBorder(I, padded, 0, m - I.rows, 0, n - I.cols, BORDER_CONSTANT, Scalar::all(0));
	return padded;
}

Mat get_dft(Mat I){
	Mat image;
	I.convertTo(image, CV_32F);
	Mat planes[] = {Mat_<float>(image), Mat::zeros(image.size(), CV_32F)};
	Mat complexI;
	merge(planes, 2, complexI);
	dft(complexI, complexI);
	return complexI;
}

Mat get_spectrum(Mat I){
	Mat complexI = get_dft(I);
	Mat planes[2];
	split(complexI, planes);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
	magnitude(planes[0], planes[1], planes[0]);// planes[0] = magnitude
	Mat magI = planes[0];
	multiply(magI,magI,magI);
	return magI;
}



Mat rand_noise(Mat I, int stddev){
	Mat noise = Mat::zeros(I.rows, I.cols, CV_32F);
	randn(noise,Scalar::all(0), Scalar::all(stddev));
	return noise;
}

Mat with_noise(Mat image, int stddev){
	Mat noise(image.rows, image.cols, CV_8U);
	rand_noise(image, stddev).convertTo(noise, CV_8U);
	Mat noisy = image.clone();
	noisy += noise;
	return noisy;
}


Mat wiener2(Mat I, Mat image_spectrum, int noise_stddev){
	Mat padded = padd_image(I);
	Mat noise = rand_noise(padded, noise_stddev);
	Mat noise_spectrum = get_spectrum(noise);

	Scalar padded_mean = mean(padded);

	Mat planes[2];
	Mat complexI = get_dft(padded);
	split(complexI, planes);	// planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))

	Mat factor = image_spectrum / (image_spectrum + noise_spectrum);
	multiply(planes[0],factor,planes[0]);
	multiply(planes[1],factor,planes[1]);


	merge(planes, 2, complexI);
	idft(complexI, complexI);
	split(complexI, planes);
//	normalize(planes[0], planes[0], 0, 128, CV_MINMAX );
	Scalar enhanced_mean = mean(planes[0]);
	double norm_factor =  padded_mean.val[0] / enhanced_mean.val[0];
	multiply(planes[0],norm_factor, planes[0]);
	Mat normalized;
	planes[0].convertTo(normalized, CV_8UC1);
	return normalized;
}



int main(){

	int noise_stddev;
	cout << "Choose noise standard dev" << endl;
	cin >> noise_stddev;

	Mat image = imread("lena.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Image", image);
	//Mat blurredImage1;
	Mat blurredImage;
	GaussianBlur(image, blurredImage, Size (7,7), 15, 0, 4);
	imshow("Blurred Image", blurredImage);
	
	
	Mat padded;
	int m = getOptimalDFTSize(image.rows);
	int n = getOptimalDFTSize(image.cols);
	copyMakeBorder(image, padded, 0, m-image.rows, 0, n-image.cols, BORDER_CONSTANT, Scalar::all(0));
	Mat planes [] = {Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F)};
	Mat complexI;
	merge(planes,2,complexI);
	dft(complexI,complexI);
	split(complexI, planes);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
	magnitude(planes[0], planes[1], planes[0]);// planes[0] = magnitude
	Mat magI = planes[0];
	magI += Scalar::all(1);
	log(magI,magI);
	magI = magI(Rect(0, 0, magI.cols & -2, magI.rows & -2));
	int cx = magI.cols/2;	
	int cy = magI.rows/2;

	Mat q0(magI, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
	Mat q1(magI, Rect(cx, 0, cx, cy));  // Top-Right
	Mat q2(magI, Rect(0, cy, cx, cy));  // Bottom-Left
	Mat q3(magI, Rect(cx, cy, cx, cy)); // Bottom-Right

	Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
	q2.copyTo(q1);
	tmp.copyTo(q2);

	normalize(magI, magI, 0, 1, CV_MINMAX);

	Mat padded2;
	int a = getOptimalDFTSize(blurredImage.rows);
	int b = getOptimalDFTSize(blurredImage.cols);
	copyMakeBorder(blurredImage, padded2, 0, a-blurredImage.rows, 0, b-blurredImage.cols, BORDER_CONSTANT, Scalar::all(0));
	Mat planes2 [] = {Mat_<float>(padded2), Mat::zeros(padded2.size(), CV_32F)};
	Mat complexI2;
	merge(planes2,2,complexI2);
	dft(complexI2,complexI2);
	split(complexI2, planes2);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
	magnitude(planes2[0], planes2[1], planes2[0]);// planes[0] = magnitude
	Mat magI2 = planes2[0];
	magI2 += Scalar::all(1);
	log(magI2,magI2);
	magI2 = magI2(Rect(0, 0, magI2.cols & -2, magI2.rows & -2));
	int cx2 = magI2.cols/2;	
	int cy2 = magI2.rows/2;

	Mat q02(magI2, Rect(0, 0, cx2, cy2));   // Top-Left - Create a ROI per quadrant
	Mat q12(magI2, Rect(cx2, 0, cx2, cy2));  // Top-Right
	Mat q22(magI2, Rect(0, cy2, cx2, cy2));  // Bottom-Left
	Mat q32(magI2, Rect(cx2, cy2, cx2, cy2)); // Bottom-Right

	Mat tmp2;                           // swap quadrants (Top-Left with Bottom-Right)
	q02.copyTo(tmp2);
	q32.copyTo(q02);
	tmp2.copyTo(q32);

	q12.copyTo(tmp2);                    // swap quadrant (Top-Right with Bottom-Left)
	q22.copyTo(q12);
	tmp2.copyTo(q22);

	normalize(magI2, magI2, 0, 1, CV_MINMAX);
	
	imshow("DFT of original", magI);
	imshow("DFT of blurred", magI2);
	
	Mat blurkernal;
	divide(complexI2, complexI, blurkernal);
	Mat visualBlurKernal;
	divide(magI2, magI, visualBlurKernal);
	imshow ("DFT of oringal / blurred", visualBlurKernal);
	Mat Xf;
	divide(complexI2,blurkernal, Xf);
	Mat visualXf;
	divide (magI2, visualBlurKernal, visualXf);
	imshow ("Xf", visualXf);

	Mat inverseDFT;
	idft(Xf, inverseDFT, DFT_INVERSE | DFT_REAL_OUTPUT);
	normalize(inverseDFT, inverseDFT, 0, 1, CV_MINMAX);
	//cout << inverseDFT << endl;
	imshow("Reconstructed", inverseDFT);

	Mat padded3=padd_image(image);
	Mat noisy;
	noisy = with_noise(padded, noise_stddev);

	//GaussianBlur(noisy, noisy, Size (3,3), 100, 0, 4);

	Mat sample(padded3.rows, padded3.cols, CV_8U);
	resize(image, sample, sample.size());
	Mat spectrum = get_spectrum(sample);
	Mat enhanced = wiener2(noisy, spectrum, noise_stddev);
	

	imshow("Noisy image before Wiener", noisy);
	imshow("Enhanced Image after Wiener", enhanced);



	waitKey(0);
	destroyAllWindows();
	
	system("pause");
    return 0;
}

