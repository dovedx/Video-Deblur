
#ifndef _mexutil_h
#define _mexutil_h
#include <opencv2/core/core.hpp>

#include <opencv2/opencv.hpp>
#include "patchmatch.h"
#include <iostream>
//#include "nn.h"
using namespace std;
using namespace cv;


unsigned char* ucharMat_to_array(Mat& A);
float* float_Mat_to_array(Mat& A);
int*  Mat_to_array(Mat& A);
Mat Array_to_mat(int a[], int h, int w);
BITMAP *convert_bitmap(Mat& A);
//BITMAP *convert_bitmapf(Mat& A);
BITMAP *convert_field(int h,int w,Params *p,int A[], int bw, int bh, int &nclip, int trim_patch=1);
BITMAP *float_convert_field(int h,int w,Params *p,float A[], int bw, int bh, int &nclip, int trim_patch=1);

Mat bitmap_to_Mat(BITMAP *a);
Mat float_Array_to_mat(float a[], int h, int w);
int* array_col_to_row(int* a,int h,int w);

float* float_array_col_to_row(float* a,int h,int w);
template<class T> 
inline void clip_value(T *p, T low, T high) { 
	if(*p < low) *p = low; 
	else if(*p > high) *p = high; 
}

// from MATLAB's column major to C++'s row major

#endif
