
#include "mexutil.h"
//#include "nn.h"
#include <string.h>

// stack RGB to a single 32bit integer

float* float_Mat_to_array(Mat& A)
{
  int h = A.rows;
  int w = A.cols;
  float* A_Array;
  A_Array = new float[3 * h*w];

  int num = 0, num1 = h*w, num2 = 2 * h*w;
  for (int j = 0; j < w; j++)
  {
	  for (int i = 0; i < h; i++)
	  {
		  A_Array[num] = A.at<Vec3f>(i, j)[0];

		  A_Array[num1 ] = A.at<Vec3f>(i, j)[1];

		  A_Array[num2 ] = A.at<Vec3f>(i, j)[2];

		  num++;
		  num1++;
		  num2++;

	  }

	  //cout << num << endl;
  }
  
  return A_Array;
}

int* Mat_to_array(Mat& A)
{
  int h = A.rows;
  int w = A.cols;
  int* A_Array;
  A_Array = new int[3 * h*w];

  int num = 0, num1 = h*w, num2 = 2 * h*w;
  for (int j = 0; j < w; j++)
  {
	  for (int i = 0; i < h; i++)
	  {
		  A_Array[num] = A.at<Vec3i>(i, j)[0];

		  A_Array[num1 ] = A.at<Vec3i>(i, j)[1];

		  A_Array[num2 ] = A.at<Vec3i>(i, j)[2];

		  num++;
		  num1++;
		  num2++;

	  }

	  //cout << num << endl;
  }
  
  return A_Array;
}

unsigned char* ucharMat_to_array(Mat& A)
{
	int h = A.rows;
	int w = A.cols;
	unsigned char* A_Array;
	A_Array = new unsigned char[3 * h*w];

	int num = 0, num1 = h*w, num2 = 2*h*w;
	for (int j = 0; j < w; j++)
	{
		for (int i = 0; i < h; i++)
		{
			A_Array[num] = A.at<Vec3b>(i, j)[0];
			
			A_Array[num1] = A.at<Vec3b>(i, j)[1];
		
			A_Array[num2] = A.at<Vec3b>(i, j)[2];
		
			num++;
			num1++;
			num2++;
			
		}
		
		//cout << num << endl;
	}
	return A_Array;
}

Mat Array_to_mat(int a[],int h,int w)
{
    Mat MAT_A(h,w,CV_64FC3);
    int num =0,num1=h*w,num2=2*h*w;
    for(int j=0;j<w;j++)
    {
      for(int i=0;i<h;i++)
      {
        MAT_A.at<Vec3d>(i,j)[0]=(double)a[num];
       
		MAT_A.at<Vec3d>(i, j)[1] = (double)a[num1];
       
		MAT_A.at<Vec3d>(i, j)[2] = (double)a[num2];
        num++;
		num1++;
		num2++;
      }
    }
    return MAT_A;
}


Mat float_Array_to_mat(float a[], int h, int w)
{
	Mat MAT_A(h, w, CV_32FC3);
	
	int num = 0, num1 = h*w, num2 = 2 * h*w;
	for (int j = 0; j < w; j++)
	{
		for (int i = 0; i < h; i++)
		{
			MAT_A.at<Vec3f>(i, j)[0] = a[num];

			MAT_A.at<Vec3f>(i, j)[1] =a[num1];

			MAT_A.at<Vec3f>(i, j)[2] = a[num2];
			num++;
			num1++;
			num2++;
		}
	}
	return MAT_A;
}


BITMAP *convert_bitmap( Mat& A) 
{
	int h = A.rows;
	int w = A.cols;
	int offset = h * w;
	int offset2 = offset << 1; // h * w * 2

    unsigned char *data = (unsigned char*) ucharMat_to_array(A);
    BITMAP *ans = create_bitmap(w, h);
	for (int y = 0; y < h; y++) 
    {
      int *row = (int *) ans->line[y];
	for (int x = 0; x < w; x++) 
	{
		unsigned char *base = data + y + x * h;
		int r = base[0];
		int g = base[offset];
		int b = base[offset2];
		clip_value<int>(&r, 0, 255);
		clip_value<int>(&g, 0, 255);
		clip_value<int>(&b, 0, 255);
		row[x] = r | (g << 8) | (b << 16);
	}
    }

	cout << "finished" << endl;
    return ans;
 
}

//BITMAP *convert_bitmapf(Mat& A) 
//{
//	int w = A.cols;
//	int h = A.rows;
//    unsigned char *data = (unsigned char *)Mat_to_array(A);
//    BITMAP *ans = create_bitmap(w, h);
//    for (int y = 0; y < h; y++) {
//      float *row = (float *) ans->line[y];
//      for (int x = 0; x < w; x++) {
//        row[x] = data[y+x*h];
//      }
//    }
//    return ans;
//   
//}

// stack (x, y) to an integer by concatenate "yx" (assume y and x has at most 12 bit, or 4096 values)
BITMAP *convert_field(int h,int w,Params *p,int A[], int bw, int bh, int &nclip, int trim_patch)
{
    nclip = 0;

	int offset = h * w;


  int bew = trim_patch ? (bw - p->patch_w + 1): bw;
  int beh = trim_patch ? (bh - p->patch_w + 1): bh;

  unsigned int *data = (unsigned int*) A;
  BITMAP *ann = create_bitmap(w, h);
 
  //annd = create_bitmap(w, h);
  for (int y = 0; y < h; y++) 
  {
    int *ann_row = (int *) ann->line[y];
    //int *annd_row = (int *) annd->line[y];
    for (int x = 0; x < w; x++) 
	{
	  unsigned int *base = data + y + x*h; 
      int xp = base[0];
      int yp = base[offset];
	  if ((unsigned)xp >= (unsigned)bew || (unsigned)yp >= (unsigned)beh)  
	  {
        nclip++;
		clip_value<int>(&xp, 0, bew-1);
		clip_value<int>(&yp, 0, beh-1);
      }
      //int dp = data[(y+x*h)+2*w*h];
      ann_row[x] = XY_TO_INT(xp, yp);
      //annd_row[x] = dp;
    }
  }
  int a = 1;
  return ann;
}

BITMAP *float_convert_field(int h,int w,Params *p,float A[], int bw, int bh, int &nclip, int trim_patch)
{
   nclip = 0;

   int offset = h * w;

  int bew = trim_patch ? (bw - p->patch_w + 1): bw;
  int beh = trim_patch ? (bh - p->patch_w + 1): bh;

  unsigned int *data = (unsigned int*) A;
  BITMAP *ann = create_bitmap(w, h);
 
  //annd = create_bitmap(w, h);
  for (int y = 0; y < h; y++) 
  {
    int *ann_row = (int *) ann->line[y];
    //int *annd_row = (int *) annd->line[y];
    for (int x = 0; x < w; x++) 
	{
	  unsigned int *base = data + y + x*h; 
      int xp = base[0];
      int yp = base[offset];
	  if ((unsigned)xp >= (unsigned)bew || (unsigned)yp >= (unsigned)beh)  
	  {
        nclip++;
		clip_value<int>(&xp, 0, bew-1);
		clip_value<int>(&yp, 0, beh-1);
      }
      //int dp = data[(y+x*h)+2*w*h];
      ann_row[x] = XY_TO_INT(xp, yp);
      //annd_row[x] = dp;
    }
  }
  int a = 1;
  return ann;

}

Mat bitmap_to_Mat(BITMAP *a) 
{
  // mwSize dims[3] = { a->h, a->w, 3 };
  // mxArray *ans = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
  // unsigned char *data = (unsigned char *) mxGetData(ans);
  int height =a->h;
  int width = a->w;
  int *A=new int[3*height*width];
  int* data = A;

  int offset = a->w * a->h;
  int offset2 = offset << 1;
  int *rchan = &data[0];
  int *gchan = &data[offset];
  int *bchan = &data[offset2];
  for (int y = 0; y < a->h; y++)
  {
    int *row = (int *) a->line[y];
    for (int x = 0; x < a->w; x++) 
    {
      int c = row[x];
	  int pos = y + x * a->h;
      rchan[pos] = c&255;
      gchan[pos] = (c>>8)&255;
      bchan[pos] = (c>>16);
    }
  }
  Mat MAT_A = Array_to_mat(A,height,width);
  delete []A;
  return MAT_A;
}

int* array_col_to_row(int* a,int h,int w)
{
	int* a_row = new int[h*w];
	int num = 0;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			a_row[num] = a[j*h + i];
			num++;
		}
	}
	return a_row;
}

float* float_array_col_to_row(float* a,int h,int w)
{
	float* a_row = new float[h*w];
	int num = 0;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			a_row[num] = a[j*h + i];
			num++;
		}
	}
	return a_row;

}