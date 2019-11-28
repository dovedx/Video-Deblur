
#include "mexutil.h"

#include "patchmatch.h"
using namespace std;
using namespace cv;

//void init_params(Params *p);

#define MODE_IMAGE  0


#ifndef M_PI
#define M_PI  3.1415926535897932384626433832795
#endif

extern int xform_scale_table[NUM_SCALES];

void creat_ann(Mat& A, Mat& B,Mat& ann_1,Mat& sim_ann_1,int h,int w)
{

	int aw = -1, ah = -1, bw = -1, bh = -1;
	BITMAP *a = NULL, *b = NULL, *ann_prev = NULL, *ann_window = NULL, *awinsize = NULL;
	Params *p = new Params();
	RecomposeParams *rp = new RecomposeParams();
	//BITMAP *borig = NULL;
	int mode = 0;
	if (mode == MODE_IMAGE) 
	{
		a = convert_bitmap(A);
		b = convert_bitmap(B);
	//	borig = b;
		aw = a->w; ah = a->h;//图片的长和宽
		bw = b->w; bh = b->h;
	}
	//Mat C = bitmap_to_Mat(a);
	//C.convertTo(C, CV_8UC3, 1, 0);
	//imshow("C", C);
	//waitKey(0);
	double *win_size = NULL;
	BITMAP *amask = NULL, *bmask = NULL;
	double scalemin = 0.8, scalemax = 1.2;		 // The product of these must be one.
	int sim_mode = 1;			//计算相识度的模式,加不加旋转和缩放，0不加，1加
	p->algo = 0;

	p->patch_w = 80;
	p->nn_iters = 5;
	p->rs_max = INT_MAX;
	p->rs_min = 1;
	p->rs_ratio = 0.5;
	p->rs_iters = 1;
	p->do_propagate = 1;
	p->knn = 0;
	int enrich_mode = 0;
	if (enrich_mode)
	{
		int nn_iters = p->nn_iters;
		p->enrich_iters = nn_iters / 2;
		p->nn_iters = 2;
	}
	init_params(p);
	if (sim_mode == 1)
	{
		init_xform_tables(scalemin, scalemax, 1);
	}
	RegionMasks *amaskm = amask ? new RegionMasks(p, amask) : NULL;
	BITMAP *ann = NULL; // NN field
	BITMAP *annd_final = NULL; // NN patch distance field
	BITMAP *ann_sim_final = NULL;
	if (mode == MODE_IMAGE)
	{
		if (sim_mode)
		{
			BITMAP *ann_sim = NULL;

			ann = sim_init_nn(p, a, b, ann_sim);
			BITMAP *annd = sim_init_dist(p, a, b, ann, ann_sim);
			sim_nn(p, a, b, ann, ann_sim, annd);
			if (ann_prev)
			{
				cout << "when searching over rotations+scales, previous guess is not supported" << endl;
			}
			annd_final = annd;
			ann_sim_final = ann_sim;
		}
		// else
		// {
		// 	ann = init_nn(p, a, b, bmask, NULL, amaskm, 1, ann_window, awinsize);
		// 	BITMAP *annd = init_dist(p, a, b, ann, bmask, NULL, amaskm);
		// 	nn(p, a, b, ann, annd, amaskm, bmask, 0, 0, rp, 0, 0, 0, NULL, p->cores, ann_window, awinsize);
		// 	if (ann_prev) minnn(p, a, b, ann, annd, ann_prev, bmask, 0, 0, rp, NULL, amaskm, p->cores);
		// 	annd_final = annd;
		// }
	}
	
	if (ann_sim_final)
	{
		
		float *data = new float[5 * h*w];
		float *xchan = &data[0];
		float *ychan = &data[aw*ah];
		float *dchan = &data[2 * aw*ah];
		float *tchan = &data[3 * aw*ah];
		float *schan = &data[4 * aw*ah];
		//double angle_scale = 2.0*M_PI / NUM_ANGLES;
		double angle_scale =(10.0*2.0*M_PI/360.0)/NUM_ANGLES;
		//cout<<NUM_ANGLES<<"::"<<NUM_SCALES<<endl;
			//cout<<"fff"<< endl;
		for (int y = 0; y < ah; y++)
		{
			int *ann_row = (int *)ann->line[y];
			int *annd_row = (int *)annd_final->line[y];
			int *ann_sim_row = ann_sim_final ? (int *)ann_sim_final->line[y] : NULL;
			for (int x = 0; x < aw; x++)
			{
				int pp = ann_row[x];
				int pos = y + x * ah;
				xchan[pos] = INT_TO_X(pp);
				ychan[pos] = INT_TO_Y(pp);
				dchan[pos] = annd_row[x];
				if (ann_sim_final)
				{
				int v = ann_sim_row[x];
				
				int tval = INT_TO_Y(v)&(NUM_ANGLES - 1);//取前4位
				int sval = INT_TO_X(v);//取后12位
				
				tchan[pos] = tval*angle_scale+(350*2.0*M_PI/360.0) ;
				//cout<<"tval"<<tval<<endl;
				//cout<<"sval:"<<sval<<endl;
				//cout<<sizeof(xform_scale_table)<<endl;
				schan[pos] = xform_scale_table[sval] * (1.0 / 65536.0);
				//cout<<"k"<<endl;
				}
			}
		}

		float* a1 = &data[0];
		float* a2 = &data[aw*ah];
		float* a3 = &data[2*aw*ah];
		float* a4 = &data[3 * aw*ah];
		float* a5 = &data[4 * aw*ah];

		float* a_row1 = float_array_col_to_row(a1,h,w);
		float* a_row2 = float_array_col_to_row(a2, h, w);
		float* a_row3 = float_array_col_to_row(a4, h, w);
		float* a_row4 = float_array_col_to_row(a5, h, w);

		Mat image1(h, w, CV_32FC3, Scalar(0.0, 0.0, 0.0));
		Mat image2(h, w, CV_32FC3, Scalar::all(0.0));
		vector<Mat> channels;
		vector<Mat> channels_2;
		split(image1,channels);
		split(image2,channels_2);
		//float* a_row4 = float_array_col_to_row(a5, h, w);

		Mat channel1(h,w, CV_32FC1, a_row1);
		Mat channel2(h,w, CV_32FC1, a_row2);
		Mat channel3(h,w, CV_32FC1, Scalar(0));

		Mat channel4(h, w, CV_32FC1, a_row3);
		Mat channel5(h, w, CV_32FC1, a_row4);

		channel1.copyTo(channels.at(0));
		channel2.copyTo(channels.at(1));
		channel3.copyTo(channels.at(2));

		channel4.copyTo(channels_2.at(0));
		//cout<<channel4(Range(100,110),Range(100,110));
		channel5.copyTo(channels_2.at(1));
		//cout<<channel5(Range(100,110),Range(100,110));
		channel3.copyTo(channels_2.at(2));

		merge(channels, image1);
		merge(channels_2, image2);
		ann_1 = image1;
		sim_ann_1 = image2;
		delete[]data;
		
		destroy_region_masks(amaskm);
		//return image;
	}
	// else
	// {
	// 	int *data = new int[3 * h*w];
	// 	int *xchan = &data[0];
	// 	int *ychan = &data[aw*ah];
	// 	int *dchan = &data[2 * aw*ah];
	// 	for (int y = 0; y < ah; y++)
	// 	{
	// 		int *ann_row = (int *)ann->line[y];
	// 		int *annd_row = (int *)annd_final->line[y];
	// 		for (int x = 0; x < aw; x++)
	// 		{
	// 			int pp = ann_row[x];
	// 			int pos = y + x * ah;
	// 			xchan[pos] = INT_TO_X(pp);
	// 			ychan[pos] = INT_TO_Y(pp);
	// 			dchan[pos] = annd_row[x];
	// 		}
	// 	}
	// 	int* a1 = &data[0];
	// 	int* a2 = &data[aw*ah];
	// 	int* a3 = &data[2*aw*ah];
		
	// 	int* a_row1 = array_col_to_row(a1,h,w);
	// 	int* a_row2 = array_col_to_row(a2, h, w);
	// 	int* a_row3 = array_col_to_row(a3, h, w);
	// 	Mat image(h, w, CV_32SC3, Scalar(0, 0, 0));
	// 	vector<Mat> channals;
	// 	split(image, channals);
		
	// 	Mat channel1(h, w, CV_32SC1, a_row1);
	// 	Mat channel2(h, w, CV_32SC1, a_row2);
	// 	Mat channel3(h, w, CV_32SC1, a_row3);
	// 	//cout << channel2(Range(10, 50), Range(480, 485)) << endl;
	// /*	imshow("chann", channel1);
	// 	waitKey(0);*/
	// 	channel1.copyTo(channals.at(0));
	// 	channel2.copyTo(channals.at(1));
	// 	channel3.copyTo(channals.at(2));
	// 	Mat mergeImage;
	// 	merge(channals, image);
	// 	ann_1 = image;
	// 	delete[]data;
	// 	destroy_region_masks(amaskm);
	// 	//cout << mergeImage(Range(100, 105), Range(100, 105)) << endl;
	// 	//system("pause");
	// 	//return image;
	// }
	
}

// Mat Vote(Mat& ann, Mat& sim_ann,Mat& B,int h,int w)
// {
// 	BITMAP *b = convert_bitmap(B);
// 	Params *p = new Params();
// 	p->patch_w = 50;
// 	BITMAP *bmask = NULL, *bweight = NULL, *amask = NULL, *aweight = NULL, *ainit = NULL;
// 	double coherence_weight = 1, complete_weight = 1;

// 	int aclip = 0, bclip = 0;

// 	int nclip = aclip + bclip;
// 	RegionMasks *amaskm = amask ? new RegionMasks(p, amask) : NULL;

// 	ann.convertTo(ann,CV_32SC3,1,0);
// 	sim_ann.convertTo(sim_ann,CV_32SC3,1,0);
// 	int* A2 = Mat_to_array(ann);
// 	BITMAP *ann1 = convert_field(h, w, p, A2, b->w, b->h, aclip);

// 	//if(sim_ann！=NULL)
// 	//{
// 		cout<<"aim_ann!=NULL"<<endl;
// 		int* A3 = Mat_to_array(sim_ann);
// 		BITMAP* sim_ann1 = convert_field(h,w,p,A3,b->w,b->h,aclip);
// 		BITMAP *a = sim_vote(p, b, ann1, sim_ann1,NULL, NULL,bmask, bweight, coherence_weight, complete_weight, amaskm,aweight, ainit, NULL, NULL, 1);
// 	//}
// 	//else
// 	//{
// 	//	BITMAP *a = vote(p, b, ann1, NULL, bmask, bweight, coherence_weight, complete_weight, amaskm,aweight, ainit, NULL, NULL, 1);
// 	//}
	
	
// 	Mat result = bitmap_to_Mat(a);
// 	result.convertTo(result, CV_8UC3, 1, 0);
// 	cout << result(Range(200, 210), Range(200, 203)) << endl;
	

// 	delete p;
// 	delete[]A2;
// 	destroy_bitmap(a);
// 	destroy_bitmap(b);
// 	destroy_bitmap(bmask);
// 	destroy_bitmap(bweight);
// 	destroy_bitmap(aweight);
// 	destroy_bitmap(ainit);
// 	destroy_region_masks(amaskm);

// 	return result;
// }

Mat rotateImage1(Mat img, int degree)
{
	degree = -degree;
	double angle = degree  * CV_PI / 180.; // 弧度
	double a = sin(angle), b = cos(angle);
	int width = img.cols;
	int height = img.rows;
	int width_rotate = int(height * fabs(a) + width * fabs(b));
	int height_rotate = int(width * fabs(a) + height * fabs(b));
	//Ðý×ªÊý×émap
	// [ m0  m1  m2 ] ===>  [ A11  A12   b1 ]
	// [ m3  m4  m5 ] ===>  [ A21  A22   b2 ]
	float map[6];
	Mat map_matrix = Mat(2, 3, CV_32F, map);
	// 旋转矩阵
	
	CvPoint2D32f center = cvPoint2D32f(width / 2, height / 2);
	CvMat map_matrix2 = map_matrix;
	cv2DRotationMatrix(center, degree, 1.0, &map_matrix2);
	map[2] += (width_rotate - width) / 2;
	map[5] += (height_rotate - height) / 2;
	Mat img_rotate;
	
	warpAffine(img, img_rotate, map_matrix, Size(width_rotate, height_rotate), INTER_LINEAR, BORDER_CONSTANT, Scalar(255, 255, 255));
	return img_rotate;
}

Mat rotate_construct(Mat& ann, Mat& sim_nn, Mat& Refer)
{
	int  h = ann.rows;
	int  w = ann.cols;
	int  patch_size=80;
	Mat result(h, w, CV_64FC3, Scalar::all(0));
	Mat map(h, w, CV_64FC3, Scalar::all(0));
	//cout<<"h:"<<h<<"w"<<w<<endl;
	ofstream outfile;
	outfile.open("data.txt", ios::binary | ios::app | ios::in | ios::out);
	for (int i = 0; i <= h-patch_size; i=i+50)
	{		
		for (int j = 0; j <= w-patch_size; j=j+50)
		{
			//cout<<"j:"<<j<<endl;
			int start_row = ann.at<Vec3f>(i,j)[1];
			int start_col = ann.at<Vec3f>(i,j)[0];

			
			outfile<<"i:"<<i<<"  "<<"start_row: "<<start_row<<"\n";
			outfile<<"j:"<<j<<"  "<<"start_col: "<<start_col<<"\n";
			outfile<<"\n";
			

			double degree = sim_nn.at<Vec3f>(i,j)[0];//弧度
			double scale = sim_nn.at<Vec3f>(i,j)[1];
			//cout<<degree<<"::"<<scale<<endl;
			int rotate_h = patch_size*scale;
			int rotate_w = patch_size*scale;

			// cout<<"i:"<<i<<"  "<<"start_row: "<<start_row<<endl;
			// cout<<"j:"<<j<<"  "<<"start_col: "<<start_col<<endl;
			// cout<<"degree:"<<degree<<endl;
			// cout<<"scale:"<<scale<<endl;


			Mat mask_scale(rotate_h, rotate_w, CV_8UC3, Scalar::all(0));
			Mat mask(patch_size,patch_size, CV_8UC3, Scalar::all(0));
			Mat one(patch_size, patch_size, CV_64FC3, Scalar::all(0));
			Mat one1(70, 70, CV_64FC3, Scalar::all(1));
			double a = sin(degree), b = cos(degree);
			degree = degree * 180 / CV_PI;//数字
			start_col = start_col - rotate_h * a;
			if(start_col<0)
			{
				start_col=0;
			}
			
			Mat rotate_mask = rotateImage1(mask_scale, -degree);
			
			for (int ii = 0; ii <rotate_mask.rows; ii++)
			{
				for (int jj = 0; jj <rotate_mask.cols; jj++)
				{
					if (rotate_mask.at<Vec3b>(ii, jj) != Vec3b(255, 255, 255))
					{
						rotate_mask.at<Vec3b>(ii, jj) = Refer.at<Vec3b>(start_row + ii, start_col + jj);
					}
				}
			}
			
			Mat m_recover = rotateImage1(rotate_mask,degree);

			resize(m_recover,m_recover,Size(),1/scale,1/scale);
			
			int flag = 0;
			int i1 = m_recover.rows/2; 
			int j1 = m_recover.cols/2;
			m_recover(Range(i1-35, i1 + 35), Range(j1 -35, j1 +35)).copyTo(mask(Range(5, 75), Range(5,75)));
			one1.copyTo(one(Range(5, 75), Range(5, 75)));
			// imshow("mask",mask);
			// waitKey(0);
			mask.convertTo(mask, CV_64FC3, 1, 0);
			result(Range(i, i + patch_size), Range(j, j + patch_size)) += mask;
			map(Range(i, i + patch_size), Range(j, j + patch_size)) += one;
			//cout<<endl<<"finish1:"<<endl;
		}
		
	}
	outfile.close();
	for (int i = 1; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			if (map.at<Vec3d>(i, j) == Vec3d(0, 0, 0))
			{
				map.at<Vec3d>(i, j) = Vec3d(1, 1, 1);
			}
		}
	}
	//map.convertTo(map, CV_64FC3, 1, 0);
	divide(result, map, result);
	result.convertTo(result,CV_8UC3,1,0);
	imshow("result",result);
	return result;
	
}


int main()
{	
	// Mat A11 = imread("bridge_image112.jpg", 1);
	// Mat B11 = imread("sharp_rotate_image.jpg", 1);
	
	Mat A11 = imread("book_blur117.jpg", 1);
    Mat B11 = imread("book_sharp75.jpg", 1);

    // Mat A11 = imread("book_image82.jpg", 1);
    // Mat B11 = imread("book_image79.jpg", 1);
    //Mat B11 = imread("sharp_rotate_image.jpg", 1);

	//resize(A11,A11,Size(),0.5,0.5);
	//resize(B11,B11,Size(),0.5,0.5);
	//resize(A1,A1,Size(),0.8,0.8);
	//imshow("src",A11);
	//imshow("ref",B11);
	//imwrite("result/src_blur.jpg",A11);
	//imwrite("result/ref_sharp.jpg",B11);
	cout<<"A1_size::"<<A11.size()<<endl;
	cout<<"B1_size::"<<B11.size()<<endl;
	
	ofstream outfile("data.txt",ios::trunc);
	outfile.close();

	Mat re_image(A11.size(),CV_64FC3,Scalar::all(0));
	Mat re_image_mask(A11.size(),CV_8UC3,Scalar::all(0));
	double time0=static_cast<double>(getTickCount());
	int A_patchsize=200;
	int expend_size=300;
	for(int i=0;i<=(A11.rows-A_patchsize);i=i+100)
	{
		//cout<<"i="<<i<<endl;
		for(int j=0;j<=(A11.cols-A_patchsize);j=j+100)
		{
			
			Mat A1 = A11(Range(i,i+A_patchsize),Range(j,j+A_patchsize));
			//cout<<"..."<<endl;
			int start_row= (i-expend_size)>0?(i-expend_size):0;	 	
			int start_col= (j-expend_size)>0?(j-expend_size):0;   

			int end_row=0;
			int end_col=0;

			if(i+A_patchsize+expend_size<B11.rows)
			{
				end_row = i+A_patchsize+expend_size;
			}
			else if(i+A_patchsize+expend_size>= B11.rows)
			{
				end_row = B11.rows-1;
				start_row = B11.rows-A_patchsize-expend_size;
			}
			if(j+A_patchsize+expend_size<B11.cols)
			{
				end_col = j+A_patchsize+expend_size;
			}
			else if(j+A_patchsize+expend_size >= B11.cols)
			{
				end_col = B11.cols-1;
				start_col=B11.cols-A_patchsize-expend_size;
			}


			Mat B1 = B11(Range(start_row,end_row),Range(start_col,end_col));
			//cout<<start_row<<":"<<end_row<<endl;
			//cout<<start_col<<":"<<end_col<<endl;
			//cout<<"..."<<endl;	
			// Mat A1=imread("bulr_book1.jpg",1);
			// Mat B1 = imread("sharp_book1.jpg",1);
			int h = A1.rows;  
			int w = A1.cols;
			Mat ann,sim_ann;
			double time1=static_cast<double>(getTickCount());
	

			creat_ann(A1,B1,ann,sim_ann, h, w);
			cout<<"................"<<endl;
			//cout<<"ann:"<<ann(Range(0,3),Range(0,3))<<endl;
			
			//cout<<ann.at<Vec3f>(0,0)<<endl;
			//cout<<sim_ann<<endl;
			Mat result = rotate_construct(ann, sim_ann,B1);

			// imshow("result",result);
			// waitKey(0);
			time1=((double)getTickCount()-time1)/getTickFrequency();

			cout<<"一次运行时间为："<<time1<<"秒"<<endl;
			
			
			Mat result_mask(result.size(),CV_8UC3,Scalar::all(0));

			for(int ii=0;ii<result.rows;ii++)
			{
				for(int jj=0;jj<result.cols;jj++)
				{
					if(result.at<Vec3b>(ii,jj)!=Vec3b(0,0,0))
					{
						result_mask.at<Vec3b>(ii,jj)=Vec3b(1,1,1);
					}
				}
			}
			result.convertTo(result,CV_64FC3,1,0);
			re_image(Range(i,i+A_patchsize),Range(j,j+A_patchsize))+=result;
			re_image_mask(Range(i,i+A_patchsize),Range(j,j+A_patchsize))+=result_mask;
			imshow("A1",A1);
			imshow("B1",B1);
			//imshow("result",result);
			//cout<<"result.size:"<<result(Range(0,8),Range(0,8))<<endl;
			waitKey(0);
			cout<<"::"<<A1.size()<<"::"<<B1.size()<<endl;	
			
		}
	}
	cout<<"finished1"<<endl;
	for(int ii=0;ii<re_image_mask.rows;ii++)
	{
		for(int jj=0;jj<re_image_mask.cols;jj++)
		{
			if(re_image_mask.at<Vec3b>(ii,jj)==Vec3b(0,0,0))
			{
				re_image_mask.at<Vec3b>(ii,jj)=Vec3b(1,1,1);
			}
		}
	}

	re_image_mask.convertTo(re_image_mask,CV_64FC3,1,0);

	Mat result;
	divide(re_image,re_image_mask,result);

	time0=((double)getTickCount()-time0)/getTickFrequency();
	cout<<"一次运行时间为："<<time0<<"秒"<<endl;
	result.convertTo(result,CV_8UC3,1,0);
	imshow("result.jpg",result);
	
	
	imwrite("result/LL_building_patch_80_rotate_scale_result_100_10_0.95_60_iter80.jpg",result);
	waitKey(0);

	
}


