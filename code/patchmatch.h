#ifndef patchmatch_h_
#define patchmatch_h_

#include "variable.h"
#include <stdio.h>
#include <vector>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <iostream>
using namespace std;
using namespace cv;

#define ANGLE_SHIFT XY_TO_INT_SHIFT
#define NUM_ANGLES (1<<ANGLE_SHIFT)
#define SCALE_SHIFT XY_TO_INT_SHIFT
#define NUM_SCALES (1<<SCALE_SHIFT)
//#define RAND_MAX 4294967295U

#define RANDI(u) (18000 * ((u) & 65535) + ((u) >> 16))

#define SAMPLE_NN             1
#define SAMPLE_BILINEAR_EXACT 0

struct BITMAP 
{
  int w, h;
  unsigned char **line;
  unsigned char *data;
  BITMAP(int ww = 1, int hh = 1): w(ww), h(hh) {}
};


class Params 
{ 
public:
  int algo;              /* Algorithm to use, one of ALGO_*. */
  
  /* Randomized NN algorithm parameters. */
  int patch_w;           /* Width and height of square patch. */
  int vec_len;           /* Length of vector if using vectorized NN algorithms (vecnn.h), for non-square patches and feature descriptors. */
  int nn_iters;          /* Iters of randomized NN algorithm. */
  int rs_max;            /* Maximum width for RS. */
  int rs_min;            /* Minimum width for RS. */
  double rs_ratio;       /* Ratio (< 1) of successive RS sizes. */
  double rs_iters;       /* RS iters per pixel.  1 => search all scales once. */
  int do_propagate;      /* Do propagation. */
  int gpu_prop;          /* Maximum propagation distance for GPU algo. */
  int xmin, ymin;        /* Min update region coord, or -1 for no box. */
  int xmax, ymax;        /* Max update region coord, or -1 for no box. */
  int resample_seamcarv; /* Resample via seam carving. */
  int vote_algo;         /* Vote algorithm, one of VOTE_*. */
  int prefer_coherent;   /* Prefer coherent regions, bool, default false. */
  int allow_coherent;    /* This must be enabled for the previous flag to take effect. */
  int cores;             /* If > 1, use OpenMP. */
  int window_w;          /* Constraint search window width. */
  int window_h;          /* Constraint search window height. */
  int weight_r;          /* Multiplicative weights for R, G, B in distance computation. */
  int weight_g;
  int weight_b;
  
  /* Tree NN algorithm parameters. */
  int pca_dim;           /* Tree PCA dim, INT_MAX for no PCA. */
  double pca_var;        /* Fraction of total variance, e.g. 0.95 for 95%, negative to not use this param. */
  float eps;             /* Tree epsilon. */
  int kcoherence_k;      /* k for kcoherence. */
  int kcoherence_iters;  /* Iters of kcoherence "propagation", 2 iters is Lefebre '95. */
  int kcoherence_neighbors; /* k-coherence neighbors, 4 or 8. */

  int knn;
  int knn_algo;
  int restarts;
  int enrich_iters;
  int enrich_times;
  int do_inverse_enrich;
  int do_enrich;

  /* Defaults. */
  Params()
    :algo(0),
     patch_w(80),
     vec_len(0),
     nn_iters(5),
     rs_max(INT_MAX),
     rs_min(1),
     rs_ratio(0.5),
     rs_iters(1),
     do_propagate(1),
     gpu_prop(8),
     xmin(-1), ymin(-1),
     xmax(-1), ymax(-1),
     resample_seamcarv(0),
     vote_algo(0),
     pca_dim(25),
     pca_var(-1),
     eps(2),
     prefer_coherent(0),
     allow_coherent(0),
     cores(1),
     window_w(INT_MAX),
     window_h(INT_MAX),
     weight_r(1),
     weight_g(1),
     weight_b(1),
     kcoherence_k(2),
     kcoherence_iters(2),
     kcoherence_neighbors(8),
     knn(0),
     knn_algo(0),
     restarts(1),
     enrich_iters(0),
     enrich_times(1),
     do_inverse_enrich(1),
     do_enrich(1)
     { }
};

class Box { 
public:
  int xmin, ymin, xmax, ymax;
};

class RegionMasks { 
public:
  BITMAP *bmp;
  Box box[256];
  RegionMasks(Params *p, BITMAP *region_masks, int full=0, BITMAP *bmask=NULL);
};

class RecomposeParams 
{ public:
	  int minnn_optp_nn_iters;    /* Optimized params: NN iters for previous offsets. */
	  int minnn_optp_rs_max;      /* Optimized params: Max RS for previous offsets. */

	  RecomposeParams()
	    :minnn_optp_nn_iters(2),
	     minnn_optp_rs_max(1) {}
};

class XFORM { 
public:
  int x0, y0, dxdu, dydu, dxdv, dydv;   /* Coords left shifted by 16. */
};


void destroy_region_masks(RegionMasks *m);

Box get_abox(Params *p, BITMAP *a, RegionMasks *amask, int trim_patch=1);

BITMAP *create_bitmap(int w, int h);

void init_params(Params *p);

void init_openmp(Params *p);

void clear(BITMAP *bmp);
void clear_to_color(BITMAP *bmp, int c);

int fixmul(int a, int b);

void init_xform_tables(double SCALE_MIN=0.8, double SCALE_MAX=1.2, int force_init=0);

XFORM get_xform(Params *p, int x, int y, int scale, int theta);  /* x and y not left shifted. */
		
void getnn(BITMAP *ann, int x, int y, int &xp, int &yp);
void getnn1(BITMAP *ann, int x, int y, int &xp, int &yp);
XFORM get_xform(Params *p, int x, int y, int scale, int theta);  /* x and y not left shifted. */

BITMAP *init_nn(Params *p, BITMAP *a, BITMAP *b, BITMAP *bmask=NULL, RegionMasks *region_masks=NULL, RegionMasks *amask=NULL, int trim_patch=1, BITMAP *ann_window=NULL, BITMAP *awinsize=NULL);
//BITMAP *init_dist(Params *p, BITMAP *a, BITMAP *b, BITMAP *ann, BITMAP *bmask=NULL, RegionMasks *region_masks=NULL, RegionMasks *amask=NULL);

BITMAP *sim_init_nn(Params *p, BITMAP *a, BITMAP *b, BITMAP *&ann_sim);
BITMAP *sim_init_dist(Params *p, BITMAP *a, BITMAP *b, BITMAP *ann, BITMAP *ann_sim);

void sim_nn(Params *p, BITMAP *a, BITMAP *b,
            BITMAP *ann, BITMAP *ann_sim, BITMAP *annd,
            RegionMasks *amask=NULL, BITMAP *bmask=NULL,
            int level=0, int em_iter=0, RecomposeParams *rp=NULL, int offset_iter=0, int update_type=0, int cache_b=0,
            RegionMasks *region_masks=NULL, int tiles=-1);

template<int TPATCH_W, int DO_BRANCH>
int sim_fast_patch_dist(int *adata, BITMAP *b, XFORM bpos, int maxval) 
{
	/* Do bounds checking outside the inner loop */
	int ul_x = bpos.x0;
	int ul_y = bpos.y0;
	int ur_x = bpos.x0 + bpos.dxdu*(TPATCH_W - 1);
	int ur_y = bpos.y0 + bpos.dydu*(TPATCH_W - 1);
	int ll_x = bpos.x0 + bpos.dxdv*(TPATCH_W - 1);
	int ll_y = bpos.y0 + bpos.dydv*(TPATCH_W - 1);
	int lr_x = ll_x + bpos.dxdu*(TPATCH_W - 1);
	int lr_y = ll_y + bpos.dydu*(TPATCH_W - 1);
	int bw16 = (b->w - 1) << 16, bh16 = (b->h - 1) << 16;
	if ((unsigned)ul_x >= (unsigned)bw16 ||
		(unsigned)ul_y >= (unsigned)bh16 ||
		(unsigned)ur_x >= (unsigned)bw16 ||
		(unsigned)ur_y >= (unsigned)bh16 ||
		(unsigned)ll_x >= (unsigned)bw16 ||
		(unsigned)ll_y >= (unsigned)bh16 ||
		(unsigned)lr_x >= (unsigned)bw16 ||
		(unsigned)lr_y >= (unsigned)bh16) 
  {
		return INT_MAX - 4096;
	}

  	int ans = 0;
  	int bx_row = bpos.x0, by_row = bpos.y0;
    #if 1
    	bx_row += 32768;
    	by_row += 32768;
    #endif
	for (int dy = 0; dy < TPATCH_W; dy++) 
  {
		int bx = bx_row, by = by_row;     
		for (int dx = 0; dx < TPATCH_W; dx++) 
    {
			unsigned int c1 = adata[dx];
			int r2, g2, b2;

			int c2 = ((int *)b->line[(by) >> 16])[(bx) >> 16];
     // cout<<"by:"<<(by>>16)<<"bx:"<<(bx>>16)<<endl;

			r2 = (c2 & 255);
			g2 = (c2 >> 8) & 255;
			b2 = (c2 >> 16);
//#endif
			int dr = (c1 & 255) - r2;
			int dg = ((c1 >> 8) & 255) - g2;
			int db = (c1 >> 16) - b2;
      //cout<<dr<<":"<<dg<<":"<<db<<endl;
			ans += (dr*dr+dg*dg+db*db);
      // cout<<dr<<":"<<dg<<":"<<db<<endl;
      
			if (DO_BRANCH && ans > maxval) { return ans; }
  			bx += bpos.dxdu;
  			by += bpos.dydu;
		}
		adata += TPATCH_W;
   // cout<<"......."<<endl;
   // Mat ss(2,2,CV_8UC3,Scalar::all(0));
   // imshow("ss",ss);
   //    waitKey(0);
		bx_row += bpos.dxdv;
		by_row += bpos.dydv;
	}
   // cout<<"finished"<<endl;
  //cout<<"bpos.dxdu:"<< ((bpos.dxdu)>>16)<<"bpos.dydu:"<<((bpos.dydu)>>16)<<endl;
  //cout<<"bpos.dxdv:"<< ((bpos.dxdv)>>16)<<"bpos.dydv:"<<((bpos.dydv)>>16)<<endl;
 
  //cout<<"ans:"<<ans<<endl;
	return ans;
}


#endif