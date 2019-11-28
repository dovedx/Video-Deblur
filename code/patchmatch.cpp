#include "patchmatch.h"



#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

int xform_cos_table[NUM_ANGLES];
  int xform_sin_table[NUM_ANGLES];
int xform_scale_table[NUM_SCALES];


template<int TPATCH_W, int DO_BRANCH>
int sim_patch_dist_ab(Params *p, BITMAP *a, int ax, int ay, BITMAP *b, int bx, int by, int bs, int bt, int maxval)
{
  int adata[TPATCH_W*TPATCH_W];
  int *ptr = adata;
  for (int dy = 0; dy < TPATCH_W; dy++) 
  {
    int *row = ((int *) a->line[ay+dy])+ax;
    for (int dx = 0; dx < TPATCH_W; dx++) 
    {
      *ptr++ = *row++;
    }
  }
  XFORM bpos = get_xform(p, bx, by, bs, bt);
  return sim_fast_patch_dist<TPATCH_W, DO_BRANCH>(adata, b, bpos, maxval);
}

void check_offset(Params *p, BITMAP *b, int x, int y, int xp, int yp) 
{
	//int h = p->patch_w/2;
	//if ((unsigned) (xp-h) >= (unsigned) BEW ||
	//    (unsigned) (yp-h) >= (unsigned) BEH) { fprintf(stderr, "offset (%d, %d) => (%d, %d) out of range b: %dx%d\n", x, y, xp, yp, b->w, b->h); exit(1); }
	if ((unsigned)xp >= (unsigned)BEW ||(unsigned)yp >= (unsigned)BEH) 
  {
		fprintf(stderr, "offset (%d, %d) => (%d, %d) out of range b: %dx%d\n", x, y, xp, yp, b->w, b->h); exit(1);
	}
}

template<int PATCH_W>
void sim_attempt_n(int &err, int &xbest, int &ybest, int &sbest, int &tbest, int *adata, BITMAP *b, XFORM bpos, int bx, int by, int bs, int bt, Params *p) {
	//int h = PATCH_W/2;
	if ((bx != xbest || by != ybest || bs != sbest || bt != tbest) &&(unsigned)(bx) < (unsigned)(b->w - PATCH_W + 1) &&(unsigned)(by) < (unsigned)(b->h - PATCH_W + 1))
	{
		//XFORM bpos = get_xform(p, bx, by, bs, bt);
		int current = sim_fast_patch_dist<PATCH_W, 1>(adata, b, bpos, err);
		if (current < err)
    {
			err = current;
			xbest = bx;
			ybest = by;
			sbest = bs;
			tbest = bt;
		}
	}
}

template<int PATCH_W>
BITMAP *sim_init_dist_n(Params *p, BITMAP *a, BITMAP *b, BITMAP *ann, BITMAP *ann_sim) 
{
 // init_xform_tables();
  BITMAP *ans = create_bitmap(a->w, a->h);
  clear_to_color(ans, 0);
  for (int y = 0; y < AEH; y++) 
  {
    int *row = (int *) ans->line[y];
    for (int x = 0; x < AEW; x++)
    {
      int bx, by, bs, bt;
      getnn(ann, x, y, bx, by);
      getnn1(ann_sim, x, y, bs, bt);
      row[x] = sim_patch_dist_ab<PATCH_W, 0>(p, a, x, y, b, bx, by, bs, bt, INT_MAX);
    }
  }
  return ans;
}


void destroy_bitmap(BITMAP *bmp) {
	if (bmp) {
		delete[] bmp->line;
		delete[] bmp->data;
		delete bmp;
	}
}
void destroy_region_masks(RegionMasks *m)
{
  if (!m) 
  { 
  	return; 
  }
  destroy_bitmap(m->bmp);
  delete m;
}


Box get_abox(Params *p, BITMAP *a, RegionMasks *amask, int trim_patch)
{
	if (!amask) 
	{
    Box ans;
    ans.xmin = ans.ymin = 0;
		ans.xmax = trim_patch ? (a->w - p->patch_w + 1) : a->w; 
		ans.ymax = trim_patch ? (a->h - p->patch_w + 1) : a->h;
    return ans;
  }
  Box ans = amask->box[0];
  //save_bitmap("amask.bmp", amask->bmp, NULL);
  if (ans.xmin < 0 || ans.ymin < 0 || ans.xmax > a->w-p->patch_w+1 || ans.ymax > a->h-p->patch_w+1) { fprintf(stderr, "box out of range %d %d %d %d (%d %d %d %d)\n", ans.xmin, ans.ymin, ans.xmax, ans.ymax, 0, 0, a->w-p->patch_w+1, a->h-p->patch_w+1); exit(1); }
  if (ans.xmin >= ans.xmax || ans.ymin >= ans.ymax) { ans.xmin = ans.ymin = 0; ans.xmax = ans.ymax = 1; } //fprintf(stderr, "box size 0 (%d %d %d %d)\n", ans.xmin, ans.ymin, ans.xmax, ans.ymax); exit(1); }
  // FIXME: Instead, set box size to 1 at (0,0) if it has size zero
  printf("get_abox (%dx%d) => %d %d %d %d\n", a->w, a->h, ans.xmin, ans.ymin, ans.xmax, ans.ymax);
  return ans;
}


RegionMasks::RegionMasks(Params *p, BITMAP *region_masks, int full, BITMAP *bmask)
{
  bmp = region_masks;
  if (!bmp) { fprintf(stderr, "Region mask is NULL\n"); exit(1); }
  for (int i = 0; i < 256; i++) 
  {
    box[i].xmin = box[i].ymin = INT_MAX;
    box[i].xmax = box[i].ymax = -INT_MAX;
  }
  if (full) 
  {
    box[0].xmin = box[0].ymin = 0;
    box[0].xmax = region_masks->w - p->patch_w + 1;
    box[0].ymax = region_masks->h - p->patch_w + 1;
  } 
  else 
  {
    for (int y = 0; y <= region_masks->h - p->patch_w; y++) 
    {
      int *row = (int *) region_masks->line[y];
      for (int x = 0; x <= region_masks->w - p->patch_w; x++) 
      {
        int i = row[x]&255;
        //if (bmask && ((int *) bmask->line[y])[x]) { continue; }
        //if ((unsigned) i > (unsigned) 255) { fprintf(stderr, "Error: region mask index not a uint8 (%d at %d,%d).\n", i, x, y); exit(1); }
        if (x < box[i].xmin) { box[i].xmin = x; }
        if (x >= box[i].xmax) { box[i].xmax = x+1; }
        if (y < box[i].ymin) { box[i].ymin = y; }
        if (y >= box[i].ymax) { box[i].ymax = y+1; }
      }
    }
  }
  for (int i = 0; i < 256; i++) 
  {
    if (box[i].xmin != INT_MAX) 
    {
      printf("%d => %d %d %d %d\n", i, box[i].xmin, box[i].ymin, box[i].xmax, box[i].ymax);
    }
  }
}


BITMAP *create_bitmap(int w, int h)
{
	BITMAP *ans = new BITMAP(w, h);
  ans->data = new unsigned char[4*w*h]; // always 32bit
  ans->line = new unsigned char*[h];
  for (int y = 0; y < h; y++) 
    ans->line[y] = &ans->data[y*4*w];
  return ans;
}

void init_params(Params *p)
{
	init_openmp(p);
}

void init_openmp(Params *p)
{
 int last_cores = 0;
 if (p->cores != last_cores) 
 {
      last_cores = p->cores;
  #if USE_OPENMP
      omp_set_num_threads(p->cores);
      omp_set_nested(1);
      omp_set_dynamic(0);
  #endif
    }
}

void clear(BITMAP *bmp)
{
	clear_to_color(bmp, 0);
}
void clear_to_color(BITMAP *bmp, int c)
{
	for (int y = 0; y < bmp->h; y++) 
	{
    int *row = (int *) bmp->line[y];
    for (int x = 0; x < bmp->w; x++) 
    {		
      row[x] = c;
    }
  }
}

int fixmul(int a0, int b0)
{
  long long a = a0;
  long long b = b0;
  return (int) ((a*b)>>16);
}
int xform_tables_init = 0;
void init_xform_tables(double SCALE_MIN, double SCALE_MAX, int force_init)
{
  // if (xform_tables_init && !force_init) 
  // { 
  // 	return; 
  // }
  // xform_tables_init = 1;
  for (int i = 0; i < NUM_ANGLES; i++) 
  {
    // double theta = i*(2.0*M_PI/NUM_ANGLES);
    // xform_cos_table[i] = int(cos(theta)*65536.0);
    // xform_sin_table[i] = int(sin(theta)*65536.0);
    // 
     double theta = i*(((10.0)*2.0*M_PI/360.0)/NUM_ANGLES); 
     xform_cos_table[i] = int(cos(2*M_PI-theta )*65536.0);
     xform_sin_table[i] = int(sin(2*M_PI-theta)*65536.0);

  }
  double a = (0.8), b = (1.2);
  for (int i = 0; i < NUM_SCALES; i++) 
  {
    double scale = (a+(b-a)*i*1.0/NUM_SCALES);
    xform_scale_table[i] = int(scale*65536.0);
  }
}

void getnn(BITMAP *ann, int x, int y, int &xp, int &yp)
{
	int dest = ((int *) ann->line[y])[x];
	xp = INT_TO_X(dest);	
	yp = INT_TO_Y(dest);
}

void getnn1(BITMAP *ann, int x, int y, int &xp, int &yp)
{
  int dest = ((int *) ann->line[y])[x];
  xp = INT_TO_X(dest);  
  yp = INT_TO_Y(dest);
}
XFORM get_xform(Params *p, int x, int y, int scale, int theta)  /* x and y not left shifted. */
{
  int c = xform_cos_table[theta];
  int s = xform_sin_table[theta];
  if ((unsigned) scale >= (unsigned) NUM_SCALES) { fprintf(stderr, "scale out of range: %d (%d)\n", scale, NUM_SCALES); exit(1); }
  int scalef = xform_scale_table[scale];
  c = fixmul(c, scalef);
  s = fixmul(s, scalef);
  XFORM ans;
  ans.dxdu = c;
  ans.dydu = -s;
  ans.dxdv = s;
  ans.dydv = c;
  int h = p->patch_w/2;
  int xc = x+h, yc = y+h;
  ans.x0 = (xc<<16)-ans.dxdu*h-ans.dxdv*h;
  ans.y0 = (yc<<16)-ans.dydu*h-ans.dydv*h;
  return ans;
}


BITMAP *init_nn(Params *p, BITMAP *a, BITMAP *b, BITMAP *bmask, RegionMasks *region_masks, RegionMasks *amask, int trim_patch, BITMAP *ann_window, BITMAP *awinsize)
{
	BITMAP *bmp = create_bitmap(a->w, a->h);
  clear(bmp);
  //cout<<"p->patch_w :"<<p->patch_w <<endl;
	int ew = trim_patch ? (b->w - p->patch_w + 1) : b->w;
	int eh = trim_patch ? (b->h - p->patch_w + 1) : b->h;
 
 	printf("init_nn: ew=%d, eh=%d\n", ew, eh);

 	Box box = get_abox(p, a, amask);

 
    #pragma omp parallel for schedule(static, 8)
    for (int y = box.ymin; y < box.ymax; y++) 
    {
      //unsigned int seed = rand();
      int *row = (int *) bmp->line[y];
      int *arow = amask ? (int *) amask->bmp->line[y]: NULL;
      for (int x = box.xmin; x < box.xmax; x++) 
      {
        if (amask && arow[x]) 
        { 
        	continue; 
        }
        int xp = 0, yp = 0;
        xp = rand()%ew;
        yp = rand()%eh;
        //cout<<ew<<"::"<<eh<<":"<<xp<<":"<<yp<<endl;
        // xp = seed % ew;
        // seed = RANDI(seed);
        // yp = seed % eh;
        // seed = RANDI(seed);
        //if (iter >= MAX_NN_GUESS_ITERS) { fprintf(stderr, "Warning: too many iters at %d,%d\n", x, y); }
        //if (x == 1 && y == 1) { printf("1, 1 => %d %d\n", xp, yp); }
        row[x] = XY_TO_INT(xp, yp);

      }
    }

  return bmp;
}
// BITMAP *init_dist(Params *p, BITMAP *a, BITMAP *b, BITMAP *ann, BITMAP *bmask=NULL, RegionMasks *region_masks=NULL, RegionMasks *amask=NULL)
// {

// }

BITMAP *sim_init_nn(Params *p, BITMAP *a, BITMAP *b, BITMAP *&ann_sim)
{
 // init_xform_tables();
  BITMAP *ann = init_nn(p, a, b);
  ann_sim = create_bitmap(ann->w, ann->h);
  clear(ann_sim);
  //int h = p->patch_w/2;
  for (int y = 0; y < AEH; y++) 
  {
    //int *ann_row = (int *) ann->line[y];
    int *row = (int *) ann_sim->line[y];
    for (int x = 0; x < AEW; x++) 
    {
      //int xp, yp;
      //getnn(ann, x, y, xp, yp);
      //xp += h;
      //yp += h;
      //ann_row[x] = XY_TO_INT(xp, yp);
      row[x] = XY_TO_INT(rand()%NUM_SCALES, rand()%NUM_ANGLES);
      //check_offset(p, b, x, y, xp, yp);
      //cout<<rand()%NUM_SCALES<<"...."<<rand()%NUM_ANGLES<<endl;
    }
  }
  return ann;
}

BITMAP *sim_init_dist(Params *p, BITMAP *a, BITMAP *b, BITMAP *ann, BITMAP *ann_sim)
{
	return sim_init_dist_n<80>(p, a, b, ann, ann_sim); 
}


template<int PATCH_W>
void sim_nn_n(Params *p, BITMAP *a, BITMAP *b,
            BITMAP *ann, BITMAP *ann_sim, BITMAP *annd,
            RegionMasks *amask=NULL, BITMAP *bmask=NULL,
            int level=0, int em_iter=0, RecomposeParams *rp=NULL, int offset_iter=0, int update_type=0, int cache_b=0,
            RegionMasks *region_masks=NULL, int tiles=-1) 
{
 // init_xform_tables();
  if (tiles < 0) { tiles = p->cores; }
  printf("in sim_nn_n, tiles=%d, rs_max=%d\n", tiles, p->rs_max);
  Box box = get_abox(p, a, amask);
  int nn_iter = 0;
  for (; nn_iter < 5; nn_iter++) 
  {
    cout<<"nn_iter"<<nn_iter<<endl;
    unsigned int iter_seed = rand();

    //#pragma omp parallel num_threads(tiles)
    //{
#if SYNC_WRITEBACK
      int *ann_writeback = new int[a->w];
      int *annd_writeback = new int[a->w];
      int *ann_sim_writeback = new int[a->w];
#endif
#if USE_OPENMP
      int ithread = omp_get_thread_num();
#elses
      int ithread = 0;
#endif
      int xmin = box.xmin, xmax = box.xmax;
     // int ymin = box.ymin + (box.ymax-box.ymin)*ithread/tiles;
     // int ymax = box.ymin + (box.ymax-box.ymin)*(ithread+1)/tiles;
     int ymin = box.ymin;
     int ymax = box.ymax;

      int ystart = ymin, yfinal = ymax, ychange=1;
      int xstart = xmin, xfinal = xmax, xchange=1;
      if ((nn_iter + offset_iter) % 2 == 1) 
      {
        xstart = xmax-1; xfinal = xmin-1; xchange=-1;
        ystart = ymax-1; yfinal = ymin-1; ychange=-1;
      }
      int dx = -xchange, dy = -ychange; 

      int bew = b->w-PATCH_W, beh = b->h-PATCH_W;
      int max_mag = max(b->w, b->h);
      int rs_ipart = int(p->rs_iters);
      double rs_fpart = p->rs_iters - rs_ipart;
      int rs_max = p->rs_max;
      if (rs_max > max_mag) { rs_max = max_mag; }

      int adata[PATCH_W*PATCH_W];
      for (int y = ystart; y != yfinal; y += ychange) 
      {
        int *annd_row = (int *) annd->line[y];
        for (int x = xstart; x != xfinal; x += xchange) 
        {
          for (int dy0 = 0; dy0 < PATCH_W; dy0++) 
          {
            int *drow = ((int *) a->line[y+dy0])+x;
            int *adata_row = adata+(dy0*PATCH_W);
            for (int dx0 = 0; dx0 < PATCH_W; dx0++) 
            {
              adata_row[dx0] = drow[dx0];
            }
          }
          
          int xbest, ybest, sbest, tbest;
          getnn(ann, x, y, xbest, ybest);
          getnn1(ann_sim, x, y, sbest, tbest);
          check_offset(p, b, x, y, xbest, ybest);

          int err = annd_row[x];
          if (err == 0) 
          {
#if SYNC_WRITEBACK
            if (y+ychange == yfinal) 
            {
              ann_writeback[x] = XY_TO_INT(xbest, ybest);
              ann_sim_writeback[x] = XY_TO_INT(sbest, tbest);
              annd_writeback[x] = err;
            }
#endif      
            continue;
          }

          /* Propagate */
          if (p->do_propagate) 
          {
            /* Propagate x */
            if ((unsigned) (x+dx) <(unsigned) (ann->w-PATCH_W))
            {
              //cout<<"pppppppppppppp"<<endl;
              int xpp, ypp, spp, tpp;
              getnn(ann, x+dx, y, xpp, ypp);
              getnn1(ann_sim, x+dx, y, spp, tpp);
              //xpp -= dx;
              XFORM bpos = get_xform(p, xpp, ypp, spp, tpp);
              xpp -= (bpos.dxdu*dx)>>16;
              ypp -= (bpos.dydu*dx)>>16;
              bpos = get_xform(p, xpp, ypp, spp, tpp);

              sim_attempt_n<PATCH_W>(err, xbest, ybest, sbest, tbest, adata, b, bpos, xpp, ypp, spp, tpp, p);
              check_offset(p, b, x, y, xbest, ybest);
            }

            /* Propagate y */
            if ((unsigned) (y+dy) <(unsigned) (ann->h-PATCH_W)) 
            {
               //cout<<"pppppppppp"<<endl;
              int xpp, ypp, spp, tpp;
              getnn(ann, x, y+dy, xpp, ypp);
              getnn1(ann_sim, x, y+dy, spp, tpp);
              //xpp -= dx;
              XFORM bpos = get_xform(p, xpp, ypp, spp, tpp);
              xpp -= (bpos.dxdv*dy)>>16;
              ypp -= (bpos.dydv*dy)>>16;
              bpos = get_xform(p, xpp, ypp, spp, tpp);

              sim_attempt_n<PATCH_W>(err, xbest, ybest, sbest, tbest, adata, b, bpos, xpp, ypp, spp, tpp, p);
              check_offset(p, b, x, y, xbest, ybest);
            }
          }    
        
          /* Random search */
          unsigned int seed = (x | (y<<11)) ^ iter_seed;
          seed = RANDI(seed);
          int rs_iters = 1-(seed*(1.0/(RAND_MAX-1))) < rs_fpart ? rs_ipart + 1: rs_ipart;
  //      int rs_iters = 1-random() < rs_fpart ? rs_ipart + 1: rs_ipart;
  //      
         // cout<<"rs_iters:"<<rs_iters<<endl;     

          int rs_max_curr = rs_max;

          // int h = p->patch_w/2;
          // int ymin_clamp = h, xmin_clamp = h;
          // int ymax_clamp = BEH+h, xmax_clamp = BEW+h;

          for (int mag = rs_max_curr; mag >= p->rs_min; mag = int(mag*0.5)) 
          {
          //  cout<<"rs_max_curr"<<rs_max_curr<<endl;
            int smag = NUM_SCALES*mag/rs_max_curr;
            int tmag = (NUM_SCALES == NUM_ANGLES) ? smag: (NUM_ANGLES*mag/rs_max_curr);  // FIXME: This should be divided by 2
            for (int rs_iter = 0; rs_iter < 1; rs_iter++)
            {
              int ymin = max(ybest-mag,0), ymax = min(ybest+mag+1,beh);
              int xmin = max(xbest-mag,0), xmax = min(xbest+mag+1,bew);
              // if (xmin < xmin_clamp) { xmin = xmin_clamp; }
              // if (ymin < ymin_clamp) { ymin = ymin_clamp; }
              // if (xmax > xmax_clamp) { xmax = xmax_clamp; }
              // if (ymax > ymax_clamp) { ymax = ymax_clamp; }

              int smin = sbest-smag, smax = sbest+smag+1; 
              int tmin = tbest-tmag, tmax = tbest+tmag+1;
              if (smin < 0) { smin = 0; } 
              if (smax > NUM_SCALES) { smax = NUM_SCALES; }
              if (tmin < 0) { tmin = 0; } 
              if (tmax > NUM_SCALES) { tmax = NUM_SCALES; }
              //fprintf(stderr, "RS: xbest: %d, ybest: %d, err: %d, mag: %d, bew: %d, beh: %d, smag: %d, tmag: %d, xmin: %d, xmax: %d, ymin: %d, ymax: %d, smin: %d, smax: %d, tmin: %d, tmax: %d\n", xbest, ybest, err, mag, bew, beh, smag, tmag, xmin, xmax, ymin, ymax, smin, smax, tmin, tmax); fflush(stderr);
              // cout<<smin<<"--"<<smax<<endl;
              // cout<<tmin<<"---"<<tmax<<endl;

              //seed = RANDI(seed);
              int xpp = xmin+rand()%(xmax-xmin);
             // seed = RANDI(seed);
              int ypp = ymin+rand()%(ymax-ymin);
             // seed = RANDI(seed);
              int spp = smin+rand()%(smax-smin);
             // seed = RANDI(seed);
              int tpp = tmin+rand()%(tmax-tmin);
             
          
              XFORM bpos = get_xform(p, xpp, ypp, spp, tpp);
              sim_attempt_n<PATCH_W>(err, xbest, ybest, sbest, tbest, adata, b, bpos, xpp, ypp, spp, tpp, p);
              //check_offset(p, b, x, y, xbest, ybest);
            }
          }

#if SYNC_WRITEBACK
          if (y+ychange != yfinal) 
          {     
#endif
          ((int *) ann->line[y])[x] = XY_TO_INT(xbest, ybest);
          ((int *) ann_sim->line[y])[x] = XY_TO_INT(sbest, tbest);
          ((int *) annd->line[y])[x] = err;
#if SYNC_WRITEBACK
          } 
          else 
          {
            ann_writeback[x] = XY_TO_INT(xbest, ybest);
            annd_writeback[x] = err;
            ann_sim_writeback[x] = XY_TO_INT(sbest, tbest);
          }
#endif
        } // x
      } // y

#if SYNC_WRITEBACK
      #pragma omp barrier
      int ywrite = yfinal-ychange;
      if (ymin < ymax && (unsigned) ywrite < (unsigned) AEH) 
      {
        int *ann_line = (int *) ann->line[ywrite];
        int *annd_line = (int *) annd->line[ywrite];
        int *ann_sim_line = (int *) ann_sim->line[ywrite];
        for (int x = xmin; x < xmax; x++) 
        {
          ann_line[x] = ann_writeback[x];
          annd_line[x] = annd_writeback[x];
          ann_sim_line[x] = ann_sim_writeback[x];
        }
      }
      delete[] ann_writeback;
      delete[] annd_writeback;
      delete[] ann_sim_writeback;
#endif
    //} // parallel
   // fprintf(stderr, "done with %d iters\n", nn_iter);
  } // nn_iter

  printf("done sim_nn_n, %d iters, rs_max=%d\n", nn_iter, p->rs_max);
  //cout<<"fff"<< endl;
}


void sim_nn(Params *p, BITMAP *a, BITMAP *b,
            BITMAP *ann, BITMAP *ann_sim, BITMAP *annd,
            RegionMasks *amask, BITMAP *bmask,int level, int em_iter, RecomposeParams *rp, 
            int offset_iter, int update_type, int cache_b,
            RegionMasks *region_masks, int tiles)
{
	sim_nn_n<80>(p, a, b, ann, ann_sim, annd, amask, bmask, level, em_iter, rp, offset_iter, update_type, cache_b, region_masks, tiles);

}
