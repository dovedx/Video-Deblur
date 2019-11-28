#ifndef variable_h_
#define variable_h_

#define AEW (a->w - p->patch_w + 1)
#define AEH (a->h - p->patch_w + 1)
#define BEW (b->w - p->patch_w + 1)
#define BEH (b->h - p->patch_w + 1)

#define XY_TO_INT(x, y) (((y)<<12)|(x))
#define INT_TO_X(v) ((v)&((1<<12)-1))
#define INT_TO_Y(v) ((v)>>12)
#define XY_TO_INT_SHIFT 12
	


#endif	