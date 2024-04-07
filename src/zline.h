{
 
 

 
#ifdef INTERP_RGB
    #define RGB(x) x
#endif
#ifdef INTERP_RGB 
    #define RGBPIXEL *pp = RGB_TO_PIXEL(r >> 8,g >> 8,b >> 8)
#endif


#ifndef INTERP_RGB /* INTERP_RGB */
    #define RGB(x)
#endif /* INTERP_RGB */
#ifndef INTERP_RGB  
    #define RGBPIXEL *pp = color
#endif

#ifdef INTERP_Z
#define ZZ(x) x
#define PUTPIXEL() 				\
  {						\
    zz=z >> ZB_POINT_Z_FRAC_BITS;		\
    if (ZCMP(zz,*pz))  { 			\
      RGBPIXEL;	\
      *pz=zz; 					\
    }						\
  }
#else /* INTERP_Z */
#define ZZ(x)
#define PUTPIXEL() RGBPIXEL
#endif /* INTERP_Z */

#define DRAWLINE(dx,dy,inc_1,inc_2) \
    n=dx;\
    ZZ(zinc=(p2->z-p1->z)/n);\
    RGB(rinc=((p2->r-p1->r) << 8)/n;\
        ginc=((p2->g-p1->g) << 8)/n;\
        binc=((p2->b-p1->b) << 8)/n);\
    a=2*dy-dx;\
    dy=2*dy;\
    dx=2*dx-dy;\
    pp_inc_1 = (inc_1) * PSZB;\
    pp_inc_2 = (inc_2) * PSZB;\
    do {\
        PUTPIXEL();\
        ZZ(z+=zinc);\
        RGB(r+=rinc;g+=ginc;b+=binc);\
        if (a>0) { pp=(PIXEL *)((char *)pp + pp_inc_1); ZZ(pz+=(inc_1));  a-=dx; }\
	else { pp=(PIXEL *)((char *)pp + pp_inc_2); ZZ(pz+=(inc_2)); a+=dy; }\
    } while (--n >= 0);

/* fin macro */

    if (dx == 0 && dy == 0) {
	PUTPIXEL();
    } else if (dx > 0) {
	if (dx >= dy) {
	    DRAWLINE(dx, dy, sx + 1, 1);
	} else {
	    DRAWLINE(dy, dx, sx + 1, sx);
	}
    } else {
	dx = -dx;
	if (dx >= dy) {
	    DRAWLINE(dx, dy, sx - 1, -1);
	} else {
	    DRAWLINE(dy, dx, sx - 1, sx);
	}
    }
}

#undef INTERP_Z
#undef INTERP_RGB

/* internal defines */
#undef DRAWLINE
#undef PUTPIXEL
#undef ZZ
#undef RGB
#undef RGBPIXEL 
