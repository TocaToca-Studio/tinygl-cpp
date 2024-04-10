#ifndef _tgl_zbuffer_h_
#define _tgl_zbuffer_h_

/*
 * Z buffer
 */

#include "zfeatures.h" 

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
   
#define ZB_Z_BITS 16

#define ZB_POINT_Z_FRAC_BITS 14

#define ZB_POINT_S_MIN ( (1<<13) )
#define ZB_POINT_S_MAX ( (1<<22)-(1<<13) )
#define ZB_POINT_T_MIN ( (1<<21) )
#define ZB_POINT_T_MAX ( (1<<30)-(1<<21) )

#define ZB_POINT_RED_MIN ( (1<<10) )
#define ZB_POINT_RED_MAX ( (1<<16)-(1<<10) )
#define ZB_POINT_GREEN_MIN ( (1<<9) )
#define ZB_POINT_GREEN_MAX ( (1<<16)-(1<<9) )
#define ZB_POINT_BLUE_MIN ( (1<<10) )
#define ZB_POINT_BLUE_MAX ( (1<<16)-(1<<10) )
 
/* 16 bit mode */
#define RGB_TO_PIXEL(r,g,b) \
  (((r) & 0xF800) | (((g) >> 5) & 0x07E0) | ((b) >> 11))
typedef uint16_t PIXEL;
#define PSZB 2 
#define PSZSH 4 

#define ZCMP(z, zpix) ((z) >= (zpix))


/* 32 bpp copy */

#define RGB16_TO_RGB32(p0, p1, v)                     \
  {                                                   \
    uint32_t g, b, gb;                            \
    g = (v & 0x07E007E0) << 5;                        \
    b = (v & 0x001F001F) << 3;                        \
    gb = g | b;                                       \
    p0 = (gb & 0x0000FFFF) | ((v & 0x0000F800) << 8); \
    p1 = (gb >> 16) | ((v & 0xF8000000) >> 8);        \
  }

static inline uint32_t PIXEL_TO_RGB(uint16_t pixel) {

}  


typedef struct {
  int32_t x,y,z;     /* integer coordinates in the zbuffer */
  int32_t s,t;       /* coordinates for the mapping */
  int32_t r,g,b;     /* color indexes */
  
  float sz,tz;   /* temporary coordinates for mapping */
} ZBufferPoint;

 
struct  ZBuffer {
    uint16_t xsize,ysize;
    /** PRECISAMOS SE LIVRAR DESSA VARIAVEL CONFUSA*/ 
    int mode;
    /* zbuffer aparentemente, a distancia do buffer em um ponto fixo de 16 bits*/
    uint16_t *zbuf;
    /* PIXEL buffer, aparentemente uma imagem em rgb16*/
    PIXEL *pbuf;
    /* indica se o pbuff e o zbuff esta alocado*/
    bool frame_buffer_allocated;
     
    PIXEL *current_texture;

  static inline ZBuffer* open(uint16_t xsize, uint16_t ysize, void* frame_buffer) {
      ZBuffer* zb;
      int size;

      zb = (ZBuffer*)gl_malloc(sizeof(ZBuffer));
      if (zb == NULL) return NULL;

      zb->xsize = xsize;
      zb->ysize = ysize; 

      size = zb->xsize * zb->ysize * sizeof(uint16_t);

      zb->zbuf = (uint16_t*)gl_malloc(size);
      if (zb->zbuf == NULL)  { 
        gl_free(zb);
        return NULL;
      }
      if (frame_buffer == NULL) {
        zb->pbuf = (PIXEL*)gl_malloc(zb->ysize * zb->linesize());
        if (zb->pbuf == NULL) {
          gl_free(zb->zbuf);
          if (zb->zbuf == NULL)  { 
            gl_free(zb);
            return NULL;
          }
        }
        zb->frame_buffer_allocated = 1;
      } else {
        zb->frame_buffer_allocated = 0;
        zb->pbuf = (PIXEL*)frame_buffer;
      }
    
      zb->current_texture = NULL;

      return zb; 
    }

  static inline void close(ZBuffer* zb) {
    if (zb->frame_buffer_allocated) gl_free(zb->pbuf);

    gl_free(zb->zbuf);
    gl_free(zb);
  }
  inline void put_pixel(uint16_t x,uint16_t y,PIXEL color) {
    pbuf[x+(y*xsize)]=color;
  }
  inline PIXEL * pixel_pointer(uint16_t x,uint16_t y) {
    return &pbuf[x+(y*xsize)];
  }
  inline void put_z(uint16_t x,uint16_t y,uint16_t z) {
    zbuf[x+(y*xsize)]=z;
  }
  inline uint16_t get_z(uint16_t x,uint16_t y) {
    return zbuf[x+(y*xsize)];
  }
  inline uint32_t linesize() {
    return xsize*sizeof(PIXEL);
  }
  inline void clear(bool clear_z, uint16_t z, bool clear_color, uint8_t r, uint8_t g, uint8_t b) {
    PIXEL color = RGB_TO_PIXEL(r, g, b); 
    if (clear_z) memset(zbuf, z, xsize * ysize * sizeof(uint16_t));
    if (clear_color) { 
      for (uint16_t y = 0; y < ysize; y++) {
        PIXEL* pp=pbuf+ (y*xsize);
        for (uint16_t x = 0; x < xsize; x++) {
          pp[x]=color;
        }
      }
    }
  }
  
  inline void setTexture(PIXEL *texture) {
    current_texture = texture;
  }
  inline void resize(void* frame_buffer, uint16_t w, uint16_t h) {
    xsize = w;
    ysize = h; 

    size_t size = xsize * ysize * sizeof(uint16_t);

    gl_free(zbuf);
    zbuf = (uint16_t*)gl_malloc(size);

    if (frame_buffer_allocated) {
      gl_free(pbuf);
      frame_buffer_allocated=false;
      pbuf=NULL;
    }

    if (frame_buffer == NULL) {
      pbuf = (PIXEL*)gl_malloc(ysize * linesize());
      frame_buffer_allocated = 1;
    } else {
      pbuf = (PIXEL*)frame_buffer;
      frame_buffer_allocated = 0;
    }
  }
  inline void copyFrameBuffer(void* buf, int _linesize) {
  
    uint32_t *p, *p1, v, w0, w1;
    int y, n;

    uint16_t* q = pbuf;
    p1 = (uint32_t*)buf;

    for (y = 0; y < ysize; y++) {
      p = p1;
      n = xsize >> 2;
      do {
        v = *(uint32_t*)q;
        #if BYTE_ORDER == BIG_ENDIAN
              RGB16_TO_RGB32(w1, w0, v);
        #else
              RGB16_TO_RGB32(w0, w1, v);
        #endif
        p[0] = w0;
        p[1] = w1;

        v = *(uint32_t*)(q + 2);
        #if BYTE_ORDER == BIG_ENDIAN
              RGB16_TO_RGB32(w1, w0, v);
        #else
              RGB16_TO_RGB32(w0, w1, v);
        #endif
        p[2] = w0;
        p[3] = w1;

        q += 4;
        p += 4;
      } while (--n > 0);

      p1 = (uint32_t*)((char*)p1 + _linesize);
    }
  }
  inline bool zcompare(uint16_t x,uint16_t y,int32_t z) {
    return (z >> ZB_POINT_Z_FRAC_BITS) >= get_z(x,y);
  }
  inline void plot(ZBufferPoint *p) { 
    if (zcompare(p->x,p->y,p->z)) {
        put_pixel(p->x,p->y,RGB_TO_PIXEL(p->r, p->g, p->b)); 
    }
}
};
  
 
/* zline.c */
 
void ZB_line(ZBuffer *zb,ZBufferPoint *p1,ZBufferPoint *p2);
void ZB_line_z(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2);

/* ztriangle.c */
 

void ZB_fillTriangleFlat(ZBuffer *zb,
		 ZBufferPoint *p1,ZBufferPoint *p2,ZBufferPoint *p3);

void ZB_fillTriangleSmooth(ZBuffer *zb,
		   ZBufferPoint *p1,ZBufferPoint *p2,ZBufferPoint *p3);

void ZB_fillTriangleMapping(ZBuffer *zb,
		    ZBufferPoint *p1,ZBufferPoint *p2,ZBufferPoint *p3);

void ZB_fillTriangleMappingPerspective(ZBuffer *zb,
                    ZBufferPoint *p0,ZBufferPoint *p1,ZBufferPoint *p2);


typedef void (*ZB_fillTriangleFunc)(ZBuffer  *,
	    ZBufferPoint *,ZBufferPoint *,ZBufferPoint *);

 
#endif /* _tgl_zbuffer_h_ */
