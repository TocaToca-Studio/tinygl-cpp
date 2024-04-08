#include <stdlib.h>
#include "zbuffer.hpp"

#define ZCMP(z,zpix) ((z) >= (zpix))

void ZB_plot(ZBuffer * zb, ZBufferPoint * p)
{
    uint16_t *pz;
    PIXEL *pp;
    int zz;

    pz = zb->zbuf + (p->y * zb->xsize + p->x);
    pp = (PIXEL *) ((char *) zb->pbuf + zb->linesize * p->y + p->x * PSZB);
    zz = p->z >> ZB_POINT_Z_FRAC_BITS;
    if (ZCMP(zz, *pz)) {
        *pp = RGB_TO_PIXEL(p->r, p->g, p->b);
        *pz = zz;
    }
}

 
#define INTERP_Z
static void ZB_line_flat_z(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2, 
                           int color)
{
    int n, dx, dy, sx, pp_inc_1, pp_inc_2;
    register int a;
    register PIXEL *pp;

    register uint16_t *pz;
    int zinc;
    register int z, zz;

    if (p1->y > p2->y || (p1->y == p2->y && p1->x > p2->x)) {
        ZBufferPoint *tmp;
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    sx = zb->xsize;
    pp = (PIXEL *) ((char *) zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

    pz = zb->zbuf + (p1->y * sx + p1->x);
    z = p1->z;

    dx = p2->x - p1->x;
    dy = p2->y - p1->y;

#include "zline.h"
}



/* line with color interpolation */
#define INTERP_Z
#define INTERP_RGB
static void ZB_line_interp_z(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2,bool interp_z=true,bool interp_rgb=true)
{
    int n, dx, dy, sx, pp_inc_1, pp_inc_2;
    register int a;
    register PIXEL *pp;
    
    register uint32_t r, g, b;
    
    register uint32_t rinc, ginc, binc;
    
    register uint16_t *pz;
    int zinc;

    register int z, zz;

    if (p1->y > p2->y || (p1->y == p2->y && p1->x > p2->x)) {
        ZBufferPoint *tmp;
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    sx = zb->xsize;
    pp = (PIXEL *) ((char *) zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

    pz = zb->zbuf + (p1->y * sx + p1->x);
    z = p1->z;

    dx = p2->x - p1->x;
    dy = p2->y - p1->y;


    r = p2->r << 8;
    g = p2->g << 8;
    b = p2->b << 8;

#include "zline.h"
}

/* no Z interpolation */

static void ZB_line_flat(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2, 
                             int color)
{
        int n, dx, dy, sx, pp_inc_1, pp_inc_2;
    register int a;
    register PIXEL *pp;

    if (p1->y > p2->y || (p1->y == p2->y && p1->x > p2->x)) {
        ZBufferPoint *tmp;
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    sx = zb->xsize;
    pp = (PIXEL *) ((char *) zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

    dx = p2->x - p1->x;
    dy = p2->y - p1->y;

#include "zline.h"
}

#define INTERP_RGB
static void ZB_line_interp(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2)
{
        int n, dx, dy, sx, pp_inc_1, pp_inc_2;
    register int a;
    register PIXEL *pp;
    
    register uint32_t r, g, b;
    

    register uint32_t rinc, ginc, binc;
    
    if (p1->y > p2->y || (p1->y == p2->y && p1->x > p2->x)) {
        ZBufferPoint *tmp;
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    sx = zb->xsize;
    pp = (PIXEL *) ((char *) zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

    dx = p2->x - p1->x;
    dy = p2->y - p1->y;


    r = p2->r << 8;
    g = p2->g << 8;
    b = p2->b << 8;
#include "zline.h"
}

void ZB_line_z(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2)
{
    int color1, color2;

    color1 = RGB_TO_PIXEL(p1->r, p1->g, p1->b);
    color2 = RGB_TO_PIXEL(p2->r, p2->g, p2->b);

    /* choose if the line should have its color interpolated or not */
    if (color1 == color2) {
        ZB_line_flat_z(zb, p1, p2, color1);
    } else {
        ZB_line_interp_z(zb, p1, p2);
    }
}

void ZB_line(ZBuffer * zb, ZBufferPoint * p1, ZBufferPoint * p2)
{
    int color1, color2;

    color1 = RGB_TO_PIXEL(p1->r, p1->g, p1->b);
    color2 = RGB_TO_PIXEL(p2->r, p2->g, p2->b);

    /* choose if the line should have its color interpolated or not */
    if (color1 == color2) {
        ZB_line_flat(zb, p1, p2, color1);
    } else {
        ZB_line_interp(zb, p1, p2);
    }
}
