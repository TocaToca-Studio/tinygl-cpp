#include <stdlib.h>

#include "zbuffer.hpp"

#define ZCMP(z, zpix) ((z) >= (zpix))

void ZB_plot(ZBuffer *zb, ZBufferPoint *p) {
  uint16_t *pz;
  PIXEL *pp;
  int zz;

  pz = zb->zbuf + (p->y * zb->xsize + p->x);
  pp = (PIXEL *)((char *)zb->pbuf + zb->linesize * p->y + p->x * PSZB);
  zz = p->z >> ZB_POINT_Z_FRAC_BITS;
  if (ZCMP(zz, *pz)) {
    *pp = RGB_TO_PIXEL(p->r, p->g, p->b);
    *pz = zz;
  }
}

static void ZB_line_flat_z(ZBuffer *zb, ZBufferPoint *p1, ZBufferPoint *p2,
                           int color) {
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
  pp = (PIXEL *)((char *)zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

  pz = zb->zbuf + (p1->y * sx + p1->x);
  z = p1->z;

  dx = p2->x - p1->x;
  dy = p2->y - p1->y;

  if (dx == 0 && dy == 0) {
    {
      zz = z >> ZB_POINT_Z_FRAC_BITS;
      if (((zz) >= (*pz))) {
        *pp = color;
        *pz = zz;
      }
    };
  } else if (dx > 0) {
    if (dx >= dy) {
      n = dx;
      zinc = (p2->z - p1->z) / n;
      ;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (1) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = color;
            *pz = zz;
          }
        };
        z += zinc;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx + 1);
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (1);
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      zinc = (p2->z - p1->z) / n;
      ;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = color;
            *pz = zz;
          }
        };
        z += zinc;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx + 1);
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (sx);
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  } else {
    dx = -dx;
    if (dx >= dy) {
      n = dx;
      zinc = (p2->z - p1->z) / n;
      ;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (-1) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = color;
            *pz = zz;
          }
        };
        z += zinc;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx - 1);
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (-1);
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      zinc = (p2->z - p1->z) / n;
      ;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = color;
            *pz = zz;
          }
        };
        z += zinc;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx - 1);
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (sx);
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  }
}

/* line with color interpolation */
static void ZB_line_interp_z(ZBuffer *zb, ZBufferPoint *p1, ZBufferPoint *p2,
                             bool interp_z = true, bool interp_rgb = true) {
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
  pp = (PIXEL *)((char *)zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

  pz = zb->zbuf + (p1->y * sx + p1->x);
  z = p1->z;

  dx = p2->x - p1->x;
  dy = p2->y - p1->y;

  r = p2->r << 8;
  g = p2->g << 8;
  b = p2->b << 8;

  if (dx == 0 && dy == 0) {
    {
      zz = z >> ZB_POINT_Z_FRAC_BITS;
      if (((zz) >= (*pz))) {
        *pp = RGB_TO_PIXEL(r, g, b);
        *pz = zz;
      }
    };
  } else if (dx > 0) {
    if (dx >= dy) {
      n = dx;
      zinc = (p2->z - p1->z) / n;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (1) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = RGB_TO_PIXEL(r, g, b);
            *pz = zz;
          }
        };
        z += zinc;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx + 1);
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (1);
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      zinc = (p2->z - p1->z) / n;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = RGB_TO_PIXEL(r, g, b);
            *pz = zz;
          }
        };
        z += zinc;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx + 1);
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (sx);
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  } else {
    dx = -dx;
    if (dx >= dy) {
      n = dx;
      zinc = (p2->z - p1->z) / n;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (-1) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = RGB_TO_PIXEL(r, g, b);
            *pz = zz;
          }
        };
        z += zinc;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx - 1);
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (-1);
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      zinc = (p2->z - p1->z) / n;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        {
          zz = z >> ZB_POINT_Z_FRAC_BITS;
          if (((zz) >= (*pz))) {
            *pp = RGB_TO_PIXEL(r, g, b);
            *pz = zz;
          }
        };
        z += zinc;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          pz += (sx - 1);
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          pz += (sx);
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  }
}

/* no Z interpolation */

static void ZB_line_flat(ZBuffer *zb, ZBufferPoint *p1, ZBufferPoint *p2,
                         int color) {
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
  pp = (PIXEL *)((char *)zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

  dx = p2->x - p1->x;
  dy = p2->y - p1->y;

  if (dx == 0 && dy == 0) {
    *pp = color;
  } else if (dx > 0) {
    if (dx >= dy) {
      n = dx;
      ;
      ;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (1) * 2;
      do {
        *pp = color;
        ;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      ;
      ;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        *pp = color;
        ;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  } else {
    dx = -dx;
    if (dx >= dy) {
      n = dx;
      ;
      ;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (-1) * 2;
      do {
        *pp = color;
        ;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      ;
      ;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        *pp = color;
        ;
        ;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  }
}

static void ZB_line_interp(ZBuffer *zb, ZBufferPoint *p1, ZBufferPoint *p2) {
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
  pp = (PIXEL *)((char *)zb->pbuf + zb->linesize * p1->y + p1->x * PSZB);

  dx = p2->x - p1->x;
  dy = p2->y - p1->y;

  r = p2->r << 8;
  g = p2->g << 8;
  b = p2->b << 8;

  if (dx == 0 && dy == 0) {
    *pp = RGB_TO_PIXEL(r, g, b);
  } else if (dx > 0) {
    if (dx >= dy) {
      n = dx;
      ;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (1) * 2;
      do {
        *pp = RGB_TO_PIXEL(r, g, b);
        ;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      ;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx + 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        *pp = RGB_TO_PIXEL(r, g, b);
        ;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  } else {
    dx = -dx;
    if (dx >= dy) {
      n = dx;
      ;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dy - dx;
      dy = 2 * dy;
      dx = 2 * dx - dy;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (-1) * 2;
      do {
        *pp = RGB_TO_PIXEL(r, g, b);
        ;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dx;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dy;
        }
      } while (--n >= 0);
      ;
    } else {
      n = dy;
      ;
      rinc = ((p2->r - p1->r) << 8) / n;
      ginc = ((p2->g - p1->g) << 8) / n;
      binc = ((p2->b - p1->b) << 8) / n;
      a = 2 * dx - dy;
      dx = 2 * dx;
      dy = 2 * dy - dx;
      pp_inc_1 = (sx - 1) * 2;
      pp_inc_2 = (sx) * 2;
      do {
        *pp = RGB_TO_PIXEL(r, g, b);
        ;
        r += rinc;
        g += ginc;
        b += binc;
        if (a > 0) {
          pp = (PIXEL *)((char *)pp + pp_inc_1);
          ;
          a -= dy;
        } else {
          pp = (PIXEL *)((char *)pp + pp_inc_2);
          ;
          a += dx;
        }
      } while (--n >= 0);
      ;
    }
  }
}

void ZB_line_z(ZBuffer *zb, ZBufferPoint *p1, ZBufferPoint *p2) {
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

void ZB_line(ZBuffer *zb, ZBufferPoint *p1, ZBufferPoint *p2) {
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
