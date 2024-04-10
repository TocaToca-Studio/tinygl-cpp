#include <stdlib.h>

#include "zbuffer.hpp"

#define ZCMP(z, zpix) ((z) >= (zpix))

void ZB_fillTriangleFlat(ZBuffer *zb, ZBufferPoint *p0, ZBufferPoint *p1,
                         ZBufferPoint *p2) {
  int color;

  {
    ZBufferPoint *t, *pr1, *pr2, *l1, *l2;
    float fdx1, fdx2, fdy1, fdy2, fz, d1, d2;
    uint16_t *pz1;
    PIXEL *pp1;
    int part, update_left, update_right;

    int nb_lines, dx1, dy1, tmp, dx2, dy2;

    int error, derror;
    int x1, dxdy_min, dxdy_max;

    int x2, dx2dy2;

    int z1, dzdx, dzdy, dzdl_min, dzdl_max;
    if (p1->y < p0->y) {
      t = p0;
      p0 = p1;
      p1 = t;
    }
    if (p2->y < p0->y) {
      t = p2;
      p2 = p1;
      p1 = p0;
      p0 = t;
    } else if (p2->y < p1->y) {
      t = p1;
      p1 = p2;
      p2 = t;
    }

    fdx1 = p1->x - p0->x;
    fdy1 = p1->y - p0->y;

    fdx2 = p2->x - p0->x;
    fdy2 = p2->y - p0->y;

    fz = fdx1 * fdy2 - fdx2 * fdy1;
    if (fz == 0) return;
    fz = 1.0 / fz;

    fdx1 *= fz;
    fdy1 *= fz;
    fdx2 *= fz;
    fdy2 *= fz;

    d1 = p1->z - p0->z;
    d2 = p2->z - p0->z;
    dzdx = (int)(fdy2 * d1 - fdy1 * d2);
    dzdy = (int)(fdx1 * d2 - fdx2 * d1);
    pp1 = (PIXEL *)((char *)zb->pbuf + zb->linesize() * p0->y);
    pz1 = zb->zbuf + p0->y * zb->xsize;

    {
      color = RGB_TO_PIXEL(p2->r,p2->g,p2->b); 
    };

    for (part = 0; part < 2; part++) {
      if (part == 0) {
        if (fz > 0) {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p2;
          pr1 = p0;
          pr2 = p1;
        } else {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p1;
          pr1 = p0;
          pr2 = p2;
        }
        nb_lines = p1->y - p0->y;
      } else {
        if (fz > 0) {
          update_left = 0;
          update_right = 1;
          pr1 = p1;
          pr2 = p2;
        } else {
          update_left = 1;
          update_right = 0;
          l1 = p1;
          l2 = p2;
        }
        nb_lines = p2->y - p1->y + 1;
      }

      if (update_left) {
        dy1 = l2->y - l1->y;
        dx1 = l2->x - l1->x;
        if (dy1 > 0)
          tmp = (dx1 << 16) / dy1;
        else
          tmp = 0;
        x1 = l1->x;
        error = 0;
        derror = tmp & 0x0000ffff;
        dxdy_min = tmp >> 16;
        dxdy_max = dxdy_min + 1;

        z1 = l1->z;
        dzdl_min = (dzdy + dzdx * dxdy_min);
        dzdl_max = dzdl_min + dzdx;
      }

      if (update_right) {
        dx2 = (pr2->x - pr1->x);
        dy2 = (pr2->y - pr1->y);
        if (dy2 > 0)
          dx2dy2 = (dx2 << 16) / dy2;
        else
          dx2dy2 = 0;
        x2 = pr1->x << 16;
      }

      while (nb_lines > 0) {
        nb_lines--;

        {
          register PIXEL *pp;
          register int n;

          register uint16_t *pz;
          register uint32_t z, zz;

          n = (x2 >> 16) - x1;
          pp = (PIXEL *)((char *)pp1 + x1 * 2);

          pz = pz1 + x1;
          z = z1;

          while (n >= 3) {
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                pp[0] = color;
                pz[0] = zz;
              }
              z += dzdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[1]))) {
                pp[1] = color;
                pz[1] = zz;
              }
              z += dzdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[2]))) {
                pp[2] = color;
                pz[2] = zz;
              }
              z += dzdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[3]))) {
                pp[3] = color;
                pz[3] = zz;
              }
              z += dzdx;
            };

            pz += 4;

            pp = (PIXEL *)((char *)pp + 4 * 2);
            n -= 4;
          }
          while (n >= 0) {
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                pp[0] = color;
                pz[0] = zz;
              }
              z += dzdx;
            };

            pz += 1;

            pp = (PIXEL *)((char *)pp + 2);
            n -= 1;
          }
        }

        error += derror;
        if (error > 0) {
          error -= 0x10000;
          x1 += dxdy_max;

          z1 += dzdl_max;

        } else {
          x1 += dxdy_min;

          z1 += dzdl_min;
        }

        x2 += dx2dy2;

        pp1 = (PIXEL *)((char *)pp1 + zb->linesize());
        pz1 += zb->xsize;
      }
    }
  }
}

void ZB_fillTriangleSmooth(ZBuffer *zb, ZBufferPoint *p0, ZBufferPoint *p1,
                           ZBufferPoint *p2) {
  int _drgbdx;

  {
    ZBufferPoint *t, *pr1, *pr2, *l1, *l2;
    float fdx1, fdx2, fdy1, fdy2, fz, d1, d2;
    uint16_t *pz1;
    PIXEL *pp1;
    int part, update_left, update_right;

    int nb_lines, dx1, dy1, tmp, dx2, dy2;

    int error, derror;
    int x1, dxdy_min, dxdy_max;

    int x2, dx2dy2;

    int z1, dzdx, dzdy, dzdl_min, dzdl_max;

    int r1, drdx, drdy, drdl_min, drdl_max;
    int g1, dgdx, dgdy, dgdl_min, dgdl_max;
    int b1, dbdx, dbdy, dbdl_min, dbdl_max;

    if (p1->y < p0->y) {
      t = p0;
      p0 = p1;
      p1 = t;
    }
    if (p2->y < p0->y) {
      t = p2;
      p2 = p1;
      p1 = p0;
      p0 = t;
    } else if (p2->y < p1->y) {
      t = p1;
      p1 = p2;
      p2 = t;
    }

    fdx1 = p1->x - p0->x;
    fdy1 = p1->y - p0->y;

    fdx2 = p2->x - p0->x;
    fdy2 = p2->y - p0->y;

    fz = fdx1 * fdy2 - fdx2 * fdy1;
    if (fz == 0) return;
    fz = 1.0 / fz;

    fdx1 *= fz;
    fdy1 *= fz;
    fdx2 *= fz;
    fdy2 *= fz;

    d1 = p1->z - p0->z;
    d2 = p2->z - p0->z;
    dzdx = (int)(fdy2 * d1 - fdy1 * d2);
    dzdy = (int)(fdx1 * d2 - fdx2 * d1);

    d1 = p1->r - p0->r;
    d2 = p2->r - p0->r;
    drdx = (int)(fdy2 * d1 - fdy1 * d2);
    drdy = (int)(fdx1 * d2 - fdx2 * d1);

    d1 = p1->g - p0->g;
    d2 = p2->g - p0->g;
    dgdx = (int)(fdy2 * d1 - fdy1 * d2);
    dgdy = (int)(fdx1 * d2 - fdx2 * d1);

    d1 = p1->b - p0->b;
    d2 = p2->b - p0->b;
    dbdx = (int)(fdy2 * d1 - fdy1 * d2);
    dbdy = (int)(fdx1 * d2 - fdx2 * d1);

    pp1 = (PIXEL *)((char *)zb->pbuf + zb->linesize() * p0->y);
    pz1 = zb->zbuf + p0->y * zb->xsize;

    {
      _drgbdx = ((drdx / (1 << 6)) << 22) & 0xFFC00000;
      _drgbdx |= (dgdx / (1 << 5)) & 0x000007FF;
      _drgbdx |= ((dbdx / (1 << 7)) << 12) & 0x001FF000;
    };

    for (part = 0; part < 2; part++) {
      if (part == 0) {
        if (fz > 0) {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p2;
          pr1 = p0;
          pr2 = p1;
        } else {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p1;
          pr1 = p0;
          pr2 = p2;
        }
        nb_lines = p1->y - p0->y;
      } else {
        if (fz > 0) {
          update_left = 0;
          update_right = 1;
          pr1 = p1;
          pr2 = p2;
        } else {
          update_left = 1;
          update_right = 0;
          l1 = p1;
          l2 = p2;
        }
        nb_lines = p2->y - p1->y + 1;
      }

      if (update_left) {
        dy1 = l2->y - l1->y;
        dx1 = l2->x - l1->x;
        if (dy1 > 0)
          tmp = (dx1 << 16) / dy1;
        else
          tmp = 0;
        x1 = l1->x;
        error = 0;
        derror = tmp & 0x0000ffff;
        dxdy_min = tmp >> 16;
        dxdy_max = dxdy_min + 1;

        z1 = l1->z;
        dzdl_min = (dzdy + dzdx * dxdy_min);
        dzdl_max = dzdl_min + dzdx;

        r1 = l1->r;
        drdl_min = (drdy + drdx * dxdy_min);
        drdl_max = drdl_min + drdx;

        g1 = l1->g;
        dgdl_min = (dgdy + dgdx * dxdy_min);
        dgdl_max = dgdl_min + dgdx;

        b1 = l1->b;
        dbdl_min = (dbdy + dbdx * dxdy_min);
        dbdl_max = dbdl_min + dbdx;
      }

      if (update_right) {
        dx2 = (pr2->x - pr1->x);
        dy2 = (pr2->y - pr1->y);
        if (dy2 > 0)
          dx2dy2 = (dx2 << 16) / dy2;
        else
          dx2dy2 = 0;
        x2 = pr1->x << 16;
      }

      while (nb_lines > 0) {
        nb_lines--;

        {
          register uint16_t *pz;
          register PIXEL *pp;
          register uint32_t tmp, z, zz, rgb, drgbdx;
          register int n;
          n = (x2 >> 16) - x1;
          pp = pp1 + x1;
          pz = pz1 + x1;
          z = z1;
          rgb = (r1 << 16) & 0xFFC00000;
          rgb |= (g1 >> 5) & 0x000007FF;
          rgb |= (b1 << 5) & 0x001FF000;
          drgbdx = _drgbdx;
          while (n >= 3) {
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                tmp = rgb & 0xF81F07E0;
                pp[0] = tmp | (tmp >> 16);
                pz[0] = zz;
              }
              z += dzdx;
              rgb = (rgb + drgbdx) & (~0x00200800);
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[1]))) {
                tmp = rgb & 0xF81F07E0;
                pp[1] = tmp | (tmp >> 16);
                pz[1] = zz;
              }
              z += dzdx;
              rgb = (rgb + drgbdx) & (~0x00200800);
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[2]))) {
                tmp = rgb & 0xF81F07E0;
                pp[2] = tmp | (tmp >> 16);
                pz[2] = zz;
              }
              z += dzdx;
              rgb = (rgb + drgbdx) & (~0x00200800);
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[3]))) {
                tmp = rgb & 0xF81F07E0;
                pp[3] = tmp | (tmp >> 16);
                pz[3] = zz;
              }
              z += dzdx;
              rgb = (rgb + drgbdx) & (~0x00200800);
            };
            pz += 4;
            pp += 4;
            n -= 4;
          }
          while (n >= 0) {
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                tmp = rgb & 0xF81F07E0;
                pp[0] = tmp | (tmp >> 16);
                pz[0] = zz;
              }
              z += dzdx;
              rgb = (rgb + drgbdx) & (~0x00200800);
            };
            pz += 1;
            pp += 1;
            n -= 1;
          }
        };

        error += derror;
        if (error > 0) {
          error -= 0x10000;
          x1 += dxdy_max;

          z1 += dzdl_max;

          r1 += drdl_max;
          g1 += dgdl_max;
          b1 += dbdl_max;

        } else {
          x1 += dxdy_min;

          z1 += dzdl_min;

          r1 += drdl_min;
          g1 += dgdl_min;
          b1 += dbdl_min;
        }

        x2 += dx2dy2;

        pp1 = (PIXEL *)((char *)pp1 + zb->linesize());
        pz1 += zb->xsize;
      }
    }
  }
}


void ZB_fillTriangleMapping(ZBuffer *zb, ZBufferPoint *p0, ZBufferPoint *p1,
                            ZBufferPoint *p2) {
  PIXEL *texture;

  {
    ZBufferPoint *t, *pr1, *pr2, *l1, *l2;
    float fdx1, fdx2, fdy1, fdy2, fz, d1, d2;
    uint16_t *pz1;
    PIXEL *pp1;
    int part, update_left, update_right;

    int nb_lines, dx1, dy1, tmp, dx2, dy2;

    int error, derror;
    int x1, dxdy_min, dxdy_max;

    int x2, dx2dy2;

    int z1, dzdx, dzdy, dzdl_min, dzdl_max;

    int s1, dsdx, dsdy, dsdl_min, dsdl_max;
    int t1, dtdx, dtdy, dtdl_min, dtdl_max;

    if (p1->y < p0->y) {
      t = p0;
      p0 = p1;
      p1 = t;
    }
    if (p2->y < p0->y) {
      t = p2;
      p2 = p1;
      p1 = p0;
      p0 = t;
    } else if (p2->y < p1->y) {
      t = p1;
      p1 = p2;
      p2 = t;
    }

    fdx1 = p1->x - p0->x;
    fdy1 = p1->y - p0->y;

    fdx2 = p2->x - p0->x;
    fdy2 = p2->y - p0->y;

    fz = fdx1 * fdy2 - fdx2 * fdy1;
    if (fz == 0) return;
    fz = 1.0 / fz;

    fdx1 *= fz;
    fdy1 *= fz;
    fdx2 *= fz;
    fdy2 *= fz;

    d1 = p1->z - p0->z;
    d2 = p2->z - p0->z;
    dzdx = (int)(fdy2 * d1 - fdy1 * d2);
    dzdy = (int)(fdx1 * d2 - fdx2 * d1);

    d1 = p1->s - p0->s;
    d2 = p2->s - p0->s;
    dsdx = (int)(fdy2 * d1 - fdy1 * d2);
    dsdy = (int)(fdx1 * d2 - fdx2 * d1);

    d1 = p1->t - p0->t;
    d2 = p2->t - p0->t;
    dtdx = (int)(fdy2 * d1 - fdy1 * d2);
    dtdy = (int)(fdx1 * d2 - fdx2 * d1);

    pp1 = (PIXEL *)((char *)zb->pbuf + zb->linesize() * p0->y);
    pz1 = zb->zbuf + p0->y * zb->xsize;

    { texture = zb->current_texture; };

    for (part = 0; part < 2; part++) {
      if (part == 0) {
        if (fz > 0) {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p2;
          pr1 = p0;
          pr2 = p1;
        } else {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p1;
          pr1 = p0;
          pr2 = p2;
        }
        nb_lines = p1->y - p0->y;
      } else {
        if (fz > 0) {
          update_left = 0;
          update_right = 1;
          pr1 = p1;
          pr2 = p2;
        } else {
          update_left = 1;
          update_right = 0;
          l1 = p1;
          l2 = p2;
        }
        nb_lines = p2->y - p1->y + 1;
      }

      if (update_left) {
        dy1 = l2->y - l1->y;
        dx1 = l2->x - l1->x;
        if (dy1 > 0)
          tmp = (dx1 << 16) / dy1;
        else
          tmp = 0;
        x1 = l1->x;
        error = 0;
        derror = tmp & 0x0000ffff;
        dxdy_min = tmp >> 16;
        dxdy_max = dxdy_min + 1;

        z1 = l1->z;
        dzdl_min = (dzdy + dzdx * dxdy_min);
        dzdl_max = dzdl_min + dzdx;

        s1 = l1->s;
        dsdl_min = (dsdy + dsdx * dxdy_min);
        dsdl_max = dsdl_min + dsdx;

        t1 = l1->t;
        dtdl_min = (dtdy + dtdx * dxdy_min);
        dtdl_max = dtdl_min + dtdx;
      }

      if (update_right) {
        dx2 = (pr2->x - pr1->x);
        dy2 = (pr2->y - pr1->y);
        if (dy2 > 0)
          dx2dy2 = (dx2 << 16) / dy2;
        else
          dx2dy2 = 0;
        x2 = pr1->x << 16;
      }

      while (nb_lines > 0) {
        nb_lines--;

        {
          register PIXEL *pp;
          register int n;

          register uint16_t *pz;
          register uint32_t z, zz;

          register uint32_t s, t;

          n = (x2 >> 16) - x1;
          pp = (PIXEL *)((char *)pp1 + x1 * 2);

          pz = pz1 + x1;
          z = z1;

          s = s1;
          t = t1;

          while (n >= 3) {
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                pp[0] = texture[((t & 0x3FC00000) | s) >> ZB_POINT_Z_FRAC_BITS];
                pz[0] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[1]))) {
                pp[1] = texture[((t & 0x3FC00000) | s) >> ZB_POINT_Z_FRAC_BITS];
                pz[1] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[2]))) {
                pp[2] = texture[((t & 0x3FC00000) | s) >> ZB_POINT_Z_FRAC_BITS];
                pz[2] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[3]))) {
                pp[3] = texture[((t & 0x3FC00000) | s) >> ZB_POINT_Z_FRAC_BITS];
                pz[3] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };

            pz += 4;

            pp = (PIXEL *)((char *)pp + 4 * 2);
            n -= 4;
          }
          while (n >= 0) {
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                pp[0] = texture[((t & 0x3FC00000) | s) >> ZB_POINT_Z_FRAC_BITS];
                pz[0] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };

            pz += 1;

            pp = (PIXEL *)((char *)pp + 2);
            n -= 1;
          }
        }

        error += derror;
        if (error > 0) {
          error -= 0x10000;
          x1 += dxdy_max;

          z1 += dzdl_max;

          s1 += dsdl_max;
          t1 += dtdl_max;

        } else {
          x1 += dxdy_min;

          z1 += dzdl_min;

          s1 += dsdl_min;
          t1 += dtdl_min;
        }

        x2 += dx2dy2;

        pp1 = (PIXEL *)((char *)pp1 + zb->linesize());
        pz1 += zb->xsize;
      }
    }
  }
}

void ZB_fillTriangleMappingPerspective(ZBuffer *zb, ZBufferPoint *p0,
                                       ZBufferPoint *p1, ZBufferPoint *p2) {
  PIXEL *texture;
  float fdzdx, fndzdx, ndszdx, ndtzdx;

  {
    ZBufferPoint *t, *pr1, *pr2, *l1, *l2;
    float fdx1, fdx2, fdy1, fdy2, fz, d1, d2;
    uint16_t *pz1;
    PIXEL *pp1;
    int part, update_left, update_right;

    int nb_lines, dx1, dy1, tmp, dx2, dy2;

    int error, derror;
    int x1, dxdy_min, dxdy_max;

    int x2, dx2dy2;

    int z1, dzdx, dzdy, dzdl_min, dzdl_max;

    float sz1, dszdx, dszdy, dszdl_min, dszdl_max;
    float tz1, dtzdx, dtzdy, dtzdl_min, dtzdl_max;

    if (p1->y < p0->y) {
      t = p0;
      p0 = p1;
      p1 = t;
    }
    if (p2->y < p0->y) {
      t = p2;
      p2 = p1;
      p1 = p0;
      p0 = t;
    } else if (p2->y < p1->y) {
      t = p1;
      p1 = p2;
      p2 = t;
    }

    fdx1 = p1->x - p0->x;
    fdy1 = p1->y - p0->y;

    fdx2 = p2->x - p0->x;
    fdy2 = p2->y - p0->y;

    fz = fdx1 * fdy2 - fdx2 * fdy1;
    if (fz == 0) return;
    fz = 1.0 / fz;

    fdx1 *= fz;
    fdy1 *= fz;
    fdx2 *= fz;
    fdy2 *= fz;

    d1 = p1->z - p0->z;
    d2 = p2->z - p0->z;
    dzdx = (int)(fdy2 * d1 - fdy1 * d2);
    dzdy = (int)(fdx1 * d2 - fdx2 * d1);

    {
      float zz;
      zz = (float)p0->z;
      p0->sz = (float)p0->s * zz;
      p0->tz = (float)p0->t * zz;
      zz = (float)p1->z;
      p1->sz = (float)p1->s * zz;
      p1->tz = (float)p1->t * zz;
      zz = (float)p2->z;
      p2->sz = (float)p2->s * zz;
      p2->tz = (float)p2->t * zz;

      d1 = p1->sz - p0->sz;
      d2 = p2->sz - p0->sz;
      dszdx = (fdy2 * d1 - fdy1 * d2);
      dszdy = (fdx1 * d2 - fdx2 * d1);

      d1 = p1->tz - p0->tz;
      d2 = p2->tz - p0->tz;
      dtzdx = (fdy2 * d1 - fdy1 * d2);
      dtzdy = (fdx1 * d2 - fdx2 * d1);
    }

    pp1 = (PIXEL *)((char *)zb->pbuf + zb->linesize() * p0->y);
    pz1 = zb->zbuf + p0->y * zb->xsize;

    {
      texture = zb->current_texture;
      fdzdx = (float)dzdx;
      fndzdx = 8 * fdzdx;
      ndszdx = 8 * dszdx;
      ndtzdx = 8 * dtzdx;
    };

    for (part = 0; part < 2; part++) {
      if (part == 0) {
        if (fz > 0) {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p2;
          pr1 = p0;
          pr2 = p1;
        } else {
          update_left = 1;
          update_right = 1;
          l1 = p0;
          l2 = p1;
          pr1 = p0;
          pr2 = p2;
        }
        nb_lines = p1->y - p0->y;
      } else {
        if (fz > 0) {
          update_left = 0;
          update_right = 1;
          pr1 = p1;
          pr2 = p2;
        } else {
          update_left = 1;
          update_right = 0;
          l1 = p1;
          l2 = p2;
        }
        nb_lines = p2->y - p1->y + 1;
      }

      if (update_left) {
        dy1 = l2->y - l1->y;
        dx1 = l2->x - l1->x;
        if (dy1 > 0)
          tmp = (dx1 << 16) / dy1;
        else
          tmp = 0;
        x1 = l1->x;
        error = 0;
        derror = tmp & 0x0000ffff;
        dxdy_min = tmp >> 16;
        dxdy_max = dxdy_min + 1;

        z1 = l1->z;
        dzdl_min = (dzdy + dzdx * dxdy_min);
        dzdl_max = dzdl_min + dzdx;

        sz1 = l1->sz;
        dszdl_min = (dszdy + dszdx * dxdy_min);
        dszdl_max = dszdl_min + dszdx;

        tz1 = l1->tz;
        dtzdl_min = (dtzdy + dtzdx * dxdy_min);
        dtzdl_max = dtzdl_min + dtzdx;
      }

      if (update_right) {
        dx2 = (pr2->x - pr1->x);
        dy2 = (pr2->y - pr1->y);
        if (dy2 > 0)
          dx2dy2 = (dx2 << 16) / dy2;
        else
          dx2dy2 = 0;
        x2 = pr1->x << 16;
      }

      while (nb_lines > 0) {
        nb_lines--;

        {
          register uint16_t *pz;
          register PIXEL *pp;
          register uint32_t s, t, z, zz;
          register int n, dsdx, dtdx;
          float sz, tz, fz, zinv;
          n = (x2 >> 16) - x1;
          fz = (float)z1;
          zinv = 1.0 / fz;
          pp = (PIXEL *)((char *)pp1 + x1 * 2);
          pz = pz1 + x1;
          z = z1;
          sz = sz1;
          tz = tz1;
          while (n >= (8 - 1)) {
            {
              float ss, tt;
              ss = (sz * zinv);
              tt = (tz * zinv);
              s = (int)ss;
              t = (int)tt;
              dsdx = (int)((dszdx - ss * fdzdx) * zinv);
              dtdx = (int)((dtzdx - tt * fdzdx) * zinv);
              fz += fndzdx;
              zinv = 1.0 / fz;
            }
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                pp[0] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[0] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[1]))) {
                pp[1] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[1] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[2]))) {
                pp[2] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[2] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[3]))) {
                pp[3] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[3] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[4]))) {
                pp[4] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[4] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[5]))) {
                pp[5] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[5] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[6]))) {
                pp[6] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[6] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[7]))) {
                pp[7] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[7] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            pz += 8;
            pp = (PIXEL *)((char *)pp + 8 * 2);
            n -= 8;
            sz += ndszdx;
            tz += ndtzdx;
          }
          {
            float ss, tt;
            ss = (sz * zinv);
            tt = (tz * zinv);
            s = (int)ss;
            t = (int)tt;
            dsdx = (int)((dszdx - ss * fdzdx) * zinv);
            dtdx = (int)((dtzdx - tt * fdzdx) * zinv);
          }
          while (n >= 0) {
            {
              zz = z >> ZB_POINT_Z_FRAC_BITS;
              if (((zz) >= (pz[0]))) {
                pp[0] = *(PIXEL *)((char *)texture +
                                   (((t & 0x3FC00000) | (s & 0x003FC000)) >>
                                    (17 - 4)));
                pz[0] = zz;
              }
              z += dzdx;
              s += dsdx;
              t += dtdx;
            };
            pz += 1;
            pp = (PIXEL *)((char *)pp + 2);
            n -= 1;
          }
        };

        error += derror;
        if (error > 0) {
          error -= 0x10000;
          x1 += dxdy_max;

          z1 += dzdl_max;

          sz1 += dszdl_max;
          tz1 += dtzdl_max;

        } else {
          x1 += dxdy_min;

          z1 += dzdl_min;

          sz1 += dszdl_min;
          tz1 += dtzdl_min;
        }

        x2 += dx2dy2;

        pp1 = (PIXEL *)((char *)pp1 + zb->linesize());
        pz1 += zb->xsize;
      }
    }
  }
}
