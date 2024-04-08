#ifndef __ZMATH__
#define __ZMATH__

#include <math.h>

struct vec3_t {
  float X, Y, Z;
  inline vec3_t() {
    X = 0;
    Y = 0;
    Z = 0;
  }
  inline vec3_t(float x, float y, float z) {
    X = x;
    Y = y;
    Z = z;
  }
  inline float Len() { return sqrt(X * X + Y * Y + Z * Z); }
  inline void Norm() {
    float n;
    n = Len();
    if (n == 0) return;
    X /= n;
    Y /= n;
    Z /= n;
  }
};

/* vector arithmetic */

struct vec4_t {
  float X, Y, Z, W;
  inline vec4_t() {
    X = 0;
    Y = 0;
    Z = 0;
    W = 0;
  }
  inline vec4_t(float x, float y, float z, float w) {
    X = x;
    Y = y;
    Z = z;
    W = w;
  }
};

/* Inversion of a general nxn matrix.
   Note : m is destroyed */

static inline int Matrix_Inv(float *r, float *m, int n) {
  int i, j, k, l;
  float max, tmp, t;

  /* identit�e dans r */
  for (i = 0; i < n * n; i++) r[i] = 0;
  for (i = 0; i < n; i++) r[i * n + i] = 1;

  for (j = 0; j < n; j++) {
    /* recherche du nombre de plus grand module sur la colonne j */
    max = m[j * n + j];
    k = j;
    for (i = j + 1; i < n; i++)
      if (fabs(m[i * n + j]) > fabs(max)) {
        k = i;
        max = m[i * n + j];
      }

    /* non intersible matrix */
    if (max == 0) return 1;

    /* permutation des lignes j et k */
    if (k != j) {
      for (i = 0; i < n; i++) {
        tmp = m[j * n + i];
        m[j * n + i] = m[k * n + i];
        m[k * n + i] = tmp;

        tmp = r[j * n + i];
        r[j * n + i] = r[k * n + i];
        r[k * n + i] = tmp;
      }
    }

    /* multiplication de la ligne j par 1/max */
    max = 1 / max;
    for (i = 0; i < n; i++) {
      m[j * n + i] *= max;
      r[j * n + i] *= max;
    }

    for (l = 0; l < n; l++)
      if (l != j) {
        t = m[l * n + j];
        for (i = 0; i < n; i++) {
          m[l * n + i] -= m[j * n + i] * t;
          r[l * n + i] -= r[j * n + i] * t;
        }
      }
  }

  return 0;
}

/* Matrix & Vertex */

struct mat4_t {
  float m[4][4];

  inline void Id() {
    int i, j;
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
        if (i == j)
          m[i][j] = 1.0;
        else
          m[i][j] = 0.0;
  }
  inline bool isId() {
    int i, j;
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++) {
        if (i == j) {
          if (m[i][j] != 1.0) return false;
        } else if (m[i][j] != 0.0)
          return false;
      }
    return true;
  }
  inline void Rotate(float t, int u) {
    float s, c;
    int v, w;
    if ((v = u + 1) > 2) v = 0;
    if ((w = v + 1) > 2) w = 0;
    s = sin(t);
    c = cos(t);
    Id();
    m[v][v] = c;
    m[v][w] = -s;
    m[w][v] = s;
    m[w][w] = c;
  }
  inline vec4_t Mulvec4_t(const vec4_t &c) {
    vec4_t a;
    a.X = m[0][0] * c.X + m[0][1] * c.Y + m[0][2] * c.Z + m[0][3] * c.W;
    a.Y = m[1][0] * c.X + m[1][1] * c.Y + m[1][2] * c.Z + m[1][3] * c.W;
    a.Z = m[2][0] * c.X + m[2][1] * c.Y + m[2][2] * c.Z + m[2][3] * c.W;
    a.W = m[3][0] * c.X + m[3][1] * c.Y + m[3][2] * c.Z + m[3][3] * c.W;
    return a;
  }
  /* inversion of a 4x4 matrix */

  inline void Inv(mat4_t *b) {
    mat4_t tmp = *b;
    /*tmp=*b;*/
    Matrix_Inv(&m[0][0], &tmp.m[0][0], 4);
  }
  inline void Inv() { Inv(this); }

  static inline mat4_t Mul(mat4_t *a, mat4_t *b) {
    mat4_t c;
    int i, j, k;
    float s;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        s = 0.0;
        for (k = 0; k < 4; k++) s += a->m[i][k] * b->m[k][j];
        c.m[i][j] = s;
      }
    }
    return c;
  }

  /* c=c*a */
  static inline void MulLeft(mat4_t *c, mat4_t *b) {
    int i, j, k;
    float s;
    mat4_t a = *c;

    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        s = 0.0;
        for (k = 0; k < 4; k++) s += a.m[i][k] * b->m[k][j];
        c->m[i][j] = s;
      }
    }
  }
  inline void Mult(mat4_t *b) { MulLeft(this, b); }

  //* transposition of a 4x4 matrix */
  inline void Transpose(mat4_t *b) {
    m[0][0] = b->m[0][0];
    m[0][1] = b->m[1][0];
    m[0][2] = b->m[2][0];
    m[0][3] = b->m[3][0];

    m[1][0] = b->m[0][1];
    m[1][1] = b->m[1][1];
    m[1][2] = b->m[2][1];
    m[1][3] = b->m[3][1];

    m[2][0] = b->m[0][2];
    m[2][1] = b->m[1][2];
    m[2][2] = b->m[2][2];
    m[2][3] = b->m[3][2];

    m[3][0] = b->m[0][3];
    m[3][1] = b->m[1][3];
    m[3][2] = b->m[2][3];
    m[3][3] = b->m[3][3];
  }
  inline void Transpose() { Transpose(this); }

  static inline mat4_t Frustrum(float left, float right, float bottom,
                                float top, float near, float farp) {
    float x, y, A, B, C, D;
    float *r;
    mat4_t m;
    x = (2.0 * near) / (right - left);
    y = (2.0 * near) / (top - bottom);
    A = (right + left) / (right - left);
    B = (top + bottom) / (top - bottom);
    C = -(farp + near) / (farp - near);
    D = -(2.0 * farp * near) / (farp - near);

    r = &m.m[0][0];
    r[0] = x;    r[1] = 0;    r[2] = A;    r[3] = 0;
    r[4] = 0;    r[5] = y;    r[6] = B;    r[7] = 0;
    r[8] = 0;    r[9] = 0;    r[10] = C;    r[11] = D;
    r[12] = 0;    r[13] = 0;    r[14] = -1;    r[15] = 0;
    return m;
  }
};

struct M3 {
  float m[3][3];
  inline float Det() {
    return m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1] -
           m[1][0] * m[0][1] * m[2][2] + m[1][0] * m[0][2] * m[2][1] +
           m[2][0] * m[0][1] * m[1][2] - m[2][0] * m[0][2] * m[1][1];
  };

  /* inverse of a 3x3 matrix */
  inline M3 Inv() {
    float det;
    M3 a;
    det = Det();

    a.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
    a.m[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / det;
    a.m[0][2] = -(-m[0][1] * m[1][2] + m[0][2] * m[1][1]) / det;

    a.m[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / det;
    a.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
    a.m[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) / det;

    a.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
    a.m[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) / det;
    a.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
    return a;
  }
};

struct rgba_t {
  float R, G, B, A;
  rgba_t() {
    R = 0;
    G = 0;
    B = 0;
    A = 0;
  }
  rgba_t(float r, float g, float b, float a) {
    R = r;
    G = g;
    B = b;
    A = a;
  }
};

#endif /* __ZMATH__ */
