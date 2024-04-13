#ifndef __ZMATH__
#define __ZMATH__

#include <math.h>


const float TH_PI = 3.14159265358979323846;
const float TH_HALF_PI = 1.57079632679489661923;
const float TH_TWO_PI = 6.28318530717958647692;
const float TH_DEG_TO_RAD = 0.017453292519943295769236907684886;
const float TH_RAD_TO_DEG = 57.295779513082320876798154814105;
const float TH_EULER = 2.718281828459045235360287471352;

static inline float DEGTORAD(float deg) { return deg * TH_DEG_TO_RAD;}
static inline float RADTODEG(float rad) { return rad * TH_RAD_TO_DEG;}
static inline float SQUARE(float x) { return x * x; }
static inline float CLAMP(float x, float min, float max) { return x < min ? min : (x > max ? max : x); }
static inline int32_t CLAMP(int32_t x, int32_t min, int32_t max) { return x < min ? min : (x > max ? max : x); }
static inline float FAST_LERP(float a,float b,float t) { return a + (b - a) * t; }
static inline float LEFP(float a,float b,float t) { return ((1.0f-t)*a+(t*b)); }
static inline float ABS(float x) { return x < 0 ? -x : x; } 
static inline float FLOAT_EQUALS(float a, float b) { return abs(a - b) < 1e-6; }
static inline float MIN(float a, float b) { return a < b ? a : b; }
static inline float MAX(float a, float b) { return a > b ? a : b; }
static inline float SIGN(float x) { return x < 0 ? -1 : 1; }
static inline uint8_t FLOAT_TO_BYTE(float f) {return (uint8_t)(CLAMP(f,0.0f,1.0f)*255.0);}
static inline float BYTE_TO_FLOAT(uint8_t b) { return b/255.0f; }
static inline uint8_t BYTE_INTERPOLATE(uint8_t a,uint8_t b,float factor) {
    return (uint8_t) FAST_LERP((float) a,(float) b,factor);
}
/** EXPERIMENTAL THIS MAY NOT WORK PROPERLY*/
static inline float COS_FROM_SIN(float _sin, float _angle) {
    float _cos=sqrt(1-SQUARE(_sin));
    float a=_angle+TH_HALF_PI;
    float b = a - ((int)(a / TH_TWO_PI) * TH_TWO_PI);
    if (b < 0) b += TH_TWO_PI;
    return (b >= TH_PI) ? -_cos : _cos;
}


static inline uint32_t LCM(uint32_t a,uint32_t b) {
    for(int i=2;i<=MIN(a,b);i++)  {
        if((a % i == 0) && (b%i == 0))
        return i;
    }
    return 1;
}


struct vec3f_t {
    float x, y, z;
    inline vec3f_t() { x = y = z = 0; }
    inline vec3f_t(float x, float y, float z) { this->x = x; this->y = y; this->z = z; }

    /* overload operators */
    vec3f_t operator+(const vec3f_t& v) const { return vec3f_t(x + v.x, y + v.y, z + v.z); }
    vec3f_t operator-(const vec3f_t& v) const { return vec3f_t(x - v.x, y - v.y, z - v.z); }
    vec3f_t operator*(const vec3f_t& v) const { return vec3f_t(x * v.x, y * v.y, z * v.z); }
    vec3f_t operator/(const vec3f_t& v) const { return vec3f_t(x / v.x, y / v.y, z / v.z); }
    vec3f_t operator+=(const vec3f_t v) { this->x+=v.x; this->y+=v.y; this->z+=v.z; return *this; }
    vec3f_t operator-=(const vec3f_t v) { this->x-= v.x; this->y-= v.y; this->z-= v.z; return *this; }
    vec3f_t operator*=(const vec3f_t v) { this->x*= v.x; this->y*= v.y; this->z*= v.z; return *this; }
    vec3f_t operator/=(const vec3f_t v) { this->x/= v.x; this->y/= v.y; this->z/= v.z; return *this; }
    vec3f_t operator-() const { return vec3f_t(-x, -y, -z); }
    vec3f_t operator+(float s) const { return vec3f_t(x + s, y + s, z + s); }
    vec3f_t operator-(float s) const { return vec3f_t(x - s, y - s, z - s); }
    vec3f_t operator*(float s) const { return vec3f_t(x * s, y * s, z * s); }
    vec3f_t operator/(float s) const { return vec3f_t(x / s, y / s, z / s); }
    vec3f_t operator+=(float s)  { this->x+=s; this->y+=s; this->z+=s; return *this; }
    vec3f_t operator-=(float s)  { this->x-=s; this->y-=s; this->z-=s; return *this;}
    vec3f_t operator*=(float s)  { this->x*=s; this->y*=s; this->z*=s; return *this; }
    vec3f_t operator/=(float s)  { this->x/=s; this->y/=s; this->z/=s; return *this; }
    vec3f_t operator+(double s) const { return vec3f_t(x + s, y + s, z + s); }
    vec3f_t operator-(double s) const { return vec3f_t(x - s, y - s, z - s); }
    vec3f_t operator*(double s) const { return vec3f_t(x * s, y * s, z * s); }
    vec3f_t operator/(double s) const { return vec3f_t(x / s, y / s, z / s); }
    vec3f_t operator+=(double s)  { this->x+=s; this->y+=s; this->z+=s; return *this; }
    vec3f_t operator-=(double s)  { this->x-=s; this->y-=s; this->z-=s; return *this;}
    vec3f_t operator*=(double s)  { this->x*=s; this->y*=s; this->z*=s; return *this; }
    vec3f_t operator/=(double s)  { this->x/=s; this->y/=s; this->z/=s; return *this; }
    vec3f_t operator+(int s) const { return vec3f_t(x + s, y + s, z + s); }
    vec3f_t operator-(int s) const { return vec3f_t(x - s, y - s, z - s); }
    vec3f_t operator*(int s) const { return vec3f_t(x * s, y * s, z * s); }
    vec3f_t operator/(int s) const { return vec3f_t(x / s, y / s, z / s); }
    vec3f_t operator+=(int s)  { this->x+=s; this->y+=s; this->z+=s; return *this; }
    vec3f_t operator-=(int s)  { this->x-=s; this->y-=s; this->z-=s; return *this;}
    vec3f_t operator*=(int s)  { this->x*=s; this->y*=s; this->z*=s; return *this; }
    vec3f_t operator/=(int s)  { this->x/=s; this->y/=s; this->z/=s; return *this; }

    bool operator==(const vec3f_t& v) const { return FLOAT_EQUALS(x, v.x) && FLOAT_EQUALS(y, v.y) && FLOAT_EQUALS(z, v.z); }
    bool operator!=(const vec3f_t& v) const { return !FLOAT_EQUALS(x, v.x) || !FLOAT_EQUALS(y, v.y) || !FLOAT_EQUALS(z, v.z); }

    inline vec3f_t opposite() const { return vec3f_t(-x, -y, -z); }
    inline vec3f_t inverse() const { return vec3f_t(-x, -y, -z); }
    inline float dot() const { return x * x + y * y + z * z; }
    inline float dot(vec3f_t v) const { return x * v.x + y * v.y + z * v.z; }
    inline float length() const { return sqrt(dot()); }
    inline float magnitude() const { return length(); }
    inline vec3f_t normalized() const { return (*this)* (1.0f/magnitude()); }
    inline void normalize() {
      float n=1.0f/length();
      x*=n; y*=n;z*=n;
    }
    inline vec3f_t abs() const { return vec3f_t(ABS(x), ABS(y), ABS(z)); }
    //inline color4f_t to_color() { return color4f_t(x,y,z,1.0f);}
    inline vec3f_t project(vec3f_t b) {
        return project(*this,b);
    }
    static inline vec3f_t cross(vec3f_t a,vec3f_t b ) {
        return vec3f_t(
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x
        );
    }
    static inline float dot(vec3f_t a,vec3f_t b) {
        return a.x * b.x + a.y * b.y + a.x * b.z;
    }
    // Funcao para calcular a projeçao de um vetor b sobre um vetor a
    static inline vec3f_t project(vec3f_t a,vec3f_t b) { 
        return a * (dot(a, b) / a.dot()); 
    }
    static inline vec3f_t ZERO() { return vec3f_t(0, 0, 0); }
    static inline vec3f_t UNIT() { return vec3f_t(1, 1, 1); }
    static inline vec3f_t X() { return vec3f_t(1, 0, 0); }
    static inline vec3f_t Y()       { return vec3f_t(0, 1, 0); }
    static inline vec3f_t Z()       { return vec3f_t(0, 0, 1); }
    static inline vec3f_t UP()      { return vec3f_t(0, 1, 0); }
    static inline vec3f_t DOWN()    { return vec3f_t(0, -1, 0); }
    static inline vec3f_t LEFT()    { return vec3f_t(-1, 0, 0); }
    static inline vec3f_t RIGHT()   { return vec3f_t(1, 0, 0); }
    static inline vec3f_t FORWARD() { return vec3f_t(0, 0, 1); }
    static inline vec3f_t BACK()    { return vec3f_t(0, 0, -1); }
    static inline float distance(vec3f_t a,vec3f_t b) { return (a-b).magnitude(); }
    // Function to calculate the barycentric coordinates of a point with respect to a triangle in r^3
    static inline vec3f_t barycentric(vec3f_t p, vec3f_t a, vec3f_t b, vec3f_t c) {
        vec3f_t v0 = b - a;
        vec3f_t v1 = c - a;
        vec3f_t v2 = p - a;

        float d00 = v0.dot(v0);
        float d01 = v0.dot(v1);
        float d11 = v1.dot(v1);
        float d20 = v2.dot(v0);
        float d21 = v2.dot(v1);

        float denom = d00 * d11 - d01 * d01;
        
        float v = (d11 * d20 - d01 * d21) / denom;
        float w = (d00 * d21 - d01 * d20) / denom;
        float u = 1.0f - v - w;

        return vec3f_t(u, v, w);
    }
};
/* vector arithmetic */

struct vec4f_t {
  float x, y, z, w;
  inline vec4f_t() {
    x = 0;
    y = 0;
    z = 0;
    w = 0;
  }
  inline vec4f_t(float vx, float vy, float vz, float vw) {
    x = vx;
    y = vy;
    z = vz;
    w = vw;
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
  inline vec4f_t Mulvec4_t(const vec4f_t &c) {
    vec4f_t a;
    a.x = m[0][0] * c.x + m[0][1] * c.y + m[0][2] * c.z + m[0][3] * c.w;
    a.y = m[1][0] * c.x + m[1][1] * c.y + m[1][2] * c.z + m[1][3] * c.w;
    a.z = m[2][0] * c.x + m[2][1] * c.y + m[2][2] * c.z + m[2][3] * c.w;
    a.w = m[3][0] * c.x + m[3][1] * c.y + m[3][2] * c.z + m[3][3] * c.w;
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
    float x, y, a, b, C, D;
    float *r;
    mat4_t m;
    x = (2.0 * near) / (right - left);
    y = (2.0 * near) / (top - bottom);
    a = (right + left) / (right - left);
    b = (top + bottom) / (top - bottom);
    C = -(farp + near) / (farp - near);
    D = -(2.0 * farp * near) / (farp - near);

    r = &m.m[0][0];
    r[0] = x;    r[1] = 0;    r[2] = a;    r[3] = 0;
    r[4] = 0;    r[5] = y;    r[6] = b;    r[7] = 0;
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

struct color4f_t {
    float r,g,b,a;
    inline color4f_t() { r = g = b = a = 1; }
    inline color4f_t(float red, float green, float blue, float alpha=1.0) { r = red; g = green; b = blue; a = alpha; }

    /*inline rgba_t to_bytes() {
        return rgba_t(FLOAT_TO_BYTE(r),FLOAT_TO_BYTE(g),FLOAT_TO_BYTE(b),FLOAT_TO_BYTE(a));
    }
    static inline color4f_t from_bytes(rgba_t bytes) {
        return color4f_t(BYTE_TO_FLOAT(bytes.r),BYTE_TO_FLOAT(bytes.g),BYTE_TO_FLOAT(bytes.b),BYTE_TO_FLOAT(bytes.a));
    }*/
    static inline color4f_t WHITE()     { return color4f_t(1,1,1,1); }
    static inline color4f_t BLACK()     { return color4f_t(0,0,0,1); }
    static inline color4f_t RED()       { return color4f_t(1,0,0,1); }
    static inline color4f_t GREEN()     { return color4f_t(0,1,0,1); }
    static inline color4f_t BLUE()      { return color4f_t(0,0,1,1); }
    static inline color4f_t YELLOW()    { return color4f_t(1,1,0,1); }
    static inline color4f_t CYAN()      { return color4f_t(0,1,1,1); }
    static inline color4f_t MAGENTA()   { return color4f_t(1,0,1,1); }
    static inline color4f_t GRAY()      { return color4f_t(0.5,0.5,0.5,1); }
    static inline color4f_t LIGHT_GRAY(){ return color4f_t(0.75,0.75,0.75,1); }
    static inline color4f_t DARK_GRAY() { return color4f_t(0.25,0.25,0.25,1); }
    static inline color4f_t ORANGE()    { return color4f_t(1,0.5,0,1); }
    static inline color4f_t PURPLE()    { return color4f_t(0.5,0,0.5,1); }
    static inline color4f_t PINK()      { return color4f_t(1,0.75,0.75,1); }
    static inline color4f_t BROWN()     { return color4f_t(0.5,0.25,0,1); }
    static inline color4f_t LIME()      { return color4f_t(0,1,0,1); }
    /*static inline color4f_t interpolate(rgba_t a, rgba_t b,float factor) {
        return color4f_t(
            FAST_LERP(a.r,b.r,factor),
            FAST_LERP(a.g,b.g,factor),
            FAST_LERP(a.b,b.b,factor),
            FAST_LERP(a.a,b.a,factor)
        );
    }*/
};


#endif /* __ZMATH__ */
