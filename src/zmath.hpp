#ifndef __ZMATH__
#define __ZMATH__

/* Some simple mathematical functions. Don't look for some logic in
   the function names :-) */

#include <stdlib.h>
#include <string.h>
#include <math.h> 


/* Matrix & Vertex */

struct M4 {
	float m[4][4]; 
	inline void Id() { 
		int i,j;
		for(i=0;i<4;i++)
		for(j=0;j<4;j++) 
		if (i==j) m[i][j]=1.0; else m[i][j]=0.0;
	}
	inline bool isId() {
		int i,j;
		for(i=0;i<4;i++)
		for(j=0;j<4;j++) {
			if (i==j) { 
				if (m[i][j] != 1.0) return false;
			} else if (m[i][j] != 0.0) return false;
		}
		return true;
	}
};

struct M3 {
	float m[3][3];
};

struct M34 {
	 float m[3][4];
};

 

struct V3 {
	float X,Y,Z;
	V3() {
		X=0;Y=0;Z=0;
	}
	inline float v(unsigned int index) {
		if(index==0) return X;
		else if(index==1) return Y;
		else if(index==2) return Z;
		else return 0;
	} 
	inline float v(unsigned int i,float value) {
		if(i==0) X=value;
		else if(i==1) Y=value;
		else if(i==2) Z=value;
	} 
	V3(float x,float y,float z) {
	 X=x;Y=y;Z=z;
	}
	inline float Len() {
		return sqrt(X*X+Y*Y+Z*Z);
	}
	inline void Norm() {
		float n;
		n=Len();
		if (n==0) return;
		X/=n;
		Y/=n;
		Z/=n;
	}
} ;
																								
/* vector arithmetic */
 
struct V4 {
	float X,Y,Z,W;
	inline float v(unsigned int index) {
		if(index==0) return X;
		else if(index==1) return Y;
		else if(index==2) return Z;
		else if(index==3) return W;
		else return 0;
	} 
	inline float v(unsigned int i,float value) {
		if(i==0) X=value;
		else if(i==1) Y=value;
		else if(i==2) Z=value;
		else if(i==3) W=value;
	} 
	V4() {
		X=0;Y=0;Z=0;W=0;
	}
	V4(float x,float y,float z,float w) {
	 X=x;Y=y;Z=z;W=w;
	}
} ;
struct COLOR4 {
	float R,G,B,A;
	COLOR4() {
		R=0;G=0;B=0;A=0;
	}
	COLOR4(float r,float g,float b,float a) {
	 R=r;G=g;B=b;A=a;
	}
};
	 
 


/* ******* Gestion des matrices 4x4 ****** */


static inline void gl_M4_Mul(M4 *c,M4 *a,M4 *b)
{
  int i,j,k;
  float s;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++) {
      s=0.0;
      for(k=0;k<4;k++) s+=a->m[i][k]*b->m[k][j];
      c->m[i][j]=s;
    }
}

/* c=c*a */
static inline void gl_M4_MulLeft(M4 *c,M4 *b)
{
  int i,j,k;
  float s;
  M4 a;

  /*memcpy(&a, c, 16*sizeof(float));
  */
  a=*c;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++) {
      s=0.0;
      for(k=0;k<4;k++) s+=a.m[i][k]*b->m[k][j];
      c->m[i][j]=s;
    }
}

static inline void gl_M4_Move(M4 *a,M4 *b)
{
	memcpy(a,b,sizeof(M4));
}

static inline void gl_MoveV3(V3 *a,V3 *b)
{
	memcpy(a,b,sizeof(V3));
}


static inline void gl_MulM4V3(V3 *a,M4 *b,V3 *c)
{
	 a->X=b->m[0][0]*c->X+b->m[0][1]*c->Y+b->m[0][2]*c->Z+b->m[0][3];
	 a->Y=b->m[1][0]*c->X+b->m[1][1]*c->Y+b->m[1][2]*c->Z+b->m[1][3];
	 a->Z=b->m[2][0]*c->X+b->m[2][1]*c->Y+b->m[2][2]*c->Z+b->m[2][3];
}

static inline void gl_MulM3V3(V3 *a,M4 *b,V3 *c)
{
	 a->X=b->m[0][0]*c->X+b->m[0][1]*c->Y+b->m[0][2]*c->Z;
	 a->Y=b->m[1][0]*c->X+b->m[1][1]*c->Y+b->m[1][2]*c->Z;
	 a->Z=b->m[2][0]*c->X+b->m[2][1]*c->Y+b->m[2][2]*c->Z;
}

static inline void gl_M4_MulV4(V4 *a,M4 *b,V4 *c)
{
	 a->X=b->m[0][0]*c->X+b->m[0][1]*c->Y+b->m[0][2]*c->Z+b->m[0][3]*c->W;
	 a->Y=b->m[1][0]*c->X+b->m[1][1]*c->Y+b->m[1][2]*c->Z+b->m[1][3]*c->W;
	 a->Z=b->m[2][0]*c->X+b->m[2][1]*c->Y+b->m[2][2]*c->Z+b->m[2][3]*c->W;
	 a->W=b->m[3][0]*c->X+b->m[3][1]*c->Y+b->m[3][2]*c->Z+b->m[3][3]*c->W;
}
	
/* transposition of a 4x4 matrix */
static inline void gl_M4_Transpose(M4 *a,M4 *b)
{
  a->m[0][0]=b->m[0][0]; 
  a->m[0][1]=b->m[1][0]; 
  a->m[0][2]=b->m[2][0]; 
  a->m[0][3]=b->m[3][0]; 

  a->m[1][0]=b->m[0][1]; 
  a->m[1][1]=b->m[1][1]; 
  a->m[1][2]=b->m[2][1]; 
  a->m[1][3]=b->m[3][1]; 

  a->m[2][0]=b->m[0][2]; 
  a->m[2][1]=b->m[1][2]; 
  a->m[2][2]=b->m[2][2]; 
  a->m[2][3]=b->m[3][2]; 

  a->m[3][0]=b->m[0][3]; 
  a->m[3][1]=b->m[1][3]; 
  a->m[3][2]=b->m[2][3]; 
  a->m[3][3]=b->m[3][3]; 
}

/* inversion of an orthogonal matrix of type Y=M.X+P */ 
static inline void gl_M4_InvOrtho(M4 *a,M4 b)
{
	int i,j;
	float s;
	for(i=0;i<3;i++)
	for(j=0;j<3;j++) a->m[i][j]=b.m[j][i];
	a->m[3][0]=0.0; a->m[3][1]=0.0; a->m[3][2]=0.0; a->m[3][3]=1.0;
	for(i=0;i<3;i++) {
		s=0;
		for(j=0;j<3;j++) s-=b.m[j][i]*b.m[j][3];
		a->m[i][3]=s;
	}
}

/* Inversion of a general nxn matrix.
   Note : m is destroyed */

static inline int Matrix_Inv(float *r,float *m,int n)
{
	 int i,j,k,l;
	 float max,tmp,t;

	 /* identit�e dans r */
	 for(i=0;i<n*n;i++) r[i]=0;
	 for(i=0;i<n;i++) r[i*n+i]=1;
	 
	 for(j=0;j<n;j++) {
			
			/* recherche du nombre de plus grand module sur la colonne j */
			max=m[j*n+j];
			k=j;
			for(i=j+1;i<n;i++)
				if (fabs(m[i*n+j])>fabs(max)) {
					 k=i;
					 max=m[i*n+j];
				}

      /* non intersible matrix */
      if (max==0) return 1;

			
			/* permutation des lignes j et k */
			if (k!=j) {
				 for(i=0;i<n;i++) {
						tmp=m[j*n+i];
						m[j*n+i]=m[k*n+i];
						m[k*n+i]=tmp;
						
						tmp=r[j*n+i];
						r[j*n+i]=r[k*n+i];
						r[k*n+i]=tmp;
				 }
			}
			
			/* multiplication de la ligne j par 1/max */
			max=1/max;
			for(i=0;i<n;i++) {
				 m[j*n+i]*=max;
				 r[j*n+i]*=max;
			}
			
			
			for(l=0;l<n;l++) if (l!=j) {
				 t=m[l*n+j];
				 for(i=0;i<n;i++) {
						m[l*n+i]-=m[j*n+i]*t;
						r[l*n+i]-=r[j*n+i]*t;
				 }
			}
	 }

	 return 0;
}


/* inversion of a 4x4 matrix */

static inline void gl_M4_Inv(M4 *a,M4 *b)
{
  M4 tmp;
  memcpy(&tmp, b, 16*sizeof(float));
  /*tmp=*b;*/
  Matrix_Inv(&a->m[0][0],&tmp.m[0][0],4);
}

static inline void gl_M4_Rotate(M4 *a,float t,int u)
{
	 float s,c;
	 int v,w;
   if ((v=u+1)>2) v=0;
	 if ((w=v+1)>2) w=0;
	 s=sin(t);
	 c=cos(t);
	 a->Id();
	 a->m[v][v]=c;	a->m[v][w]=-s;
	 a->m[w][v]=s;	a->m[w][w]=c;
}
	

/* inverse of a 3x3 matrix */
static inline void gl_M3_Inv(M3 *a,M3 *m)
{
	 float det;
	 
	 det = m->m[0][0]*m->m[1][1]*m->m[2][2]-m->m[0][0]*m->m[1][2]*m->m[2][1]-
		 m->m[1][0]*m->m[0][1]*m->m[2][2]+m->m[1][0]*m->m[0][2]*m->m[2][1]+
		 m->m[2][0]*m->m[0][1]*m->m[1][2]-m->m[2][0]*m->m[0][2]*m->m[1][1];

	 a->m[0][0] = (m->m[1][1]*m->m[2][2]-m->m[1][2]*m->m[2][1])/det;
	 a->m[0][1] = -(m->m[0][1]*m->m[2][2]-m->m[0][2]*m->m[2][1])/det;
	 a->m[0][2] = -(-m->m[0][1]*m->m[1][2]+m->m[0][2]*m->m[1][1])/det;
	 
	 a->m[1][0] = -(m->m[1][0]*m->m[2][2]-m->m[1][2]*m->m[2][0])/det;
	 a->m[1][1] = (m->m[0][0]*m->m[2][2]-m->m[0][2]*m->m[2][0])/det;
	 a->m[1][2] = -(m->m[0][0]*m->m[1][2]-m->m[0][2]*m->m[1][0])/det;

	 a->m[2][0] = (m->m[1][0]*m->m[2][1]-m->m[1][1]*m->m[2][0])/det;
	 a->m[2][1] = -(m->m[0][0]*m->m[2][1]-m->m[0][1]*m->m[2][0])/det;
	 a->m[2][2] = (m->m[0][0]*m->m[1][1]-m->m[0][1]*m->m[1][0])/det;
}

		


#endif  /* __ZMATH__ */
