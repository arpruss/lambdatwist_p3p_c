/*
P3P
Copyright (C) 2021 Mikael Persson, Klas Nordberg and Alexander Pruss

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>

typedef double Vector2[2];
typedef double Vector3[3];
typedef double Matrix33[3][3];

#define NUMERIC_LIMIT 1e-13
#define KLAS_P3P_CUBIC_SOLVER_ITER 50

static int root2double(double b, double c,double* r1p, double* r2p){
    double v=b*b -4*c;
    if(v<0){
        *r1p=*r2p=0.5*b;
        return 0;
    }

    double y=sqrt(v);
    if(b<0){
        *r1p= 0.5*(-b+y);
        *r2p= 0.5*(-b-y);
    }else{
        *r1p= 2.0*c/(-b+y);
        *r2p= 2.0*c/(-b-y);
    }
    return 1;
}


static double cubick(double b, double c, double d){
    /* Choose initial solution */
    double r0;
    
    // not monotonic
    if (b*b  >= 3*c){
        // h has two stationary points, compute them
        //T t1 = t - std::sqrt(diff);
        double v=sqrt(b*b -3*c);
        double t1 = (-b - v)/3;

        // Check if h(t1) > 0, in this case make a 2-order approx of h around t1
        double k = ((t1+b)*t1+c)*t1+d;

        if (k > 0) {
            //Find leftmost root of 0.5*(r0 -t1)^2*(6*t1+2*b) +  k = 0
            r0 = t1 - sqrt(-k/(3*t1 + b));
            // or use the linear comp too
            //r0=t1 -
        } else {
            double t2 = (-b + v)/3;
            k = ((t2+b)*t2+c)*t2+d;
            //Find rightmost root of 0.5*(r0 -t2)^2*(6*t2+2*b) +  k1 = 0
            r0 = t2 + sqrt(-k/(3*t2 + b));
        }
    }
    else{
        r0=-b/3;
        if(fabs(((3*r0+2*b)*r0+c))<1e-4) r0+=1;
    }

    /* Do ITER Newton-Raphson iterations */
    /* Break if position of root changes less than 1e.16 */
    //T starterr=std::abs(r0*(r0*(r0 + b) + c) + d);
    double fx,fpx;

    for (unsigned int cnt = 0; cnt < KLAS_P3P_CUBIC_SOLVER_ITER; ++cnt){

        //(+ (* r0 (+  c (* (+ r0 b) r0) )) d )
        fx=(((r0+b)*r0+c)*r0+d);


        if((cnt<7 || fabs(fx)>NUMERIC_LIMIT)  ){
            fpx=((3*r0+2*b)*r0+c);

            r0-= fx/fpx;
        }
        else
            break;
    }
    //timer_b.toc();

    return r0;
}

static void invert3(Matrix33 M, double a00, double a01, double a02,
                                double a10, double a11, double a12,
                                double a20, double a21, double a22) {
    double idet; // Determinant

    M[0][0] = a11 * a22 - a12 * a21;
    M[0][1] = a02 * a21 - a01 * a22;
    M[0][2] = a01 * a12 - a02 * a11;

    M[1][0] = a12 * a20 - a10 * a22;
    M[1][1] = a00 * a22 - a02 * a20;
    M[1][2] = a02 * a10 - a00 * a12;

    M[2][0] = a10 * a21 - a11 * a20;
    M[2][1] = a01 * a20 - a00 * a21;
    M[2][2] = a00 * a11 - a01 * a10;

    idet =1/(a00 * M[0][0] + a01 * M[1][0] + a02 * M[2][0]);
    
    for(int i=0;i<3;i++) for(int j=0;j<3;j++)
        M[i][j] *= idet;
}

static inline double max(double a,double b) {
    return a<b ? b : a;
}

static inline double normSquared3(Vector3 v) {
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

static inline void normalize3(Vector3 out,Vector3 in) {
    double n = 1/sqrt(normSquared3(in));
    out[0] = in[0] * n;
    out[1] = in[1] * n;
    out[2] = in[2] * n;
}

static inline void homogenize23(Vector3 out,Vector2 in) {
    double n = 1/sqrt(in[0]*in[0]+in[1]*in[1]+1);
    out[0] = in[0] * n;
    out[1] = in[1] * n;
    out[2] = n;
}


static inline void set3(Vector3 out,double x0,double x1,double x2) {
    out[0] = x0;
    out[1] = x1;
    out[2] = x2;
}

static inline void sub3(Vector3 out,Vector3 a,Vector3 b) {
    out[0] = a[0]-b[0];
    out[1] = a[1]-b[1];
    out[2] = a[2]-b[2];
}

static inline void add3(Vector3 out,Vector3 a,Vector3 b) {
    out[0] = a[0]+b[0];
    out[1] = a[1]+b[1];
    out[2] = a[2]+b[2];
}

static inline double dot3(Vector3 a,Vector3 b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

static inline void cross3(Vector3 out, Vector3 a, Vector3 b) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];    
}

static inline void scale3(Vector3 out, double a, Vector3 b) {
    out[0] = a * b[0];
    out[1] = a * b[1];
    out[2] = a * b[2];
}

static void mult33(Matrix33 m, Matrix33 a, Matrix33 b) {
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            m[i][j] = a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
}

static void apply3(Vector3 out, Matrix33 m, Vector3 a) {
    for (int i=0;i<3;i++)
        out[i] = m[i][0]*a[0]+m[i][1]*a[1]+m[i][2]*a[2];
}

static void apply3T(Vector3 out, Matrix33 m, Vector3 a) {
    for (int i=0;i<3;i++)
        out[i] = m[0][i]*a[0]+m[1][i]*a[1]+m[2][i]*a[2];
}

/**
 * @brief refineL
 * @param L
 * @param a12
 * @param a13
 * @param a23
 * @param b12
 * @param b13
 * @param b23
 *
 * Gauss-Newton Solver
 * For unknown reasons it always works for the correct solution, but not always for the other solutions!
 *
 */
static void gauss_newton_refineL(Vector3 L,
                          double a12, double a13, double a23,
                          double b12, double b13, double b23,
                          int iterations
                          ){

    // const expr makes it easier for the compiler to unroll
    for(int i=0;i<iterations;++i){
        double l1=L[0];
        double l2=L[1];
        double l3=L[2];
        double r1=l1*l1 + l2*l2 +b12*l1*l2 -a12;
        double r2=l1*l1 + l3*l3 +b13*l1*l3 -a13;
        double r3=l2*l2 + l3*l3 +b23*l2*l3 -a23;

        if(fabs(r1) +fabs(r2) +fabs(r3)<1e-10) break;

        double dr1dl1=(2.0)*l1 +b12*l2;
        double dr1dl2=(2.0)*l2 +b12*l1;
        //T dr1dl3=0;

        double dr2dl1=(2.0)*l1 +b13*l3;
        //T dr2dl2=0;
        double dr2dl3=(2.0)*l3 +b13*l1;


        //T dr3dl1=0;
        double dr3dl2=(2.0)*l2 + b23*l3;
        double dr3dl3=(2.0)*l3 + b23*l2;

        Vector3 r;
        set3(r, r1, r2, r3);

        // or skip the inverse and make it explicit...
        {

            double v0=dr1dl1;
            double v1=dr1dl2;
            double v3=dr2dl1;
            double v5=dr2dl3;
            double v7=dr3dl2;
            double v8=dr3dl3;
            double det=(1.0)/(- v0*v5*v7 - v1*v3*v8);

            Matrix33 Ji = { { -v5*v7, -v1*v8,  v1*v5 },
                            { -v3*v8,  v0*v8, -v0*v5 },
                            { v3*v7, -v0*v7, -v1*v3 } };
            Vector3 L1;
            Vector3 z;
            apply3(z,Ji,r);
            scale3(z,det,z);
            sub3(L,L,z);
            //L=Vector3<T>(L) - det*(Ji*r);
            //%l=l - g*H\G;%inv(H)*G
            //L=L - g*J\r; //% works because the size is ok!

            {
                double l1=L1[0];
                double l2=L1[1];
                double l3=L1[2];
                double r11=l1*l1 + l2*l2 +b12*l1*l2 -a12;
                double r12=l1*l1 + l3*l3 +b13*l1*l3 -a13;
                double r13=l2*l2 + l3*l3 +b23*l2*l3 -a23;
                if(fabs(r11) +fabs(r12) + fabs(r13)>fabs(r1) +fabs(r2) +fabs(r3)){
                    // cout<<"bad step: "<< det*(Ji*r)<<++badsteps<<" "<< i<<endl;
                    break;
                }
                else
                    set3(L,L1[0],L1[1],L1[2]);
            }
        }
    }
    // cout<<i<<endl;
}

static void eigwithknown330(Matrix33 x, Matrix33 E, Vector3 L) {
    // one eigenvalue is known to be 0.

    //the known one...
    L[2]=0;

    Vector3 v3;
    set3(v3, x[1][0]*x[2][1]- x[2][0]*x[1][1],
                   x[2][0]*x[0][1]- x[2][1]*x[0][0],
                   x[1][1]*x[0][0]- x[1][0]*x[0][1]);
    double length = sqrt(normSquared3(v3));
    v3[0] /= length;
    v3[1] /= length;
    v3[2] /= length;

    double x01_squared=x[0][1]*x[0][1];
    // get the two other...
    double b=- x[0][0] - x[1][1] - x[2][2];
    double c=- x01_squared - x[0][2]*x[0][2] - x[1][2]*x[1][2] +
            x[0][0]*(x[1][1] + x[2][2]) + x[1][1]*x[2][2];
    double e1,e2;
    //roots(poly(x))
    root2double(b,c,&e1,&e2);

    if(fabs(e1)<fabs(e2)) {
        double t = e1;
        e1 = e2;
        e2 = t;
    }
    L[0]=e1;
    L[1]=e2;

    double mx0011=-x[0][0]*x[1][1];
    double prec_0 = x[0][1]*x[1][2] - x[0][2]*x[1][1];
    double prec_1 = x[0][1]*x[0][2] - x[0][0]*x[1][2];

    double e=e1;
    double tmp=1/(e*(x[0][0] + x[1][1]) + mx0011 - e*e + x01_squared);
    double a1= -(e*x[0][2] + prec_0)*tmp;
    double a2= -(e*x[1][2] + prec_1)*tmp;
    double rnorm=1/sqrt(a1*a1 +a2*a2 + 1);
    a1*=rnorm;
    a2*=rnorm;
    
    Vector3 v1;
    set3(v1,a1,a2,rnorm);

    //e=e2;
    double tmp2=1/(e2*(x[0][0] + x[1][1]) + mx0011 - e2*e2 + x01_squared);
    double a21= -(e2*x[0][2] + prec_0)*tmp2;
    double a22= -(e2*x[1][2] + prec_1)*tmp2;
    double rnorm2=1/sqrt(a21*a21 +a22*a22 +1);
    a21*=rnorm2;
    a22*=rnorm2;
    Vector3 v2;
    set3(v2,a21,a22,rnorm2);


    // optionally remove axb from v1,v2
    // costly and makes a very small difference!
    // v1=(v1-v1.dot(v3)*v3);v1.normalize();
    // v2=(v2-v2.dot(v3)*v3);v2.normalize();
    // v2=(v2-v1.dot(v2)*v2);v2.normalize();

    for (int i=0;i<3;i++) {
        E[i][0] = v1[i];
        E[i][1] = v2[i];
        E[i][2] = v3[i];
    }
}

int p3p_lambdatwist(Vector2 _y1, Vector2 _y2, Vector2 _y3, Vector3 x1, Vector3 x2, Vector3 x3,
        Matrix33* Rs, Vector3* Ts, int refinement_iterations) {
            
    Vector3 y1,y2,y3;
    
    homogenize23(y1,_y1);
    homogenize23(y2,_y2);
    homogenize23(y3,_y3);
    
    double b12=-2*dot3(y1,y2);
    double b13=-2*dot3(y1,y3);
    double b23=-2*dot3(y2,y3);

    Vector3 d12;
    Vector3 d13;
    Vector3 d23;
    Vector3 d12xd13;
    
    sub3(d12,x1,x2);
    sub3(d13,x1,x3);
    sub3(d23,x2,x3);
    
    cross3(d12xd13,d12,d13);
    
    double a12 = normSquared3(d12);
    double a13 = normSquared3(d13);
    double a23 = normSquared3(d23);

    //a*g^3 + b*g^2 + c*g + d = 0
    double c31=-0.5*b13;
    double c23=-0.5*b23;
    double c12=-0.5*b12;
    double blob=c12*c23*c31-1;

    double s31_squared=1-c31*c31;
    double s23_squared=1-c23*c23;
    double s12_squared=1-c12*c12;

    double p3 = a13*(a23*s31_squared - a13*s23_squared);

    double p2 = 2*blob*a23*a13 + a13*(2*a12 + a13)*s23_squared + a23*(a23 - a12)*s31_squared;

    double p1 = a23*(a13 - a23)*s12_squared - a12*a12*s23_squared - 2*a12*(blob*a23 + a13*s23_squared);

    double p0 = a12*(a12*s23_squared - a23*s12_squared);

    double g=0;

    //p3 is detD2 so its definietly >0 or its a degenerate case

    {
        p3=1/p3;
        p2*=p3;
        p1*=p3;
        p0*=p3;

        // get sharpest double root of above...

        g=cubick(p2,p1,p0);
    }

    // we can swap D1,D2 and the coeffs!
    // oki, Ds are:
    //D1=M12*XtX(2,2) - M23*XtX(1,1);
    //D2=M23*XtX(3,3) - M13*XtX(2,2);

    //[    a23 - a23*g,                 (a23*b12)/2,              -(a23*b13*g)/2]
    //[    (a23*b12)/2,           a23 - a12 + a13*g, (a13*b23*g)/2 - (a12*b23)/2]
    //[ -(a23*b13*g)/2, (a13*b23*g)/2 - (a12*b23)/2,         g*(a13 - a23) - a12]

    Matrix33 A;
    A[0][0] =        a23*(1-g);
    A[0][1] =         a23*b12*0.5;
    A[0][2] =        a23*b13*g*-0.5;
    A[1][1] =    a23 - a12 + a13*g;
    A[1][2] =        b23*(a13*g - a12)*0.5;
    A[2][2] =  g*(a13 - a23) - a12;
    A[1][0] = A[0][1];
    A[2][0] = A[0][2];
    A[2][1] = A[1][2];

    // get sorted eigenvalues and eigenvectors given that one should be zero...
    Matrix33 V;
    Vector3 L;
    eigwithknown330(A,V,L);

    double v=sqrt(max(0,-L[1]/L[0]));

    int valid=0;
    Vector3 Ls[4];    

    // use the t=Vl with t2,st2,t3 and solve for t3 in t2
    { //+v

        double s=v;
        double w2=1/(s*V[0][1] - V[0][0]);
        double w0=(V[1][0] - s*V[1][1])*w2;
        double w1=(V[2][0] - s*V[2][1])*w2;

        double a=1/((a13 - a12)*w1*w1 - a12*b13*w1 - a12);
        double b=(a13*b12*w1 - a12*b13*w0 - 2*w0*w1*(a12 - a13))*a;
        double c=((a13 - a12)*w0*w0 + a13*b12*w0 + a13)*a;

        if (b*b-4*c>=0) {
            double tau1,tau2;
            root2double(b,c,&tau1,&tau2);
            if (tau1>0){
                double tau=tau1;
                double d=a23/(tau*(b23 + tau) + 1);

                double l2=sqrt(d);
                double l3=tau*l2;

                double l1=w0*l2 +w1*l3;
                if(l1>=0){

                    set3(Ls[valid],l1,l2,l3);

                    ++valid;
                }

            }
            if(tau2>0){
                double tau=tau2;
                double d=a23/(tau*(b23 + tau) + 1);

                double l2=sqrt(d);
                double l3=tau*l2;
                double l1=w0*l2 +w1*l3;
                if(l1>=0){
                    set3(Ls[valid],l1,l2,l3);
                    ++valid;
                }

            }
        }
    }

    { //+v
        double s=-v;
        double w2=1/(s*V[0][1] - V[0][0]);
        double w0=(V[1][0] - s*V[1][1])*w2;
        double w1=(V[2][0] - s*V[2][1])*w2;

        double a=1/((a13 - a12)*w1*w1 - a12*b13*w1 - a12);
        double b=(a13*b12*w1 - a12*b13*w0 - 2*w0*w1*(a12 - a13))*a;
        double c=((a13 - a12)*w0*w0 + a13*b12*w0 + a13)*a;


        if(b*b-4*c>=0) {
            double tau1,tau2;

            root2double(b,c,&tau1,&tau2);
            if (tau1>0) {
                double tau=tau1;
                double d=a23/(tau*(b23 + tau) + 1);
                if(d>0){
                  double l2=sqrt(d);

                  double l3=tau*l2;

                  double l1=w0*l2 + w1*l3;
                  if(l1>=0){
                      set3(Ls[valid],l1,l2,l3);
                      ++valid;
                  }
                }
            }
            if(tau2>0){
                double tau=tau2;
                double d=a23/(tau*(b23 + tau) + 1);
                if(d>0){
                  double l2=sqrt(d);

                  double l3=tau*l2;

                  double l1=w0*l2 +w1*l3;
                  if(l1>=0){
                      set3(Ls[valid],l1,l2,l3);
                      ++valid;
                  }
                }
            }
        }
    }

    for(int i=0;i<valid;++i) {              
        gauss_newton_refineL(Ls[i],a12,a13,a23,b12,b13,b23,
            refinement_iterations);        
    }

    Vector3 ry1,ry2,ry3;
    Vector3 yd1;
    Vector3 yd2;
    Vector3 yd1xd2;
    Matrix33 X;
    
    invert3(X,d12[0],d13[0],d12xd13[0],
                    d12[1],d13[1],d12xd13[1],
                    d12[2],d13[2],d12xd13[2]);

    for(int i=0;i<valid;++i){
        // compute the rotation:
        scale3(ry1,Ls[i][0],y1);
        scale3(ry2,Ls[i][1],y2);
        scale3(ry3,Ls[i][2],y3);

        sub3(yd1,ry1,ry2);
        sub3(yd2,ry1,ry3);
        cross3(yd1xd2,yd1,yd2);

        Matrix33 Y;
        for(int j=0;j<3;j++) {
            Y[j][0] = yd1[j];
            Y[j][1] = yd2[j];
            Y[j][2] = yd1xd2[j];
        }

        mult33(Rs[i],Y,X);
        
        for(int j=0;j<3;j++) {
            Ts[i][j] = ry1[j] - dot3(Rs[i][j],x1); 
        }
    }

    return valid;
}

void applyRT(Vector2 out, Matrix33 R, Vector3 T, Vector3 in) {
    Vector3 temp;
    apply3(temp, R, in);
    add3(temp, T, temp);
    out[0] = temp[0]/temp[2];
    out[1] = temp[1]/temp[2];    
}

void toHomography(double* h, Matrix33 R, Vector3 T) {
    h[0] = R[0][0]/T[2];
    h[1] = R[0][1]/T[2];
    h[2] = T[0]/T[2];
    h[3] = R[1][0]/T[2];
    h[4] = R[1][1]/T[2];
    h[5] = T[1]/T[2];
    h[6] = R[2][0]/T[2];
    h[7] = R[2][1]/T[2];
}

void applyHomography(Vector2 out, double* h, Vector2 in) {
    double scale = 1/(in[0]*h[6]+in[1]*h[7]+1);
    
    out[0] = (h[0]*in[0]+h[1]*in[1]+h[2])*scale;
    out[1] = (h[3]*in[0]+h[4]*in[1]+h[5])*scale;
}

#ifdef TEST
#include <stdio.h>

void showVect2(Vector2 v) {
    printf("%lg %lg\n", v[0],v[1]);
}

void showVect3(Vector3 v) {
    printf("%lg %lg %lg\n", v[0],v[1],v[2]);
}

void main() {
    Vector2 y1 = { 0, 0 };
    Vector2 y2 = { 0, 7 };
    Vector2 y3 = { 10, 0 };
    Vector3 x1 = { 0, 0, 0 };
    Vector3 x2 = { 0, 20, 0 };
    Vector3 x3 = { 20, 0, 0 };
    Matrix33 Rs[4];
    Vector3 Ts[4];

    int valid = p3p_lambdatwist(y1,y2,y3,x1,x2,x3,Rs,Ts,5);
    if (valid) {
        for (int i=0; i<valid; i++) {
            printf("test %d \n",i);
            showVect3(Ts[i]);
            for (int j=0; j<3; j++)
                showVect3(Rs[i][j]);
            Vector2 out;
            double h[8];
            toHomography(h, Rs[i], Ts[i]);
            printf("x1->y1 ");
            applyRT(out, Rs[i], Ts[i], x1);
            showVect2(out);
            applyHomography(out, h, x1);
            showVect2(out);
            printf("x2->y2 ");
            applyRT(out, Rs[i], Ts[i], x2);
            showVect2(out);
            applyHomography(out, h, x2);
            showVect2(out);
            printf("x3->y3 ");
            applyRT(out, Rs[i], Ts[i], x3);
            showVect2(out);
            applyHomography(out, h, x3);
            showVect2(out);
        }
    }
}
#endif
