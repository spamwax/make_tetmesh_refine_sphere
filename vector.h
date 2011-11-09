#ifndef __vector_header
#define __vector_header

#include "init.h"

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180./PI)
#endif

void v_make(const double *pnt1, const double *pnt2, int n, double *vec) ; //      subroutine v_make(pnt1, pnt2, n, vec)
int v_norm(double *vec, int n) ; //      subroutine v_norm(vec, n)
void v_cros(const double *vec1, const double *vec2, int n, double *prod) ; //      subroutine v_cros( vec1, vec2, n, prod )
double v_magn(const double *vec, int n) ; //      real function v_magn(vec, n)
double v_dot(const double *vec1, const double *vec2, int n) ;  //    real function v_dot(vec1,vec2, n)
double v_rang(const double *vec0, const double *vec1, int n) ;  //      real function v_rang( vec0, vec1, n )
double *v_add(const double *v1, const double *v2, int n);
double v_dist( const double x[], const double y[], int n ); /* Caclculates the distance between x and y in n-dimension */
double v_angle (double vec1[], double vec2[], int degf, int dim, double ZeroTol = FloatTol);
double v_angle3(double vec1[], double vec2[], int degf, int dim, double ZeroTol = FloatTol);
double angle_rad_2d(double p1[], double p2[], double p3[], int degf, double zeroTol); // Returns the angle swept between rays measured from p1-p2 to p3-p2 ccw
#endif
