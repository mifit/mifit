#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <cmath>
#include <vector>

#include "Matrix.h"

#ifndef M_PI
#define M_PI (3.1415926535897932384626433832795028841972)
#endif

namespace chemlib
{

    inline double DotVect(const double in1[3], const double in2[3])
    {
        return ( in1[0]*in2[0] + in1[1]*in2[1] + in1[2]*in2[2] );
    }

    inline double DotVect(const double lhx,
                          const double lhy,
                          const double lhz,
                          const double rhx,
                          const double rhy,
                          const double rhz)
    {
        return ( lhx*rhx + lhy*rhy + lhz*rhz );
    }

    inline double CosVectorAngle(double v1[3], double v2[3])
    {
        double d1d2;

        d1d2 = (double)sqrt((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])
                            *(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]));

        return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) / d1d2;
    }

    inline void CrossVects(double a[3], double b[3], double c[3])
    {
        /* c := cross product of a and b */
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
    }

//Cross two vectors, then dot the cross product with a
//third vector and return the dot product
    inline float Cross_and_dot_3D(float *v0, float *v1, float *v2)
    {
        return (v0[1] * v1[2] - v0[2] * v1[1]) * v2[0]
               +(v0[2] * v1[0] - v0[0] * v1[2]) * v2[1]
               +(v0[0] * v1[1] - v0[1] * v1[0]) * v2[2];
    }

    inline float Cross_and_dot_3D(double *v0, double *v1, double *v2)
    {
        return (float)((v0[1] * v1[2] - v0[2] * v1[1]) * v2[0]
                       +(v0[2] * v1[0] - v0[0] * v1[2]) * v2[1]
                       +(v0[0] * v1[1] - v0[1] * v1[0]) * v2[2]);
    }

    template <typename T>
    inline void ScaleVect(std::vector<T> &v, double scale_factor)
    {
        v[0] = scale_factor * v[0];
        v[1] = scale_factor * v[1];
        v[2] = scale_factor * v[2];
    }

    template <typename T>
    inline void ScaleVect(T *v, double scale_factor)
    {
        v[0] = (T)(scale_factor * v[0]);
        v[1] = (T)(scale_factor * v[1]);
        v[2] = (T)(scale_factor * v[2]);
    }

    template <typename T>
    inline void ScaleVect2D(T *v, double scale_factor)
    {
        v[0] = scale_factor * v[0];
        v[1] = scale_factor * v[1];
    }

    template <typename T>
    inline double VectLength(const std::vector<T> &v)
    {
        return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

    template <typename T>
    inline double VectLength(const T *v)
    {
        return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

    template <typename T>
    inline double VectLength2D(const T *v)
    {
        return sqrt(v[0]*v[0] + v[1]*v[1]);
    }

    template <typename T>
    inline double SquaredVectLength(const std::vector<T> &v)
    {
        return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    }

    template <typename T>
    inline double SquaredVectLength(const T *v)
    {
        return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    }

    template <typename T>
    inline void NormVect(T *v)
    {
        ScaleVect(v, 1 / VectLength(v));
    }

    template <typename T>
    inline void NormVect(const T *in, T *out)
    {
        double normfactor;
        normfactor = 1./VectLength(in);

        out[0] = in[0] * normfactor;
        out[1] = in[1] * normfactor;
        out[2] = in[2] * normfactor;
    }

    template <typename T>
    inline void NormVect(std::vector<T> &v)
    {
        ScaleVect(v, 1 / VectLength(v));
    }

/*
    inline void NormVect(double in[3], double out[3]) {
        double normfactor;
        normfactor = 1./VectLength(in);

        out[0] = in[0] * normfactor;
        out[1] = in[1] * normfactor;
        out[2] = in[2] * normfactor;
    }
 */
    template <typename T>
    inline void NormVect2D(T *v)
    {
        ScaleVect2D(v, 1 / VectLength2D(v));
    }

/////////////////////////////////////////////////////////////////////////////
// Function:    AzimuthAngle
// Purpose:		Measures the angle (in radians) necessary to convert the given
//				x-y point to polar/cylindrical coordinates
// Input:       x and y values of a point
// Output:      A value between -PI and PI...the rotation angle in radians
// Requires:
/////////////////////////////////////////////////////////////////////////////
    inline double AzimuthAngle(double x, double y)
    {
        if (x == 0 && y == 0)
        {
            return 0;
        }
        else
        {
            return atan2(y, x);
        }
    }

/////////////////////////////////////////////////////////////////////////////
// Function:    TriangleAngle
// Purpose:		Applies the law of cosines to calculate an angle (in radians)
//				within a triangle, given the lengths of the three sides of the
//				triangle
// Input:       The lengths of the three sides
// Output:      A value between 0 and PI...the angle in radians
// Requires:	The side opposite the desired angled should be the first argument
/////////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline double TriangleAngle(T a, T b, T c)
    {

        double cos_alpha = ( b * b + c * c - a * a ) / ( 2 * b * c );

        if (cos_alpha > 1 || cos_alpha < -1)
        {
            std::string err = "Error in TriangleAngle, triangle inequality violation";
            throw err;
        }

        return acos(cos_alpha);
    }

    template <typename T>
    inline double CalcGradNorm(T *u, T *v, T *w)
    {
        T uv01 = u[0]*v[1] - u[1]*v[0];
        T uv02 = u[0]*v[2] - u[2]*v[0];
        T uv12 = u[1]*v[2] - u[2]*v[1];
        T uw01 = u[0]*w[1] - u[1]*w[0];
        T uw02 = u[0]*w[2] - u[2]*w[0];
        T uw12 = u[1]*w[2] - u[2]*w[1];
        T vw01 = v[0]*w[1] - v[1]*w[0];
        T vw02 = v[0]*w[2] - v[2]*w[0];
        T vw12 = v[1]*w[2] - v[2]*w[1];

        return SQR(uv01) + SQR(uv02) + SQR(uv12)
               +SQR(uw01) + SQR(uw02) + SQR(uw12)       //These three lines summarize the 3 external pts
               +SQR(vw01) + SQR(vw02) + SQR(vw12)
               + //				SQR(-uv12 + uw12 - vw12) +
                 //				SQR(uv02 - uw02 + vw02) +				//These three lines summarize the center pt
                 //				SQR(-uv01 + uw01 - vw01);
               SQR(-uv12 + uw12 - vw12)
               +SQR(uv02 - uw02 + vw02)                 //These three lines summarize the center pt
               +SQR(-uv01 + uw01 - vw01);
    }

    inline double CalcTetraVolume(double d12, double d13, double d14, double d23,
                                  double d24, double d34)
    {

        TNT::Matrix<double> mat(5, 5);

        mat[0][0] = 0;
        mat[0][1] = 1;
        mat[0][2] = 1;
        mat[0][3] = 1;
        mat[0][4] = 1;
        mat[1][0] = 1;
        mat[1][1] = 0;
        mat[1][2] = d12 * d12;
        mat[1][3] = d13 * d13;
        mat[1][4] = d14 * d14;
        mat[2][0] = 1;
        mat[2][1] = mat[1][2];
        mat[2][2] = 0;
        mat[2][3] = d23 * d23;
        mat[2][4] = d24 * d24;
        mat[3][0] = 1;
        mat[3][1] = mat[1][3];
        mat[3][2] = mat[2][3];
        mat[3][3] = 0;
        mat[3][4] = d34 * d34;
        mat[4][0] = 1;
        mat[4][1] = mat[1][4];
        mat[4][2] = mat[2][4];
        mat[4][3] = mat[3][4];
        mat[4][4] = 0;

        double det = TNT::Determinant(mat);

        if (det >= 0)
        {
            return sqrt(det / 288.0);
        }
        else
        {
            throw "Error in CalcTetraVolume.  Input distances are inconsistent with a tetrahedron in three dimensions";
            //			return sqrt( -det / 288.0);
        }
    }

    inline double CalcTetraVolume(std::vector<double> &dists)
    {
        if (dists.size() < 6)
        {
            throw "Error in CalcTetraVolume!  Called with less than six distances";
            return 0;
        }
        else
        {
            return CalcTetraVolume(dists[0], dists[1], dists[2], dists[3], dists[4], dists[5]);
        }
    }

    inline double Nonstringent_acos(double cos)
    {
        if (cos > 1)
        {
            return 0;
        }
        else if (cos < -1)
        {
            return M_PI;
        }
        else
        {
            return acos(cos);
        }
    }

} //namespace chemlib
#endif //MATH_UTIL
