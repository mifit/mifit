#include <cmath>

/* orthog = returns an orthogonalization matrix for converting
 *          from orthogonal Angstroms to fractions of unit cell;
 * uinv is a general purpose 3x3 matrix inversion routine-
 *          if the matrix from orthog is sent through uinv the
 *          new matrix will convert from fractions to orthogonal A.
 * transform multiplies coordinates by 3x3 matrix
 *          transform(ftoc,&x,&y,&z)
 *
 *  initrotate sets up a matrix for the a 3-d rotation about
 *  the point x,y,z along vector v about angle alpha
 *  call rotate to apply rotation alpha in degrees
 *           initrotate(x,y,z,v1,v2,v3,alpha,matrix[4][3])
 *  apply with
 *	     rotate(x,y,z,&xp,&yp,&zp,mat[4][3])
 *  all floats
 */

void
transform(float ut[3][3], float *x, float *y, float *z)
{
    float xp, yp, zp;
    xp = ut[0][0]*(*x)+ut[0][1]*(*y)+ut[0][2]*(*z);
    yp = ut[1][0]*(*x)+ut[1][1]*(*y)+ut[1][2]*(*z);
    zp = ut[2][0]*(*x)+ut[2][1]*(*y)+ut[2][2]*(*z);
    *x = xp;
    *y = yp;
    *z = zp;
}

int
orthog(float a, float b, float c, float alpha, float beta, float gamma, float ut[3][3])
{
    double degtor;
    double vol, al, be, ga, as, bs, cs, cosals, cosbes, cosgas;
    //double sabgs1;
    degtor = acos(-1.)/180.;
    al = alpha*degtor;
    be = beta*degtor;
    ga = gamma*degtor;
    vol = a* b* c *sqrt(1.-cos(al)*cos(al)-cos(be)*cos(be)
                        -cos(ga)*cos(ga)+2.*(cos(al)*cos(be)*cos(ga)));
    as = (b* c*sin(al))/vol;
    bs = (a* c*sin(be))/vol;
    cs = (a* b*sin(ga))/vol;
    cosals = (cos(be)*cos(ga)-cos(al))/(sin(be)*sin(ga));
    cosbes = (cos(ga)*cos(al)-cos(be))/(sin(ga)*sin(al));
    cosgas = (cos(al)*cos(be)-cos(ga))/(sin(al)*sin(be));
    /*
     *   let u be a 3x3 transformation matrix which relates a new cell
     *   with orthogonal axes of unit dimensions to the original unit
     *   cell (crystal system).
     *    Aug 5, 1991 changed to XPLOR convention was:
     *    ut[0][0]=as;
     *    ut[0][1]=0.;
     *    ut[0][2]=0.;
     *    ut[1][0]=bs*cosgas;
     *    ut[1][1]=1./b;
     *    ut[1][2]= -(1./tan(al))/b;
     *    ut[2][0]=cs*cosbes;
     *    ut[2][1]=0.;
     *    ut[2][2]=1./(c*sin(al));
     *  June 1, 1996 corrections from LFT
     */
    /*sabgs1 = sqrt(1.0-cos(al)*cos(al));*/
    ut[0][0] = 1.F/a;
    ut[0][1] = -(float)((cos(ga))/(sin(ga)*a));
    /*ut[0][2]= -(cos(ga)*sin(be)*cosals + cos(be)*sin(ga))/
       (sin(be)*sabgs1*sin(ga)*a);*/
    ut[0][2] = (float)(as*cosbes);
    ut[1][0] = 0.0F;
    ut[1][1] = 1.F/(float)(sin(ga)*b);
    //ut[1][2]=cosals/(sabgs1*sin(ga)*b);
    ut[1][2] = (float)(bs*cosals);
    ut[2][0] = 0.0F;
    ut[2][1] = 0.0F;
    //ut[2][2]=1.0/(sin(be)*sabgs1*c);
    ut[2][2] = (float)cs;
    return (9);
    /* usage:
       xf = x * ut[0][0] + y * ut[0][1] + z * ut[0][2];
       yf = x * ut[1][0] + y * ut[1][1] + z * ut[1][2];
       zf = x * ut[2][0] + y * ut[2][1] + z * ut[2][2];
     */
}

void uinv(float mat[3][3], float imat[3][3])
{
    float det = mat[0][0]*mat[1][1]*mat[2][2] + mat[1][0]*mat[2][1]*mat[0][2]
          + mat[2][0]*mat[0][1]*mat[1][2] - mat[2][2]*mat[0][1]*mat[1][0]
          - mat[0][0]*mat[2][1]*mat[1][2] - mat[2][0]*mat[1][1]*mat[0][2];
    imat[0][0] = mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2];
    imat[1][0] = mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2];
    imat[2][0] = mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0];
    imat[0][1] = mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2];
    imat[1][1] = mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0];
    imat[2][1] = mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1];
    imat[0][2] = mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1];
    imat[1][2] = mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2];
    imat[2][2] = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            imat[j][i] = imat[j][i]/det;
}

void uinvd(double mat[3][3], double imat[3][3])
{
    double det = mat[0][0]*mat[1][1]*mat[2][2] + mat[1][0]*mat[2][1]*mat[0][2]
          + mat[2][0]*mat[0][1]*mat[1][2] - mat[2][2]*mat[0][1]*mat[1][0]
          - mat[0][0]*mat[2][1]*mat[1][2] - mat[2][0]*mat[1][1]*mat[0][2];
    imat[0][0] = mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2];
    imat[1][0] = mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2];
    imat[2][0] = mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0];
    imat[0][1] = mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2];
    imat[1][1] = mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0];
    imat[2][1] = mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1];
    imat[0][2] = mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1];
    imat[1][2] = mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2];
    imat[2][2] = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            imat[j][i] = imat[j][i]/det;
}

void
initrotate(float a1, float a2, float a3,
           float v1, float v2, float v3,
           float alpha,
           float mat[4][3])
{
    float rho, theta, cal, sal, cph, sph, cth, sth, cph2, sph2, cth2, sth2, pi, cal1;
    float degtor;
    degtor = (float)acos(-1.0)/180.0F;

    alpha *= degtor;
    cal = (float)cos(alpha);
    sal = (float)sin(alpha);
    cal1 = 1.0F-cal;
    rho = (float)sqrt(v1*v1+v2*v2+v3*v3);
    pi = 4.0F * (float)atan(1.0);
    if (rho == 0.0)
    {
        theta = 0;
        cph = 1.0F;
        sph = 0;
    }
    else
    {
        if (v1 == 0.0)
        {
            theta = (v2 >= 0.0 ? 0.5F*pi : 1.5F*pi);
        }
        else
        {
            theta = (float)atan(v2/v1);
            if (v1 < 0)
            {
                theta += pi;
            }
        }
        cph = v3/rho;
        sph = (float)sqrt(1.0F - cph*cph);
    }
    cth = (float)cos(theta);
    sth = (float)sin(theta);
    cph2 = cph*cph;
    sph2 = 1.0F - cph2;
    cth2 = cth*cth;
    sth2 = 1.0F - cth2;

    mat[0][0] = (cal*cph2+sph2)*cth2+cal*sth2;
    mat[0][1] = sal*cph+cal1*sph2*cth*sth;
    mat[0][2] = sph*(cph*cth*cal1-sal*sth);
    mat[1][0] = sph2*cth*sth*cal1-sal*cph;
    mat[1][1] = sth2*(cal*cph2+sph2)+cal*cth2;
    mat[1][2] = sph*(cph*sth*cal1+sal*cth);
    mat[2][0] = sph*(cph*cth*cal1+sal*sth);
    mat[2][1] = sph*(cph*sth*cal1-sal*cth);
    mat[2][2] = cal*sph2+cph2;
    mat[3][0] = a1-a1*mat[0][0]-a2*mat[1][0]-a3*mat[2][0];
    mat[3][1] = a2-a1*mat[0][1]-a2*mat[1][1]-a3*mat[2][1];
    mat[3][2] = a3-a1*mat[0][2]-a2*mat[1][2]-a3*mat[2][2];
    /* whew! */
}

void rotate(float *xp, float *yp, float *zp, float mat[4][3])
{
    float x = *xp;
    float y = *yp;
    float z = *zp;
    *xp = x*mat[0][0] + y*mat[1][0] + z*mat[2][0] + mat[3][0];
    *yp = x*mat[0][1] + y*mat[1][1] + z*mat[2][1] + mat[3][1];
    *zp = x*mat[0][2] + y*mat[1][2] + z*mat[2][2] + mat[3][2];
}

void rotate(float v[3], float mat[4][3])
{
    float x = v[0];
    float y = v[1];
    float z = v[2];
    v[0] = x*mat[0][0] + y*mat[1][0] + z*mat[2][0] + mat[3][0];
    v[1] = x*mat[0][1] + y*mat[1][1] + z*mat[2][1] + mat[3][1];
    v[2] = x*mat[0][2] + y*mat[1][2] + z*mat[2][2] + mat[3][2];
}

#define dtor(x) ( (x) * (3.141592653589793F/180.0F))

void buildmatz(float z, float mat[3][3])
{
    float sz = (float)sin(dtor(z));
    float cz = (float)cos(dtor(z));
    mat[0][0] = cz;
    mat[0][1] = sz;
    mat[0][2] = 0;

    mat[1][0] = -sz;
    mat[1][1] = cz;
    mat[1][2] = 0;

    mat[2][0] = 0;
    mat[2][1] = 0;
    mat[2][2] = 1;
}

void buildmat(float xdegrees, float ydegrees, float zdegrees, float mat[3][3])
{
    /* small-angle approximation matrix */
    float sx, cx, sy, cy, sz, cz;
    if (xdegrees == 0)
    {
        sx = 0;
        cx = 1;
    }
    else
    {
        sx = (float)sin(dtor(xdegrees));
        cx = (float)cos(dtor(xdegrees));
    }
    if (ydegrees == 0)
    {
        sy = 0;
        cy = 1;
    }
    else
    {
        sy = (float)sin(dtor(ydegrees));
        cy = (float)cos(dtor(ydegrees));
    }
    if (zdegrees == 0)
    {
        sz = 0;
        cz = 1;
    }
    else
    {
        sz = (float)sin(dtor(zdegrees));
        cz = (float)cos(dtor(zdegrees));
    }

    mat[0][0] = cy*cz;
    mat[0][1] = sz;
    mat[0][2] = sy;

    mat[1][0] = -sz;
    mat[1][1] = cx*cz;
    mat[1][2] = -sx;

    mat[2][0] = -sy;
    mat[2][1] = sx*cy;
    mat[2][2] = cx*cy;
}

void xl_rotate(float x, float y, float z, float *xp, float *yp, float *zp, float mat[4][3])
{
    *xp = x*mat[0][0] + y*mat[1][0] + z*mat[2][0] + mat[3][0];
    *yp = x*mat[0][1] + y*mat[1][1] + z*mat[2][1] + mat[3][1];
    *zp = x*mat[0][2] + y*mat[1][2] + z*mat[2][2] + mat[3][2];
}

