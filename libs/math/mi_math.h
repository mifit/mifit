#ifndef mifit_math_math_h
#define mifit_math_math_h

#include <limits>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/**
 * Compute the distance between 2 points
 * Input:    a, b -- Points
 * Value:    Distance between the points
 */
double ml_ptpt_distance(double a[3], double b[3]);

/*
 * Determine if a point is inside a triangle.
 * Input:    pt -- Point to be tested.
 *     pta, ptb, ptc -- Triangle vertices.
 * Output:   pt_bc -- Intersection of the ray from vertex A through
 *     the point with the side BC.  This is useful for
 *     interpolating normals, colors etc. at pt.
 * Value:    TRUE if inside else FALSE
 */
int ml_pttri_inside(double pt[3], double pta[3], double ptb[3], double ptc[3], double pt_bc[3]);

/*
 * Normalize a vector (i.e. set to unit length)
 * Input:    a  vector to be normalized
 * Output:   a  normalized vector
 *     NOTE: Overwrites old values
 * Value:    FALSE if vector has almost 0 length.
 */
int ml_normalize(double a[3]);

/*
 * Calculate the cross product c = a x b
 * Input:    a,b  vectors to cross
 * Output:   c  cross product
 */
void ml_cross(double a[3], double b[3], double c[3]);

/*
 * Calculate the dot product  a (dot) b
 * Input:    a,b  vectors to dot
 * Value:    Dot product.
 */
double ml_dot(double a[3], double b[3]);

/*
 * Multiply a vector (1x3) times a matrix (3x3)
 * Input   a -- vector
 *       mat --  matrix
 * Output:   b -- transformed vector
 */
void ml_vec_mat33(double a[3], double mat[3][3], double b[3]);

/*
 * Multiply a vector (1x3) times a matrix (4x4).  The fourth
 *     element in the vector is the homogeneous element w.
 * Input   a -- vector
 *     vw -- homogeneous coordinate for the vector.
 *       mat --  matrix
 * Output:   b -- transformed vector
 */
void ml_vec_mat44(double a[3], double vw, double mat[4][4], double b[3]);

/*
 * Multiply a vector (1x4) times a matrix (4x4)
 * Input   a -- vector
 *       mat --  matrix
 * Output:   b -- transformed vector (must be [4])
 */
void ml_vec4_mat44(double a[4], double mat[4][4], double b[4]);

/*
 * Compute a plane definition from 3 points.
 * Input:    a,b,c -- Point arrays
 * Output:   plane -- normalized normal and "d"
 * Value:    FALSE if points are coincident or colinear.
 */
int ml_plane_define(double a[3], double b[3], double c[3], double plane[4]);

/*
 * Test point a to see if it lies on (close to) the plane.
 * Input:    plane -- plane equation array
 *     a -- point array
 * Value:    FALSE if point too far from the plane definition
 *     TRUE if point close enough to the plane definition.
 */
int ml_plane_test(double plane[4], double a[3]);

/*
 * Compute a line equation.
 *     (x-x0)/(x1-x0) = (y-y0)/(y1-y0) = (z-z0)/(z1-z0)
 *     0-2 = x0, y0; z0;
 *     3-5 = (x1-x0), etc. -- Normalized.
 * Input:    pta, ptb -- Endpoints for line segement.
 * Output:   lineq -- Line equation.
 * Value:    False if points are coincident.
 */
int ml_line_define(double pta[3], double ptb[3], double lineq[6]);

/*
 * Intersect a line and a plane.
 * Input:    line -- line equation.
 *     plane -- plane equation.
 * Output:   pt -- Intersection point.
 * Value:    False if parallel.
 */
int ml_lineplane_intersect(double lineq[6], double plane[4], double pt[3]);

/*
 * Intersect a line and a sphere.
 * Input:    lineq -- line equation.
 *     center -- sphere center
 *     radius -- sphere radius
 * Output:   pta -- first intersection point (min z).
 *     ptb -- second intersection point (max z).
 * Value:    False if no intersection.
 */
int ml_linesphere_intersect(double lineq[6], double center[3], double radius, double pta[3], double ptb[3]);

/*
 * Intersect two line segments defined by four points.
 * Input:    pta1, pta2 -- Endpoints for line segement a.
 *     ptb1, ptb2 -- Endpoints for line segement b.
 * Output:   pt_intersect -- Intersection point (undefined if
 *     FALSE)
 * Value:    FALSE if segments are parallel or intersection does not
 *     lie on both line segments
 */
int ml_lseglseg_intersect(double pta[3], double ptb[3], double ptc[3], double ptd[3], double np[3], double pt_intersect[3], int debug);

/*
 * Multiply two 4x4 matrices together [matc] = [mata][matb]
 * Input:    mata, matb -- matrices to be multiplied
 * Output:   matc -- resultant matrix
 * Note:   The result matrix is first accumlated in a local data structure
 *     so matc may be the same as mata or matb.  (e.g.
 *     ml_matmat (mata, matb, matb) works ok.
 */
void ml_matmat44(double mata[4][4], double matb[4][4], double matc[4][4]);

/*
 * Compute a transform matrix such that
 *     [normal][mat] = [0,0,1].
 * Input:    normal -- normal to be transformed (re-normalized)
 * Output:   mat -- tranform matrix
 */
void ml_mat_normalTransformZ(double normal[3], double mat[3][3]);

/*
 * Transpose a matrix (Note that this may be used for
 *     inverting orthonormal rotation matrices).
 * Input:    mat -- matrix to be transposed
 * Output:   mat_trans -- transposed matrix
 * Note:   The matrix is first transposed to an internal matrix so
 *     ml_mat_transpose (mat, mat) works.
 */
void ml_mat_transpose33(double mat[3][3], double mat_trans[3][3]);


//@{
// Translates a float to a character string.
//@}
const char* ftoa(float);
//@{
// Translates a long integer to a character string.
//@}
const char* ftoa(long);
//@{
// Translates a integer to a character string.
//@}
const char* ftoa(int);

template<typename T>
inline bool equals(T a, T b) {
  return abs(a - b) <= std::numeric_limits<T>::epsilon();
}

inline double toRadians(double degrees) {
  return degrees*M_PI/180.0;
}

inline double toDegrees(double radians) {
  return radians*180.0/M_PI;
}

#if defined ( _MSC_VER ) && !defined ( __NUTC__ )
inline double rint(double v) {
  double theFloor = floor(v);
  double theCeil  = ceil(v);

  if (fabs(v - theFloor) < fabs(v - theCeil) ) {
    return theFloor;
  }

  return theCeil;
}
#endif

#endif
