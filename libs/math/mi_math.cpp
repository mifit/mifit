#include  <stdio.h>
#include  "mi_math.h"
#include  <cmath>
#include <utillib.h>

double ml_plane_maxdiff = 0.01;
double ml_vector_minlen = 0.0001;
double ml_ptpt_mindist = 0.0001;
double ml_parallel_min = 0.001;

double ml_ptpt_distance(double a[3], double b[3]) {

  double sum, diff;
  int i;

  sum = 0.0;
  for (i = 0; i < 3; i++) {
    diff = a[i]-b[i];
    sum = sum + diff*diff;
  }

  if (sum > ml_ptpt_mindist) {
    return sqrt(sum);
  } else {
    return 0.0;
  }

}

int ml_pttri_inside(double pt[3], double pta[3], double ptb[3], double ptc[3], double pt_bc[3]) {

  int i, status;
  double dist, distAB, distAC, dist_pt, dist_bc;
  double pt_ray[3];
  static double normal[3] = {0.0, 0.0, 1.0};

  // Method
  //   Construct a ray from A through the point pt and test for intersection
  //   with the edge BC.  Actually we will use a line segment from A through
  //   the point which is a little longer than the the max distance from
  //   A to B or C as our "ray".
  //

  // Find the max distance for the "ray"
  distAB = ml_ptpt_distance(pta, ptb);
  distAC = ml_ptpt_distance(pta, ptc);
  if (distAB > distAC) {
    dist = distAB * 1.1;
  } else {
    dist = distAC * 1.1;
  }

  // If the point is further than this, it cannot be inside
  dist_pt = ml_ptpt_distance(pta, pt);
  if (dist_pt > dist) {
    return false;
  }

  // Compute a ray point by scaling
  for (i = 0; i < 3; i++) {
    pt_ray[i] = (pt[i] - pta[i])*(dist/dist_pt) + pta[i];
  }

  // Intersect ray with BC
  status = ml_lseglseg_intersect(pta, pt_ray, ptb, ptc, normal, pt_bc, false);
  if (!status) {
    return false;
  }
  dist_bc = ml_ptpt_distance(pta, pt_bc);
  if (dist_bc >= dist_pt) {
    return true;
  } else {
    return false;
  }

}

int ml_normalize(double a[3]) {

  double len;
  int i;

  // Find length and divide each term
  len = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  if (len < ml_vector_minlen) {
    return false;
  }

  for (i = 0; i < 3; i++) {
    a[i] = a[i]/len;
  }

  return true;
}

void ml_cross(double a[3], double b[3], double c[3]) {

  int i;
  double temp[3];

  temp[0] = a[1]*b[2] - a[2]*b[1];
  temp[1] = a[2]*b[0] - a[0]*b[2];
  temp[2] = a[0]*b[1] - a[1]*b[0];

  for (i = 0; i < 3; i++) {
    c[i] = temp[i];
  }

}

double ml_dot(double a[3], double b[3]) {

  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);

}

void ml_vec_mat33(double a[3], double mat[3][3], double b[3]) {

  int i;

  for (i = 0; i < 3; i++) {
    b[i] = a[0]*mat[0][i] + a[1]*mat[1][i] + a[2]*mat[2][i];
  }

}

void ml_vec_mat44(double a[3], double vw, double mat[4][4], double b[3]) {

  int i;

  for (i = 0; i < 3; i++) {
    b[i] = a[0]*mat[0][i] + a[1]*mat[1][i] + a[2]*mat[2][i] + vw*mat[3][i];
  }

}

void ml_vec4_mat44(double a[4], double mat[4][4], double b[4]) {

  int i;

  for (i = 0; i < 4; i++) {
    b[i] = a[0]*mat[0][i] + a[1]*mat[1][i] + a[2]*mat[2][i] +a[3]*mat[3][i];
  }

}

int ml_plane_define(double a[3], double b[3], double c[3], double plane[4]) {

  double ab[3], ac[3];
  int i;

  // Construct the difference vectors
  for (i = 0; i < 3; i++) {
    ab[i] = b[i] - a[i];
    ac[i] = c[i] - a[i];
  }

  // Find cross product and normalize
  ml_cross(ab, ac, plane);
  if (!ml_normalize(plane)) {
    return false;
  }

  // Substitue point a back in to find d
  plane[3] = -(a[0]*plane[0]+a[1]*plane[1]+a[2]*plane[2]);

  return true;
}

int ml_plane_test(double plane[4], double a[3]) {

  double diff;

  diff = a[0]*plane[0] + a[1]*plane[1] + a[2]*plane[2] + plane[3];
  if (fabs(diff) > ml_plane_maxdiff) {
    return false;
  } else {
    return true;
  }

}

int ml_line_define(double pta[3], double ptb[3], double lineq[6]) {

  int i;
  for (i = 0; i < 3; i++) {
    lineq[i] = pta[i];
    lineq[i+3] = ptb[i]-pta[i];
  }

  // Normalize dx, dy, dz
  return ml_normalize(&lineq[3]);

}

int ml_lineplane_intersect(double lineq[6], double plane[4], double pt[3]) {

  int ix, iy, iz, dx, dy, dz;
  double divisor, yterm, zterm, dydx, dzdx;
  double dmax;

  // Method
  //   The plane equation is:
  //  ax + by + cz + d = 0        (1)
  //   The line equation is
  //  x-x0/dx = y-y0/dy = z-z0/dz     (2)
  //   Solving in terms of x
  //  y = (x-x0)*(dy/dx) + y0       (3)
  //  z = (x-x0)*(dz/dx) + z0       (3)
  //
  //  y = x*(dy/dx) - x0*dy/dx + y0     (4)
  //  z = x*(dz/dx) - x0*dz/dx + z0     (4)
  //    Substituting in equation 1
  //  ax + b(x*dy/dx - x0*dy/dx + y0) + c(x*dz/dx - X0*dz/dx + z0) + d = 0;
  //    Collecting terms
  //  x*(a + b*dy/dx + c*dz/dx) = b*x0*dy/dx - b*y0 +
  //            c*x0*dz/dx - c*z0 - d
  //    Introducing the variables used in the code
  //  x*divisor = yterm + zterm - d
  //  x = (yterm + zterm - d)/divisor
  //    Substitute into equation 3 for y and z.
  //

  // Use the largest direction component as "x"
  dmax = lineq[3];
  ix = 0;
  if (fabs(dmax) < fabs(lineq[4])) {
    dmax = lineq[4];
    ix = 1;
  }
  if (fabs(dmax) < fabs(lineq[5])) {
    ix = 2;
  }

  if (ix == 0) {
    iy = 1; iz = 2;
    dx = 3; dy = 4; dz = 5;
  } else if (ix == 1) {
    iy = 2; iz = 0;
    dx = 4; dy = 5; dz = 3;
  } else { // ix == 2
    iy = 0; iz = 1;
    dx = 5; dy = 3; dz = 4;
  }

  dydx = lineq[iy]/lineq[ix];
  dzdx = lineq[iz]/lineq[ix];

  // Compute divisor and test for 0
  divisor = plane[ix] + plane[iy]*dydx + plane[iz]*dzdx;
  if (fabs(divisor) < ml_parallel_min) {
    return false;
  }

  // Compute terms
  yterm = plane[iy]*lineq[ix]*dydx - plane[iy]*lineq[iy];
  zterm = plane[iz]*lineq[ix]*dzdx - plane[iz]*lineq[iz];

  // Compute x and substitue back into the line equations (3)
  pt[ix] = (yterm + zterm - plane[3])/divisor;
  pt[iy] = (pt[ix]-lineq[ix])*dydx + lineq[iy];
  pt[iz] = (pt[ix]-lineq[ix])*dzdx + lineq[iz];

  return true;

}

int ml_linesphere_intersect(double lineq[6], double center[3],
                            double radius, double pta[3], double ptb[3]) {

  int i;
  double z2;
  double pt_off[3], pt_trans[3];
  double mat[3][3], mat_inv[3][3];

  // Transform the line so that it lines up with the Z axis to make this simple
  ml_mat_normalTransformZ(&lineq[3], mat);

  // Transform the coordinates of the line
  for (i = 0; i < 3; i++) {
    pt_off[i] = lineq[i] - center[i];
  }
  ml_vec_mat33(pt_off, mat, pt_trans);

  // Compute z**2 and return with a miss if it is less than 0
  z2 = radius*radius - pt_trans[0]*pt_trans[0] - pt_trans[1]*pt_trans[1];
  if (z2 < ml_parallel_min) {
    return false;
  }

  pt_trans[2] = sqrt(z2);

  // Transform back and offset
  ml_mat_transpose33(mat, mat_inv);

  ml_vec_mat33(pt_trans, mat_inv, pta);

  pt_trans[2] = -pt_trans[2];
  ml_vec_mat33(pt_trans, mat_inv, ptb);

  // Add the center back in
  for (i = 0; i < 3; i++) {
    pta[i] = pta[i] + center[i];
    ptb[i] = ptb[i] + center[i];
  }

  // Make sure pta is in front
  if (pta[2] > ptb[2]) {
    for (i = 0; i < 3; i++) {
      pt_off[i] = pta[i];
      pta[i] = ptb[i];
      ptb[i] = pt_off[i];
    }
  }
  return true;

}

double ml_lsegpt_param(double pta[3], double ptb[3], double ptc[3]) {

  int i, np;
  double param;
  double lineq[6];

  // Define a line to give us x0, y0, z0 and dx, dy, dz
  ml_line_define(pta, ptb, lineq);

  // Average the parametric value for cases where dx, etc is not 0
  param = 0.0;
  np = 0;
  for (i = 0; i < 3; i++) {
    if (fabs(lineq[3+i]) > ml_parallel_min) {
      param = param+(ptc[i]-lineq[i])/lineq[i+3];
      np++;
    }
  }
  if (np) {
    return param/(double)np;
  } else {
    return param;
  }

}

int ml_lseglseg_intersect(double pta[3], double ptb[3], double ptc[3],
                          double ptd[3], double*,
                          double pt_intersect[3], int debug) {
  double ad[3], ab[3], ac[3], cd[3], ca[3], cb[3];
  double np1[3], np2[3], plane[3];
  double pt_intersectab[3], pt_intersectcd[3];
  double d, parameter_ab, parameter_cd;
  double dist1, dist2;
  int i;
  double np[3];

  // Compute diffence vertors
  for (i = 0; i < 3; i++) {
    ad[i] = ptd[i] - pta[i];
    ab[i] = ptb[i] - pta[i];
    ac[i] = ptc[i] - pta[i];
    cd[i] = ptd[i] - ptc[i];
    ca[i] = -ac[i];
    cb[i] = ptb[i] - ptc[i];
  }

  // Find normals to planes of triples Note: order of cross is vital here*/
  ml_cross(ad, ab, np1);
  ml_cross(ab, ac, np2);
  if (debug) {
    printf("np1 %.4f %.4f %.4f \n", np1[0], np1[1], np1[2]);
    printf("np2 %.4f %.4f %.4f \n", np2[0], np2[1], np2[2]);
  }

  // If cross product is 0 then ABD are colinear which means no intersect
  //   for the purposes of this function
  if (!ml_normalize(np1)) {
    return false;
  }
  if (!ml_normalize(np2)) {
    return false;
  }

  // If normals point the opposite directions, then C and D are on the same
  //   side of AB and no intersection is possible  Otherwise average them
  if (ml_dot(np1, np2) < 0.0) {
    if (debug) {
      printf("Opposite normals\n");
      printf("np1 %.4f %.4f %.4f\n", np1[0], np1[1], np1[2]);
      printf("np2 %.4f %.4f %.4f\n", np2[0], np2[1], np2[2]);
    }
    return false;
  }

  ml_cross(ad, ac, np);
  if (!ml_normalize(np)) {
    return false;
  }

  // Now cross np with ab to construct a plane which contains AB
  ml_cross(ab, np, plane);
  ml_normalize(plane);
  d = -(plane[0]*pta[0] + plane[1]*pta[1] + plane[2]*pta[2]);

  // Compute the distance to the plane from C and D and use the ratio to
  // determine the total distances
  dist1 = ml_dot(plane, ptc) + d;
  dist2 = ml_dot(plane, ptd) + d;
  if (debug) {
    printf("cd dist1 %.4f dist2 %.4f\n", dist1, dist2);
  }
  if ((dist1 < 0.0) && (dist2 < 0.0)) {
    return false;
  }
  if ((dist1 > 0.0) && (dist2 > 0.0)) {
    return false;
  }
  dist1 = fabs(dist1);
  dist2 = fabs(dist2);
  parameter_cd  = dist1/(dist1+dist2);

  // Do the same process for the other line segement
  ml_cross(cb, ca, np);
  if (!ml_normalize(np)) {
    return false;
  }

  // Now cross np with cd to construct a plane which contains CD
  ml_cross(cd, np, plane);
  ml_normalize(plane);
  d = -(plane[0]*ptc[0] + plane[1]*ptc[1] + plane[2]*ptc[2]);

  // Compute the distance to the plane from A and B and use the ratio to
  //   determine the total distances
  dist1 = ml_dot(plane, pta) + d;
  dist2 = ml_dot(plane, ptb) + d;
  if (debug) {
    printf("ab dist1 %.4f dist2 %.4f\n", dist1, dist2);
  }
  if ((dist1 < 0.0) && (dist2 < 0.0)) {
    return false;
  }
  if ((dist1 > 0.0) && (dist2 > 0.0)) {
    return false;
  }
  dist1 = fabs(dist1);
  dist2 = fabs(dist2);
  parameter_ab  = dist1/(dist1+dist2);

  // Compute the point based on the parametric value for CD and compute a
  //   prameteric value for AB to test for inclusion
  for (i = 0; i < 3; i++) {
    pt_intersectcd[i] = parameter_cd*cd[i] + ptc[i];
    pt_intersectab[i] = parameter_ab*ab[i] + pta[i];
    pt_intersect[i] = (pt_intersectcd[i]+pt_intersectab[i])/2.0;
  }
  if (debug) {
    //
    //  point_draw (pt_intersectab);
    //  point_draw (pt_intersectcd);
    //  point_draw (pt_intersect);
    //  ptn_draw (pta, plane, 0.2);
    //
    printf("dist1 %.4f dist2 %.4f parameter_cd %.4f parameter_ab %.4f\n",
      dist1, dist2, parameter_cd, parameter_ab);
    printf("lseg: Enter a number to continue >");
    scanf("%d", &i);
  }

  // Make sure intersection is in both line segements -- endpoint intersections
  //   don't count
  if ((parameter_ab < 0.01) || (parameter_ab > 0.99)) {
    return false;
  }
  if ((parameter_cd < 0.01) || (parameter_cd > 0.99)) {
    return false;
  }

  return true;

}

void ml_matmat44(double mata[4][4], double matb[4][4], double matc[4][4]) {

  int i, j;
  double matr[4][4];

  // For each row of matrix a multiply by colum of matrix b
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      matr[i][j] = mata[i][0] * matb[0][j] + mata[i][1] * matb[1][j] +
                   mata[i][2] * matb[2][j] + mata[i][3] * matb[3][j];
    }
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      matc[i][j] = matr[i][j];
    }
  }

}

void ml_mat_normalTransformZ(double normal[3], double mat[3][3]) {
  int i, ix = 0, iy, iz;
  double max;
  double mat_inv[3][3];

  // Method
  //   Besides [normal][mat] = [0, 0, 1] we need two more equations.
  //   We choose these to give us the matrix equation:
  //  [ ][  mat ] = [1, 0, 0]
  //  [ ][  ] = [0, 1, 0]
  //  [normal ][  ] = [0, 0, 1]
  //   Thus we know that mat is simply the inverse (transpose) of the
  //   matrix with "normal" as its third line.  The trick is to find
  //   the first two lines.  The only requirement is that all three
  //   lines are orthonormal, since we don't care how the transform
  //   "spins" about the z axis.  So we pick some unit vector normal
  //   to line 3 as line 2 and then cross lines 2 and 3 to create line 1.
  //

  // Third row of the matrix is the normal
  for (i = 0; i < 3; i++) {
    mat_inv[2][i] = normal[i];
  }
  if (!ml_normalize(mat_inv[2])) {
    printf("ml_mat_normalTransformZ: 0 normal\n");
    return;
  }

  // Second row is some normal vector X such that X dot normal = 0
  max = 0.0;
  for (i = 0; i < 3; i++) {
    if (fabs(max) < fabs(normal[i])) {
      max = normal[i];
      ix = i;
    }
  }

  iy = ix + 1;
  if (iy > 2) {
    iy = 0;
  }

  iz = iy + 1;
  if (iz > 2) {
    iz = 0;
  }

  // ab - ba + 0 == 0
  mat_inv[1][ix] = mat_inv[2][iy];
  mat_inv[1][iy] = -mat_inv[2][ix];
  mat_inv[1][iz] = 0.0;
  ml_normalize(mat_inv[1]);

  // Now cross to get the first row
  ml_cross(mat_inv[1], mat_inv[2], mat_inv[0]);

  // Transpose
  ml_mat_transpose33(mat_inv, mat);

}

void ml_mat_transpose33(double mat[3][3], double mat_trans[3][3]) {

  int i, j;
  double mat_temp[3][3];

  // Transpose
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      mat_temp[i][j] = mat[j][i];
    }
  }

  // Copy to output
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      mat_trans[i][j] = mat_temp[i][j];
    }
  }

}

const char* ftoa(int f) {
  static std::string s;
  s = format("%d", f);
  return s.c_str();
}

const char* ftoa(long f) {
  static std::string s;
  s = format("%ld", f);
  return s.c_str();
}

const char* ftoa(float f) {
  static std::string s;
  s = format("%0.3f", f);
  return s.c_str();
}

