#include <cmath>
#include <cstdio>
#include "mathlib.h"
#include "nonguilib.h"

#include "RibbonSpan.h"

// Cubic Hermite and Cubic B-spline matrices
#define c45 0.70711
#define pi  3.1459265
/*
   static double   cmat[4][4] =   {{-0.5, 1.5, -1.5, 0.5},
          {1.0,  -2.5, 2.0, -0.5},
        {-0.5, 0.0,  0.5, 0.0},
        {0.0,  1.0,  0.0, 0.0}};
 */
#define b0 0.0
#define b1 1.0/6.0
#define b3 3.0/6.0
#define b4 4.0/6.0
#define b6 1.0
static double cmat[4][4] = {{-b1,  b3, -b3, b1},
                            { b3, -b6,  b3, b0},
                            {-b3,  b0,  b3, b0},
                            { b1,  b4,  b1, b0}};

/*
   static double cdata[RB_MP][3];
   static double ndata[RB_MP][3];
   static double cdata[8][3] =  {{0.0, 0.5, 0.0},
        {1.0, 0.25, 0.0},
        {1.5, 0.0, 0.0},
        {1.0, -0.25, 0.0},
        {0.0, -0.5, 0.0},
        {-1.0, -0.25, 0.0},
        {-1.5, 0.0, 0.0},
        {-1.0, 0.25, 0.0}};
   static double ndata[8][3] =   {{0.0, 1.0, 0.0},
          {c45, c45, 0.0},
        {1.0, 0.0, 0.0},
        {c45, -c45, 0.0},
        {0.0, -1.0, 0.0},
        {-c45, -c45, 0.0},
        {-1.0, 0.0, 0.0},
        {-c45, c45, 0.0}};
 */
static double cdata[8][3] =    {{1.5, 0.5, 0.0},
                                {1.5, 0.5, 0.0},
                                {1.5, -0.5, 0.0},
                                {1.5, -0.5, 0.0},
                                {-1.5, -0.5, 0.0},
                                {-1.5, -0.5, 0.0},
                                {-1.5, 0.5, 0.0},
                                {-1.5, 0.5, 0.0}};
static double ndata[8][3] =    {{0.0, 1.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {0.0, -1.0, 0.0},
                                {0.0, -1.0, 0.0},
                                {-1.0, 0.0, 0.0},
                                {-1.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0}};
static double arrowFlat_data[6][3] =   {{2.0, 0.5, 0.0},
                                        {2.0, -0.5, 0.0},
                                        {-2.0, -0.5, 0.0},
                                        {-2.0, 0.5, 0.0},
                                        {0.0, 0.5, 2.0},
                                        {0.0, -0.5, 2.0}};
static double arrowFlat_norm[5][3] =   {{0.0, 0.0, -1.0},
                                        {c45, 0.0, c45},
                                        {-c45, 0.0, c45},
                                        {0.0, 1.0, 0.0},
                                        {0.0, -1.0, 0.0}};

static bool ribbon_flip[4] = {false, false, false, false}; // Control for ribbon flipping to avoid 180 twists



RibbonSpan::RibbonSpan() {
  npr = 8;
  nseg = 10;
  m_pNext = 0;
}

RibbonSpan::~RibbonSpan() {

}

//void RibbonSpan::Draw (int iDrawMode, GLView *pView)
//{
//  pView->DrawPolySurf (npr, nseg, ribbon_pts, ribbon_norms);
//}

bool RibbonSpan::MakeSpan(double ptca[4][3], double ptc[4][3], double pto[4][3], bool bFirst) {
  /* --------------------------------------------------------------------
     Function: Generate one segment of a ribbon from 4 control points.
      The ribbon will either be a cubic hermite or a cubic
      B-spline depending on the setting made by ribbon_type_set().
      The number of profile points, the number of segment points,
      the profile points, and profile point normals are
      set with ribbon_profile_set().
     Input:    ptca -- The CA coordinates for four residues
      ptc -- The C coordinates for four residues
      pto -- the O coordinates for four residues.
      bFirst -- This should be true iff this is the first
        segment.  It is used to initialize direction
        flipping.
     Output:   ribbon_pts -- points returned.  They are returned in the
        order:
        profile pt1
            segment pt1
            segment pt2
            ....
        profile pt2
        ...
      ribbon_norms -- normals returend.
     Notes:    The local direction of the ribbon is from ca->c.  The
      direction of the wide side is c->o.
     ----------------------------------------------------------------------- */

  int i, j, k;
  double mat[4][4], a[3], b[3], c[3], d[3], pnts[RibbonSpan_MP][4][4];
  double norms[RibbonSpan_MP][4][4], t[4];
  double spnts[RibbonSpan_MP][RibbonSpan_MS+1][4], snorms[RibbonSpan_MP][RibbonSpan_MS+1][4];
  double x[4], tval, tinc, bold[3];

  if (bFirst) {
    ribbon_flip[1] = false;
  }

  // Update ribbon_flip
  for (i = 0; i < 3; i++) {
    ribbon_flip[i] = ribbon_flip[i+1];
  }

  // For each point construct a matrix and fill in the pnts and norms
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      a[j] = ptca[i][j] - ptc[i][j];
      b[j] = pto[i][j] - ptc[i][j];
    }
    ml_normalize(a);
    ml_normalize(b);

    // Update the flip direction for the last point or if this is the bFirst
    if ((i == 3) || ((i > 0) && bFirst)) {
      if (ml_dot(b, bold) < 0.0) {
        ribbon_flip[i] = true;
      } else {
        ribbon_flip[i] = false;
      }
    }

    // Flip the direction if needed
    if (ribbon_flip[i]) {
      for (j = 0; j < 3; j++) {
        b[j] = -b[j];
      }
    }

    // Save this direction for latter testing
    for (j = 0; j < 3; j++) {
      bold[j] = b[j];
    }

    // Find a thrird direction for a transform matrix
    ml_cross(a, b, c);
    ml_normalize(c);

    // The b direction is not necessarily orthogonal to a so do another cross
    ml_cross(a, c, b);
    ml_normalize(b);

    // Now construct the transform matrix
    for (j = 0; j < 3; j++) {
      mat[0][j] = b[j];
      mat[1][j] = c[j];
      mat[2][j] = a[j];
      mat[3][j] = 0.0;
      mat[j][3] = 0.0;
    }

    // mat will transform the points so that the x axis lies in the CA-C=O
    //   plane.  Note that normals are not translated.

    for (j = 0; j < npr; j++) {
      ml_vec_mat44(cdata[j], 0.0, mat, d);
      ml_vec_mat44(ndata[j], 0.0, mat, a);
      for (k = 0; k < 3; k++) {
        pnts[j][i][k] = d[k] + ptca[i][k];
        norms[j][i][k] = a[k];
      }
      pnts[j][i][3] = 1.0;
      norms[j][i][3] = 1.0;
    }           // End for each vertex
  }       // End for each of the 4 residues

  // Construct the cubic hermite (cardinal spline) points

  tval = 0.0;
  tinc = 1.0/(double)(nseg-1);
  t[3] = 1.0;
  for (i = 0; i < nseg+1; i++) { // Begin parameteric loop
    t[0] = tval*tval*tval;
    t[1] = tval*tval;
    t[2] = tval;
    ml_vec4_mat44(t, cmat, x);
    for (j = 0; j < npr; j++) {
      ml_vec4_mat44(x, pnts[j], spnts[j][i]);
      ml_vec4_mat44(x, norms[j], snorms[j][i]);
      ml_normalize(snorms[j][i]);
    }
    tval = tval + tinc;
  }         // End parametric loop

  // Construct vertices for the central points
  for (i = 0; i < npr; i++) { // Begin for each parametric step
    for (k = 0; k < nseg; k++) { // Begin for each point in the ribbon
      for (j = 0; j < 3; j++) {
        ribbon_pts[i*nseg+k][j] = spnts[i][k][j];
        ribbon_norms[i*nseg+k][j] = snorms[i][k][j];
      }
    }
  }    // End for each parametric step
  return true;
}

/*
   //////////////////////////////////////////////////////////////////////
   arrowFlat_new (double center[3], double normx[3], double normy[3],
   double arrow_pts[6][3], double arrow_norms[5][3])
   //////////////////////////////////////////////////////////////////////
   //Function: Generate an arrow at center at the orientation given by
   //    the orthogonal normals in orient.
   //Input: center -- Center point for the arrow
   //    normx -- X direction normal
   //    normy -- Y direction normal
   //Output: arrow_pts -- 4 base points and 2 end points
   //    3 ------- 0
   //    |   |
   //    2 ------- 1
   //    arrow_norms -- 5 normals back, right, left, top, bottom

   { // Begin arrowFlat_new()

    double    mat[4][4];
    int     i;

   // Compute a transformation matrix
    for (i=0; i<3; i++)
    {
   mat[0][i] = normx[i];
   mat[1][i] = normy[i];
   mat[3][i] = center[i];
    }
    ml_cross (normx, normy, mat[2]);

   // Transform the points and the normals
    for (i=0; i<6; i++)
   ml_vec_mat44 (arrowFlat_data[i], 1.0, mat, arrow_pts[i]);

    for (i=0; i<5; i++)
   ml_vec_mat44 (arrowFlat_norm[i], 0.0, mat, arrow_norms[i]);

   } // End arrowFlat_new()
 */
//////////////////////////////////////////////////////////////////////
int RibbonSpan::SetProfile(double dx, double dy, int square, int profile_pts,
                           int profile_segs) {
  //////////////////////////////////////////////////////////////////////
  //Function: Generate the profile for a ribbon segment and set
  //    the number of profile points and number of profile
  //    segments.
  //Input:    dx -- half axis size in X
  //    dy -- half axis size in Y
  //    square -- TRUE implies square edges
  //    profile_pts -- Number of points in the profile.
  //      For square this should be 2X since each point
  //      has an assoicated normal
  //    profile_segs -- Number of segments in a profile
  //Value:    FALSE if errors.
  //Effect:   Global parameters set.
  //////////////////////////////////////////////////////////////////////
  // Begin ribbon_segment_new
  double inc, ang, ang_inc;
  double oldx, oldy, oldz;
  int i, j;

  // Set npr and nseg
  if ((profile_pts < 4) || (profile_pts > RibbonSpan_MP)) {
    Logger::log("Ribbon profile points must be between %d and %d.  It is %d\n",
      4, RibbonSpan_MP, profile_pts);
    return false;
  }
  npr = profile_pts;

  if ((profile_segs < 2) || (profile_segs > RibbonSpan_MS)) {
    Logger::log("Ribbon segments must be between %d and %d.  It is %d\n",
      2, RibbonSpan_MS, profile_segs);
    return false;
  }
  nseg = profile_segs;

  // Set the cdata and ndata arrays
  if (square) {
    if (npr != 8) {
      Logger::log("Square option only supports 8 profile points %d\n", npr);
      return false;
    }

    // Initialize everthing to 0 to start with
    for (i = 0; i < npr; i++) {
      for (j = 0; j < 3; j++) {
        cdata[i][j] = 0.0;
        ndata[i][j] = 0.0;
      }
    }

    // Set the data
    cdata[0][0] = cdata[1][0] = dx;
    cdata[0][1] = cdata[1][1] = dy;
    cdata[2][0] = cdata[3][0] = dx;
    cdata[2][1] = cdata[3][1] = -dy;
    cdata[4][0] = cdata[5][0] = -dx;
    cdata[4][1] = cdata[5][1] = -dy;
    cdata[6][0] = cdata[7][0] = -dx;
    cdata[6][1] = cdata[7][1] = dy;

    // Now set the non-zero normals
    ndata[0][1] = 1.0;
    ndata[1][0] = 1.0;
    ndata[2][0] = 1.0;
    ndata[3][1] = -1.0;
    ndata[4][1] = -1.0;
    ndata[5][0] = -1.0;
    ndata[6][0] = -1.0;
    ndata[7][1] = 1.0;

  } else {
    inc = (double)npr;
    ang = 0.0;
    ang_inc = pi*2.0/inc;
    for (i = 0; i < npr; i++) {
      cdata[i][0] = dx*cos(ang);
      cdata[i][1] = dy*sin(ang);
      cdata[i][2] = 0.0;

      ndata[i][0] = cos(ang);
      ndata[i][1] = sin(ang);
      ndata[i][2] = 0.0;
      ang += ang_inc;
    }
  }

  // Update the arrow data
  oldx = arrowFlat_data[0][0];
  oldy = arrowFlat_data[0][1];
  oldz = arrowFlat_data[4][2];
  for (i = 0; i < 6; i++) {
    arrowFlat_data[i][0] = (arrowFlat_data[i][0]/oldx)*dx;
    arrowFlat_data[i][1] = (arrowFlat_data[i][1]/oldy)*dy;
    arrowFlat_data[i][2] = (arrowFlat_data[i][2]/oldz)*dy*4;
  }

  //Reset the flig
  for (i = 0; i < 4; i++) {
    ribbon_flip[i] = false;
  }

  return true;

}

