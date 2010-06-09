#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <chemlib/Residue.h>

#include "Helix.h"
#include "RESIDUE.h"

#include <vector>

using namespace chemlib;

Helix::Helix(double radius)
    : m_pNext(NULL),
      m_dRadius(radius)
{
}

Helix::~Helix()
{
}

bool Helix::MakeHelix(MIMoleculeBase* /* mol */, Residue *pHelixStart, Residue *pHelixStop)
{
    // Find the HELIX subset of this molecule
    if ((pHelixStart == 0) || (pHelixStop == 0))
    {
        return false;
    }

    //Count the number of residues
    std::vector<Residue*> pHelixFrags;
    MIIter<Residue> pHelixNext = Residue::getIter(pHelixStart);
    for (; Monomer::isValid(pHelixNext) && pHelixNext != pHelixStop; ++pHelixNext)
    {
        pHelixFrags.push_back(pHelixNext);
    }
    //Add the last one and one more because sequences have been shorted by 1 [old comment]
    if (Monomer::isValid(pHelixNext))
    {
        pHelixFrags.push_back(pHelixNext);
        ++pHelixNext;
        if (Monomer::isValid(pHelixNext))
        {
            pHelixFrags.push_back(pHelixNext);
        }
    }

    //There must be at least 4 residues
    int iAxisLen = pHelixFrags.size() -1;
    if (pHelixFrags.size() < 4)
    {
        return false;
    }

    //Get the values of the CA atoms for each residues
    double caA[4][3], caB[4][3];
    for (int i = 0; i < 4; i++)
    {
        MIAtom *pAtomA = atom_from_name("CA", *pHelixFrags[i]);
        MIAtom *pAtomB = atom_from_name("CA", *pHelixFrags[iAxisLen-i]);
        if (!pAtomA || !pAtomB)
        {
            return false;
        }
        caA[i][0] = pAtomA->x();
        caA[i][1] = pAtomA->y();
        caA[i][2] = pAtomA->z();
        caB[i][0] = pAtomB->x();
        caB[i][1] = pAtomB->y();
        caB[i][2] = pAtomB->z();
    }

    //Create the axis (Stored in member variables)
    CreateAxis(iAxisLen, caA, caB);

    //Create the cylinder data
    CreateCylinder();

    // Adjust the axis a bit
    for (unsigned int k = 0; k < 3; k++)
    {
        m_dCapNormA[0][k] = m_dAxisB[k]-m_dAxisA[k];
    }
    ml_normalize(m_dCapNormA[0]);

    for (unsigned int k = 0; k < 3; k++)
    {
        m_dAxisA[k] += m_dCapNormA[0][k];
        m_dAxisB[k] += m_dCapNormA[0][k];
    }

    double capNormPA[3], capNormPB[3], capNormPC[3];

    //For each segemnt create the triangle for the caps.
    for (int j = 0; j < Helix_NSEG; j++)
    {
        // Compute the normals for the capA
        for (int k = 0; k < 3; k++)
        {
            capNormPA[k] = m_dEndA[j][k] - m_dAxisA[k];
            if (j == Helix_NSEG-1)
            {
                capNormPB[k] = m_dEndA[0][k] - m_dAxisA[k];
                capNormPC[k] = m_dEndA[1][k] - m_dAxisA[k];
            }
            else if (j == Helix_NSEG-2)
            {
                capNormPB[k] = m_dEndA[j+1][k] - m_dAxisA[k];
                capNormPC[k] = m_dEndA[0][k] - m_dAxisA[k];
            }
            else
            {
                capNormPB[k] = m_dEndA[j+1][k] - m_dAxisA[k];
                capNormPC[k] = m_dEndA[j+2][k] - m_dAxisA[k];
            }
        }
        ml_normalize(capNormPA);
        ml_normalize(capNormPB);
        ml_normalize(capNormPC);
        ml_cross(capNormPA, capNormPB, m_dCapNormA[j]);
        ml_cross(capNormPB, capNormPC, m_dCapNormB[j]);

        for (int k = 0; k < 3; k++)
        {
            m_dCapNormA[j][k] = -(m_dCapNormA[j][k]+m_dCapNormB[j][k])/2.0;
        }

        // And now for cap B
        for (unsigned int k = 0; k < 3; k++)
        {
            capNormPA[k] = m_dEndB[j][k] - m_dAxisB[k];
            if (j == Helix_NSEG-1)
            {
                capNormPB[k] = m_dEndB[0][k] - m_dAxisB[k];
                capNormPC[k] = m_dEndB[1][k] - m_dAxisB[k];
            }
            else if (j == Helix_NSEG-2)
            {
                capNormPB[k] = m_dEndB[j+1][k] - m_dAxisB[k];
                capNormPC[k] = m_dEndB[0][k] - m_dAxisB[k];
            }
            else
            {
                capNormPB[k] = m_dEndB[j+1][k] - m_dAxisB[k];
                capNormPC[k] = m_dEndB[j+2][k] - m_dAxisB[k];
            }
        }
        ml_normalize(capNormPA);
        ml_normalize(capNormPB);
        ml_normalize(capNormPC);
        ml_cross(capNormPB, capNormPA, m_dCapNormB[j]);
        ml_cross(capNormPC, capNormPB, capNormPA);

        for (unsigned int k = 0; k < 3; k++)
        {
            m_dCapNormB[j][k] = -(capNormPA[k]+m_dCapNormB[j][k])/2.0;
        }
    }
    /*
          vertices[n_vertex].coord.x = m_dAxisA[i][0];
          vertices[n_vertex].coord.y = m_dAxisA[i][1];
          vertices[n_vertex].coord.z = m_dAxisA[i][2];
          vertices[n_vertex].normal.i = m_dCapNormA[0];
          vertices[n_vertex].normal.j = m_dCapNormA[1];
          vertices[n_vertex].normal.k = m_dCapNormA[2];
          vertices[n_vertex].color.red = red;
          vertices[n_vertex].color.green = green;
          vertices[n_vertex].color.blue = blue;
          vertices[n_vertex].color.alpha = alpha;
          n_vertex++;

          vertices[n_vertex].coord.x = m_dEndA[i][j][0];
          vertices[n_vertex].coord.y = m_dEndA[i][j][1];
          vertices[n_vertex].coord.z = m_dEndA[i][j][2];
          vertices[n_vertex].normal.i = m_dCapNormA[0];
          vertices[n_vertex].normal.j = m_dCapNormA[1];
          vertices[n_vertex].normal.k = m_dCapNormA[2];
          vertices[n_vertex].color.red = red;
          vertices[n_vertex].color.green = green;
          vertices[n_vertex].color.blue = blue;
          vertices[n_vertex].color.alpha = alpha;
          n_vertex++;

          vertices[n_vertex].coord.x = m_dEndA[i][j][0];
          vertices[n_vertex].coord.y = m_dEndA[i][j][1];
          vertices[n_vertex].coord.z = m_dEndA[i][j][2];
          vertices[n_vertex].normal.i = m_dNormA[i][j][0];
          vertices[n_vertex].normal.j = m_dNormA[i][j][1];
          vertices[n_vertex].normal.k = m_dNormA[i][j][2];
          vertices[n_vertex].color.red = red;
          vertices[n_vertex].color.green = green;
          vertices[n_vertex].color.blue = blue;
          vertices[n_vertex].color.alpha = alpha;
          n_vertex++;

          vertices[n_vertex].coord.x = m_dEndB[i][j][0];
          vertices[n_vertex].coord.y = m_dEndB[i][j][1];
          vertices[n_vertex].coord.z = m_dEndB[i][j][2];
          vertices[n_vertex].normal.i = m_dNormB[i][j][0];
          vertices[n_vertex].normal.j = m_dNormB[i][j][1];
          vertices[n_vertex].normal.k = m_dNormB[i][j][2];
          vertices[n_vertex].color.red = red;
          vertices[n_vertex].color.green = green;
          vertices[n_vertex].color.blue = blue;
          vertices[n_vertex].color.alpha = alpha;
          n_vertex++;

          vertices[n_vertex].coord.x = m_dEndB[i][j][0];
          vertices[n_vertex].coord.y = m_dEndB[i][j][1];
          vertices[n_vertex].coord.z = m_dEndB[i][j][2];
          vertices[n_vertex].normal.i = m_dCapNormB[0];
          vertices[n_vertex].normal.j = m_dCapNormB[1];
          vertices[n_vertex].normal.k = m_dCapNormB[2];
          vertices[n_vertex].color.red = red;
          vertices[n_vertex].color.green = green;
          vertices[n_vertex].color.blue = blue;
          vertices[n_vertex].color.alpha = alpha;
          n_vertex++;

          vertices[n_vertex].coord.x = m_dAxisB[i][0];
          vertices[n_vertex].coord.y = m_dAxisB[i][1];
          vertices[n_vertex].coord.z = m_dAxisB[i][2];
          vertices[n_vertex].normal.i = m_dCapNormB[0];
          vertices[n_vertex].normal.j = m_dCapNormB[1];
          vertices[n_vertex].normal.k = m_dCapNormB[2];
          vertices[n_vertex].color.red = red;
          vertices[n_vertex].color.green = green;
          vertices[n_vertex].color.blue = blue;
          vertices[n_vertex].color.alpha = alpha;
          n_vertex++;

            }
        }

       // Compute the quads and allocate a surface element array

        surf_elements = (surf_element_c*)malloc(nc*Helix_NSEG*3*sizeof(surf_element_c));
        n_quad = 0;
        for (i=0; i<nc; i++)
        {
       for (j=0; j<Helix_NSEG; j++)
       {
          for (k=0; k<3; k++)
          {
              surf_elements[n_quad].vertex_ref[0] = i*Helix_NSEG*6+j*6+k*2;
              surf_elements[n_quad].vertex_ref[3] = i*Helix_NSEG*6+j*6+k*2+1;
              if (j == Helix_NSEG-1)
              {
                  surf_elements[n_quad].vertex_ref[1] = i*Helix_NSEG*6+k*2;
                  surf_elements[n_quad].vertex_ref[2] = i*Helix_NSEG*6+k*2+1;
              }
              else
              {
                  surf_elements[n_quad].vertex_ref[1] = i*Helix_NSEG*6+j*6+k*2+6;
                  surf_elements[n_quad].vertex_ref[2] = i*Helix_NSEG*6+j*6+k*2+7;
              }
              n_quad++;
          }
       }
        }


     */
    return true;
}

//////////////////////////////////////////////////////////////////////
bool Helix::CreateAxis(int axis_len, double caA[4][3], double caB[4][3])
{
    /* --------------------------------------------------------------------
       Function: Create the axis from the first 4 and the last 4
        residues.
       Input:    axis_len -- Number of residues in the cylinder
        caA -- First 4 CA coordinates
        caB -- Last 4 CA coordinates
       Output: m_dAxisA -- Coordinates for one end of the cylinder axis.
        m_dAxisB -- Coordinates for the other end of the cylinder axis.
       Value:    False if errors.
       -------------------------------------------------------------------- */
    // Begin ss_axis_create()

    double ptA[3], ptB[3];

    // Average for both ends
    if (axis_len == 3)
    {
        for (int i = 0; i < 3; i++)
        {
            ptA[i] = (caA[0][i] + caA[1][i])/2.0;
            ptB[i] = (caB[0][i] + caB[1][i])/2.0;
        }
    }
    else if (axis_len == 4)
    {
        for (int i = 0; i < 3; i++)
        {
            ptA[i] = (caA[0][i] + caA[1][i] + caA[2][i])/3.0;
            ptB[i] = (caB[0][i] + caB[1][i] + caB[2][i])/3.0;
        }
    }
    else
    {
        for (int i = 0; i < 3; i++)
        {
            ptA[i] = (caA[0][i] + caA[1][i] + caA[2][i] + caA[3][i])/4.0;
            ptB[i] = (caB[0][i] + caB[1][i] + caB[2][i] + caB[3][i])/4.0;
        }
    }

    // Now project last point onto axis */
    ProjectPoint(caA[0], ptA, ptB, m_dAxisA);
    ProjectPoint(caB[0], ptA, ptB, m_dAxisB);

    return true;

} // End ss_axis_create()

//////////////////////////////////////////////////////////////////////
bool Helix::ProjectPoint(double pt[3], double m_dEndA[3], double m_dEndB[3], double ptp[3])
{
    /* --------------------------------------------------------------------
       Function: Project a point onto a line.
       Input:    pt -- Point to be projected.
        m_dEndA, m_dEndB -- End points of line
       Output:   ptp -- Projected coordinates.
       Value:    False if errors.
       -------------------------------------------------------------------- */
    // Begin ss_point_project()

    double norm[3], vec[3];
    double dista, distb, distp, dist, t;
    int i;

    dista = ml_ptpt_distance(pt, m_dEndA);
    distb = ml_ptpt_distance(pt, m_dEndB);
    dist = ml_ptpt_distance(m_dEndA, m_dEndB);
    if (dista > distb)
    {
        for (i = 0; i < 3; i++)
        {
            norm[i] = m_dEndB[i] - m_dEndA[i];
            vec[i] = pt[i] - m_dEndA[i];
        }
        distp = (ml_dot(norm, vec))/dist;
        t = distp/dist;
        for (i = 0; i < 3; i++)
        {
            ptp[i] = t*norm[i] + m_dEndA[i];
        }
    }
    else
    {
        for (i = 0; i < 3; i++)
        {
            norm[i] = m_dEndA[i] - m_dEndB[i];
            vec[i] = pt[i] - m_dEndB[i];
        }
        distp = (ml_dot(norm, vec))/dist;
        t = distp/dist;
        for (i = 0; i < 3; i++)
        {
            ptp[i] = t*norm[i] + m_dEndB[i];
        }
    }

    return true;

} // End ss_point_project()

//////////////////////////////////////////////////////////////////////
bool Helix::CreateCylinder()
{
    /* --------------------------------------------------------------------
       Function: Create a cylinder with the given axis, m_dRadius, and number
        of segments.
       Input:    m_dAxisA, m_dAxisB -- Axis endpoints
        m_dRadius -- m_dRadius for the cylinder
        Helix_NSEG -- Number of segments in the cylinder
       Output:   m_dEndA -- Coordinates for one end of the cylinder
        m_dEndB -- Coordinates for the other end of the cylinder.
        m_dNormA -- Normals for one end
        m_dNormB -- Normals for the other end.
       Value:    False if errors.
       -------------------------------------------------------------------- */
    /* Begin ss_cylinder_create() */

#define pi 3.14159265

    double axisNorm[3];
    double dist, mat[3][3], invmat[3][3];
    double angle_inc, angle, a_pt[3], b_pt[3];
    double a_center[3], b_center[3];

    /* Compute the axis and a transform matrix */
    int i;
    for (i = 0; i < 3; i++)
    {
        axisNorm[i] = m_dAxisB[i] - m_dAxisA[i];
    }
    dist = ml_ptpt_distance(m_dAxisA, m_dAxisB);
    ml_normalize(axisNorm);

    ml_mat_normalTransformZ(axisNorm, invmat);
    ml_mat_transpose33(invmat, mat);

    /* Get center points for use in calculating normals */
    a_pt[0] = b_pt[0] = 0.0;
    a_pt[1] = b_pt[1] = 0.0;
    a_pt[2] = 0.0;
    b_pt[2] = dist;
    ml_vec_mat33(a_pt, mat, a_center);
    ml_vec_mat33(b_pt, mat, b_center);

    /* Generate the top and bottom circles */
    angle_inc = pi*2.0/(double)Helix_NSEG;
    angle = 0.0;
    for (i = 0; i < Helix_NSEG; i++)
    {
        a_pt[0] = b_pt[0] = m_dRadius*cos(angle);
        a_pt[1] = b_pt[1] = m_dRadius*sin(angle);
        a_pt[2] = 0.0;
        b_pt[2] = dist;
        ml_vec_mat33(a_pt, mat, m_dEndA[i]);
        ml_vec_mat33(b_pt, mat, m_dEndB[i]);
        for (int j = 0; j < 3; j++)
        {
            m_dNormA[i][j] = m_dEndA[i][j] - a_center[j];
            m_dNormB[i][j] = m_dEndB[i][j] - b_center[j];
            m_dEndA[i][j] += m_dAxisA[j];
            m_dEndB[i][j] += m_dAxisA[j];
        }
        ml_normalize(m_dNormA[i]);
        ml_normalize(m_dNormB[i]);
        angle += angle_inc;
    }

    return true;

}

void Helix::SetColor(unsigned char red, unsigned char green, unsigned char blue)
{
    rgb[0] = red;
    rgb[1] = green;
    rgb[2] = blue;
    rgb[3] = 255;
}
