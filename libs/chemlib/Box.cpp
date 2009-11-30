#include <vector>
#include "Box.h"
#include "MIAtom.h"

namespace chemlib
{


/////////////////////////////////////////////////////////////////////////////
// Function:    Contains
// Purpose:		Evaluate whether an atom is inside the box
// Input:       Pointer to an atom
// Output:      True or false
// Requires:
/////////////////////////////////////////////////////////////////////////////

bool Box::Contains(const MIAtom &atom) const
{
    if (atom.x() < bounds[0][0])
    {
        return false;
    }
    else if (atom.x() > bounds[0][1])
    {
        return false;
    }
    else if (atom.y() < bounds[1][0])
    {
        return false;
    }
    else if (atom.y() > bounds[1][1])
    {
        return false;
    }
    else if (atom.z() < bounds[2][0])
    {
        return false;
    }
    else if (atom.z() > bounds[2][1])
    {
        return false;
    }
    else
    {
        return true;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Origin
// Purpose:		Returns the corner of the box with the minimum coordinates
// Input:       None
// Output:      vector<float> with the coordinates
// Requires:
/////////////////////////////////////////////////////////////////////////////

std::vector<float> Box::Origin()
{
    std::vector<float> origin;

    origin.push_back(bounds[0][0]);
    origin.push_back(bounds[1][0]);
    origin.push_back(bounds[2][0]);
    return origin;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Dimensions
// Purpose:		Returns the x, y, and z dimensions of the box
// Input:       None
// Output:      vector<float> with the dimensions
// Requires:
/////////////////////////////////////////////////////////////////////////////

std::vector<float> Box::Dimensions()
{
    std::vector<float> dimensions;

    dimensions.push_back(bounds[0][1] - bounds[0][0]);
    dimensions.push_back(bounds[1][1] - bounds[1][0]);
    dimensions.push_back(bounds[2][1] - bounds[2][0]);
    return dimensions;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Constructor
// Purpose:		Create the box from a vector of ptrs to residues
// Input:       Vector of residues, plus a margin to include in all six dimensions
// Output:      Constructs a box that surrounds all the atoms in the input residues
// Requires:
/////////////////////////////////////////////////////////////////////////////
Box::Box(const std::vector<Residue*> &residues, float margin)
{

    std::vector<Residue*>::const_iterator res;
    bool first = true;                                          //Flags the first time through
    float x, y, z;
    unsigned int i;

    for (res = residues.begin(); res != residues.end(); ++res)      //Loop by residue, then by atom
    {
        for (i = 0; i < (*res)->atoms().size(); ++i)
        {

            x = (*res)->atom(i)->x();                            //Temporary storage of coords
            y = (*res)->atom(i)->y();
            z = (*res)->atom(i)->z();

            if (first)                                          //On the first time through
            {
                bounds[0][0] = x;                               //initialize all the values
                bounds[0][1] = x;                               //to the coords of the first
                bounds[1][0] = y;                               //atom
                bounds[1][1] = y;
                bounds[2][0] = z;
                bounds[2][1] = z;
                first = false;
                continue;
            }

            bounds[0][0] = (x < bounds[0][0]) ? x : bounds[0][0];       //Expand the
            bounds[0][1] = (x > bounds[0][1]) ? x : bounds[0][1];       //size of the box
            bounds[1][0] = (y < bounds[1][0]) ? y : bounds[1][0];       //to contain this
            bounds[1][1] = (y > bounds[1][1]) ? y : bounds[1][1];       //atom
            bounds[2][0] = (z < bounds[2][0]) ? z : bounds[2][0];
            bounds[2][1] = (z > bounds[2][1]) ? z : bounds[2][1];
        }
    }

    for (i = 0; i < 3; ++i)                  //Complete the box by adding a margin
    {
        bounds[i][0] -= margin;     //around the input atoms in all 6 directions
        bounds[i][1] += margin;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Expand
// Purpose:		Expands the box the necessary amount accomadate an additional atom
// Input:       Ref to an atom, plus a margin to include at the edges
// Output:      Modifies the box
// Requires:
/////////////////////////////////////////////////////////////////////////////
void Box::Expand(const MIAtom &atom, float margin)
{

    if (atom.x() - margin < bounds[0][0])             //Check x dimension
    {
        bounds[0][0] = atom.x() - margin;
    }
    else if (atom.x() + margin > bounds[0][1])
    {
        bounds[0][1] = atom.x() + margin;
    }

    if (atom.y() - margin < bounds[1][0])             //Check y dimension
    {
        bounds[1][0] = atom.y() - margin;
    }
    else if (atom.y() + margin > bounds[1][1])
    {
        bounds[1][1] = atom.y() + margin;
    }

    if (atom.z() - margin < bounds[2][0])             //Check z dimension
    {
        bounds[2][0] = atom.z() - margin;
    }
    else if (atom.z() + margin > bounds[2][1])
    {
        bounds[2][1] = atom.z() + margin;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Print
// Purpose:		Prints the coordinates of the box
// Input:       None
// Output:      Prints the coordinates to stdout
// Requires:
/////////////////////////////////////////////////////////////////////////////
//void Box::Print() {
//	for(int i=0; i<3; ++i) {
//		std::cout << bounds[i][0] << " "
//				  << bounds[i][1] << endl;
//	}
//}

}
