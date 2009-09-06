#ifndef MI_MODEL_H
#define MI_MODEL_H

#include <cmath>
#include <vector>
#include <string>


namespace chemlib {

class RESIDUE;


// move these to "system"
#define MAXNAME 6

// types for edge order
#define NORMALBOND 0
#define SINGLEBOND 1
#define DOUBLEBOND 2
#define TRIPLEBOND 3
#define PARTIALDOUBLEBOND 4
#define HYDROGENBOND 5
#define IONICBOND 6
#define METALLIGANDBOND 7

//definitions of bond types - the bonds get sorted by depth
//and this flag lets us find each type again.
#define B_NORMAL    0
#define B_POINT     1
#define B_LINK      2
#define B_RIBBON    3
#define B_CONNECT   4
#define B_SYMM      5
#define B_SYMM_POINT 6

//Chiral order types for Tetrahedral and Allene-like chiral centers
#define COUNTERCLOCKWISE 1
#define CLOCKWISE 2

// types for edge stereochemistry
#define STEREO_NONE 0
#define STEREO_UP 1
#define STEREO_DOWN 2
#define STEREO_WEDGE 3
#define STEREO_DASHED 4
#define STEREO_INFER 5
#define STEREO_EITHER 6


//Chirality types
//The "default" is a temporary marker used when reading SMILES...this should always be resolved
//to a positive value during chirality processing.
#define CH_DEFAULT -1
#define CH_NONE 0
#define CH_TETRAHEDRAL 1
#define CH_ALLENE_LIKE 2
#define CH_TRIGONAL_BIPYRAMIDAL 3
#define CH_SQUARE_PLANAR 4
#define CH_OCTAHEDRAL 5

//Molecular geometry types
//(the arrangement of features around an atom)
#define LINEAR 1
#define TRIGONAL_PLANAR 2
#define TETRAHEDRAL 3
#define SQUARE_PLANAR 4
#define TRIGONAL_BIPYRAMIDAL 5
#define OCTAHEDRAL 6


// residue.linkage_type definitions
#define NTERMINUS ((unsigned short)1)
#define FIRST ((unsigned short)1)
#define MIDDLE ((unsigned short)2)
#define LAST ((unsigned short)4)
#define CTERMINUS ((unsigned short)4)
#define BRANCH ((unsigned short)8)
#define LINKAGEMASK ((unsigned short)15)
#define PEPTIDE ((unsigned short)16)
#define NUCLEIC ((unsigned short)32)

#define NELEMENTS   104

namespace AtomType {
/**
 * first byte used for radius type
 */
const unsigned int TYPEMASK = ((unsigned int)0x000000FF);
const unsigned int BONDED = ((unsigned int)0x00000100);
const unsigned int TORSIONATOM = ((unsigned int)0x00000200);
const unsigned int RIBBONATOM = ((unsigned int)0x00000400);
const unsigned int DELETEDATOM = ((unsigned int)0x00000800);
const unsigned int MAINCHAIN = ((unsigned int)0x00001000);
const unsigned int SIDECHAIN = ((unsigned int)0x00002000);
const unsigned int FITATOM = ((unsigned int)0x00004000);
const unsigned int REFIATOM = ((unsigned int)0x00008000);
const unsigned int INVISIBLE = ((unsigned int)0x00100000);
const unsigned int DUMMYATOM = ((unsigned int)0x00020000);
const unsigned int FREEZE = ((unsigned int)0x00040000);
const unsigned int SYMMATOM = ((unsigned int)0x00080000);
};


} //namespace chemlib

#endif //MI_MODEL_H
