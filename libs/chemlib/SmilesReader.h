#ifndef SMILESREADER_H
#define SMILESREADER_H

#include <vector>
#include <stack>

#include "Ligand.h"
#include "LigandPerceiver.h"
#include "mol_util.h"
#include "CovalentGeom.h"

namespace chemlib {

RESIDUE* SmilesToMol(const std::string& smiles,
                     std::vector<Bond>& bonds,
                     std::string& ErrorMessage);

bool SmilesToLig(const std::string& smiles,
                 Ligand** lig,
                 std::string& ErrorMessage);

bool ValidateSmiles(const std::string& smiles, std::string& ErrorMessage);

class SmilesReader;
class Ligand;

class SmiRingClosures {
public:
  friend class SmilesReader;
  SmiRingClosures();            //Constructor
  ~SmiRingClosures();           //Destructor
  void Create(int*);        //Form the first half of a ring-closing bond
  void Complete(int);       //Complete the second half of an r-c bond
  void ProcessSimple(char);         //Receive an r-c bond request in single-digit form
  void ProcessComplex(char*);       //Receive an r-c bond request in "%nn" form
private:
  SmilesReader* _reader;
  int _short_half_bonds[100];      //Indices of the halfbond in the bonds vector in the mol obj
  bool _open_half_bonds[100];      //Table to track current ring closures
  int _ringnum;      //Number of the current half-bond (0-9 for single-digit,10-99 for long)
};


class SmilesReader {
public:
  SmilesReader(const char*, Ligand*);
  SmilesReader(const char*, Ligand*, RESIDUE&, std::vector<Bond>& bonds);
  ~SmilesReader();
  void TraverseString();
  int CountAtoms();
private:
  void AddBasicAtom(int atomic_number, bool isaromatic = false);
  void ProcessBracketAtom();
  void AddAtom(int, bool, int, int, int, int, int);
  int GetNumber();
  int CountRepetitions(char);
  friend class SmiRingClosures;
  SmiRingClosures _ring_closer;

  int _natoms;      //Total number of atoms in the Ligand
  char* _smistr;       //Input SMILES string
  Ligand* _lig;       //Class for writing Ligand into
  std::stack <MIAtom*> _roots;       //MIAtom stack to keep track of last atom read
  //from all branchs;
  char* _cur_pos;       //Current position in string;
  unsigned char _bond_order;      //Current bond order (see Xguicryst.h)
  char _bond_stereo;      //Current bond stereo flags
  MIAtom _atom;
  MIAtomList _atoms;
  Residue* _res;
};

} //namepsace chemlb

#endif  //SMILESREADER_H
