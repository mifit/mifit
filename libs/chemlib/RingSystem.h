#ifndef RINGSYSTEM_H
#define RINGSYSTEM_H

#include "MIAtom_fwd.h"
#include "Bond.h"
#include "Residue.h"
#include "Aromatic.h"

#include <vector>
#include <utility>
#include <algorithm>

namespace chemlib {

class LigDictionary;
class Ligand;
//@{
// A class to hold and manage the atoms and bonds in a ring system
//@}
class RingSystem {
  void AssignAromaticAtoms();
public:
  void Clear();                                     //Clears the object (no effect on the parent mol)
  void SetIndex(int);
  void SetMolecule(Ligand*);                        //Set the molecule to which this rs belongs
  void AddBond(Bond*);                              //Add a bond to the ringsys
  void AddAtom(MIAtom*);                            //Add an atom to the ringsys
  bool Contains(const MIAtom*) const;               //Check for membership of given atom
  bool Contains(const Bond*) const;                 //Check for membership of given bond
  bool IsAromatic() const;                          //True if it contains an aromatic system
  bool IsAllAromatic() const;                       //True if every atom and bond is aromatic
  void SetAllAromatic();
  void Extend(MIAtom*);
  void DetectAromatics();                           //Flags and stores aromaticity info (called once)
  void AssignAromaticBonds();                       //Flags which bonds are aromatic
  void ExtendAromatic(Aromatic&, MIAtom*, bool*);
  std::string PrintAromatics() const;
  std::string Print() const;                            //Prints a summary of the ring system

  void GeneratePlanes(LigDictionary&);              //Create pln definitions and add to dictionary
  void GenerateImpropers(LigDictionary&);           //Create imp definitions and add to dictionary
  void CyclohexImpropers(LigDictionary&) const;
  void GenerateConnTable();                         //Create the connection table
  void GenerateSRData();
  void DetectFusedRing();                           //Use the conn_table to set the _fused flag
  void Flatten();                                   //Project all the atoms into 2D
  void GetNormal(std::vector<double>& normal) const;         //Provides a normal vector to the ring plane
  void FixExocyclics();
  void GetExocyclics(std::vector<std::pair <MIAtom*, MIAtom*> > *) const;
  void GetCovExocyclics(std::vector<std::pair <MIAtom*, MIAtom*> > *);
  void ExoDirection(MIAtom* source, double* v) const;
  void LSqrPlane(double[], double*) const;              //Calc a best fit plane through all the atoms
  int SmallestRing(const MIAtom*);
  int SmallestRing(const Bond*);
  int SmallestRing(const MIAtom*, const MIAtom*, const MIAtom*);
  int SmallestAliphaticRing(const MIAtom*,
                            const MIAtom*,
                            const MIAtom*,
                            const MIAtom*) const;
  Aromatic* GetAromaticSys(const MIAtom*);
  Aromatic* GetAromaticSys(const Bond*);
  Aromatic* GetAromaticSys(const MIAtom*, const MIAtom*, const MIAtom*);

  int NumAtoms() const {
    return _atoms.size();
  }

  int NumBonds() const {
    return _bonds.size();
  }

  int NumAromatics() const {
    return _aromatics.size();
  }

  bool IsPlanar(float tolerance) const;
protected:
  int GetAtomIndex(const MIAtom*) const;            //Gets the index of an atom in the _atoms array
  void ExtendPath(int*, bool*, int, int &) const;
  //	void ExtendAromPath(int *, bool *, int, int &) const;
  void ExtendAliphaticPath(int*, bool*, int, int &) const;

  Ligand* _lig;                                     //Ptr to parent molecule
  int _rsnumber;                                    //Index of this ringsys in the molecule
  bool _fused;                                      //Flag for whether this is a "simple" ring
  MIAtomList _atoms;                      //Ptrs to atoms contained in the ringsys
  std::vector<Bond*> _bonds;                        //Ptrs to bonds contained in the ringsys
  std::vector<Aromatic> _aromatics;                 //Aromatic systems--see Aromatic class
  TNT::Matrix<bool> _conn_table;                    //Table of connectivity between the atoms
};



}   //namespace chemlib

#endif //RINGSYSTEM_H
