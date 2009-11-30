#ifndef AROMATIC_H
#define AROMATIC_H

#include "Matrix.h"
#include "sequence_util.h"
#include "Dictionary.h"

namespace chemlib
{

    class MIAtom;
    class Bond;

//@{
// A class to hold and manage the atoms and bonds in an aromatic system
//@}
    class Aromatic
    {
    public:
        //Methods that modify the object
        void Clear();                               //Clears the object (no effect on the parent mol)
        void AddAtom(MIAtom*);                      //Add an atom to the ringsys
        void AddBond(Bond*);                        //Add a bond to the ringsys
        void GenerateConnTable();                   //Create the connection table
        void GenerateSRData();                      //Calc smallest_aromatic_ring values

        //Methods that access the object
        bool Contains(const MIAtom*) const;         //Check for membership of given atom
        bool Contains(const Bond*) const;           //Check for membership of given bond
        int NumAtoms() const;                       //Number of atoms in the ringsys
        int NumBonds() const;                       //Number of bonds in the ringsys

        void Print(std::string &s) const;                //Appends a summary to a string

        //Methods that generate Dictionary data
        void GeneratePlane(LigDictionary&) const;   //Adds a plane to the dictionary for this aromatic
        void GenerateImpropers(LigDictionary&) const; //Adds impropers to the dictionary

        //Ring-finding
        int SmallestRing(const MIAtom*) const;
        int SmallestRing(const Bond*) const;
        int SmallestRing(const MIAtom*, const MIAtom*, const MIAtom*) const;
        int SmallestRing(const MIAtom*, const MIAtom*, const MIAtom*, const MIAtom*) const;
    protected:
        //Data members
        MIAtomList _atoms;                    //Ptrs to atoms contained in the aromatic sys
        std::vector<Bond*> _bonds;                  //Ptrs to bonds contained in the aromatic sys
        TNT::Matrix<bool> _conn_table;              //Table of connectivity between the atoms

        //Utility methods
        int GetAtomIndex(const MIAtom*) const;      //Gets the index of an atom in the _atoms array
        void ExtendPath(int*, bool*, int, int&) const;
    };
}   //namespace chemlib
#endif //AROMATIC_H

