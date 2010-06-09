#ifndef CHIRAL_H
#define CHIRAL_H

#include "Matrix.h"
#include "Residue.h"
#include "Bond.h"
#include "CHIRALDICT.h"

#include <vector>

namespace chemlib
{
    class Ligand;

    class Chiral
    {
    public:
        void SetCenter(MIAtom*);    //Assignment Methods
        void SetOrder(int);
        void AddSub(MIAtom*);

        MIAtom *GetCenter();         //Accessor Methods
        int GetOrder();
        MIAtom *GetSub(int);
        void Clear();
        void Card(std::string&);            //Output Methods
        int Measure();

        double ideal_volume;
    protected:
        MIAtom *_center;                //Pointer to the central atom
        MIAtomList _subs;     //Pointer to the (spatially-ordered) substituent atoms
        int _order;                 //Integer representing the order of the substituents

    private: // only used by Ligand
        friend class Ligand;
        CHIRAL ToCHIRAL(const MIAtomList&, const std::vector<const MIAtom*>&);
    };

    std::string FindChiralCenters(const Monomer &res, const std::vector<Bond> &bonds);
    std::string FindChiralCenters(const Residue *res, std::vector<Bond> &bonds, bool copyChiralClasses = true);

    void GraphPotentials(const Monomer &res,
                         const std::vector<Bond> &bonds,
                         std::vector<double> &gp);
    void GraphPotentials(const Residue &res,
                         const std::vector<Bond> &bonds,
                         std::vector<double> &gp);

//	void construct_g_matrix(const Ligand &lig,
//							TNT::Matrix<double> &m);
    void construct_g_matrix(const Monomer &res,
                            const std::vector<Bond> &bonds,
                            TNT::Matrix<double> &m);
    void construct_g_matrix(const Residue &res,
                            const std::vector<Bond> &bonds,
                            TNT::Matrix<double> &m);

    void construct_c_matrix(const Ligand &lig, TNT::Matrix<double> &m);
    void construct_c_matrix(const Monomer &res, TNT::Matrix<double> &m);
    void construct_c_matrix(const Residue &res,
                            const std::vector<Bond> &bonds,
                            TNT::Matrix<double> &m);

} //namespace chemlib


#endif //CHIRAL_H
