#include <vector>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>

#include "CoordGenerator.h"
#include "ConfIterator.h"
#include "confgen.h"

using namespace chemlib;
using namespace std;

namespace conflib
{

void GenerateDictionary(Residue *res,
                        std::vector<Bond> &bonds,
                        std::vector<Bond> &BondLengths,
                        std::vector<ANGLE> &Angles,
                        std::vector<TORSION> &Torsions,
                        std::vector<TORSION> &Impropers,
                        std::vector<PLANE> &Planes,
                        std::vector<CHIRAL> &Chirals)
{

    MIAtomList orig_atoms = res->atoms();
    std::vector<Residue*> orig_res;
    orig_res.push_back(res);

    Ligand lig(*res, bonds);

    LigandPerceiver lp;
    //	lp.AssignAromaticBonds(&lig);
    lp.AssignHybridization(&lig);
    lp.AssignAtomGeom(&lig);
    lp.AssignChirality(&lig);

    lig.FindRingSystems();

    LigDictionary dict(&lig, (lig.residues.front()));
    dict.AnalyzeRings();
    dict.GenTorsions();
    dict.GenChirals();
    dict.GenDoubleBondImpropers();
    dict.FindAcycPlanes();
    dict.SetMolGeometry();

    CovalentGeometry cg(&lig, (lig.residues.front()));
    cg.AssignResidue();


    //	if (bonds[i].getAtom1() == BondLengths[j].getAtom1() &&
    //				bonds[i].getAtom2() == BondLengths[j].getAtom2() ||
    //				bonds[i].getAtom1() == BondLengths[j].getAtom2() &&
    //				bonds[i].getAtom2() == BondLengths[j].getAtom1())

    lig.GetConstraints(orig_atoms,
                       orig_res,
                       BondLengths,
                       Angles,
                       Torsions,
                       Impropers,
                       Planes,
                       Chirals);

    unsigned int j;
    vector<Bond>::iterator bnd;
    pair<MIAtom*, MIAtom*> ap;
    for (j = 0; j < BondLengths.size(); ++j)
    {
        ap.first = BondLengths[j].getAtom1();
        ap.second = BondLengths[j].getAtom2();
        bnd = find_if(bonds.begin(), bonds.end(), std::bind2nd(ContainsAtoms(), ap));
        if (bnd != bonds.end())
        {
            bnd->ideal_length = BondLengths[j].ideal_length;
        }
    }
}

int GenerateEnsemble(Residue *res,
                     std::vector<Bond> &bonds,
                     MIMolDictionary *dictionary,
                     bool replace)
{
    int nConfs = 0;

    //Create a confsaver to store the new conformations
    ConfSaver tmp(res);

    vector<TORSION> torsions;
    dictionary->GetFlexibleTorsions(torsions, res);

    conflib::ConfEnumerator ce(res, bonds, torsions);
    if (ce.NumberConfs() < 1000)
    {
        nConfs = ce.GenerateConfs(tmp, 1000);
    }
    else
    {
        conflib::ConfSampler cs(res, bonds, torsions);
        nConfs = cs.GenerateConfs(tmp, 1000);
    }

    if (tmp.NumberSets() == 0)
    {
        return 0;
    }
    /*
        if (replace) {
            dictionary->DeleteConformers(res->type());
        }
     */
    dictionary->AddConfs(tmp, replace);
    return nConfs;
}

int GenerateEnsemble(Residue *res,
                     MIMoleculeBase *model,
                     std::vector<Bond> &bonds,
                     std::vector<TORSION> &torsions,
                     GeomSaver &confs)
{
    int nConfs = 0;

    //Create a confsaver to store the new conformations
    ConfSaver tmp(res);

    //	vector<TORSION> torsions;
    //	dictionary->GetFlexibleTorsions(torsions, res);

    conflib::ConfEnumerator ce(res, bonds, torsions);
    if (ce.NumberConfs() < 1000)
    {
        nConfs = ce.GenerateConfs(tmp, 1000);
    }
    else
    {
        conflib::ConfSampler cs(res, bonds, torsions);
        nConfs = cs.GenerateConfs(tmp, 1000);
    }

    if (tmp.NumberSets() == 0)
    {
        return 0;
    }

    tmp.ConvertToGeomSaver(confs, model);
    return nConfs;
}

void GenerateCoordinates(Residue *res,
                         const std::vector<Bond> &bonds,
                         std::string&)
{

    //Construct a ligand object
    Ligand lig(*res, bonds);
    Monomer *new_res = (lig.residues.front());

    //Assign some basic features of the molecule
    LigandPerceiver lp;
    lp.AssignHybridization(&lig);
    lp.AssignAtomGeom(&lig);
    lp.AssignChirality(&lig);

    //Determine which atoms and bonds are cyclic and cluster
    //them into distinct ring systems
    lig.FindRingSystems();

    //Determine the constraints (planes, torsions, & chirals)
    //for the molecule
    LigDictionary dict(&lig, (lig.residues.front()));
    dict.AnalyzeRings();
    dict.GenTorsions();
    dict.GenChirals();
    dict.GenDoubleBondImpropers();
    dict.FindAcycPlanes();
    dict.SetMolGeometry();

    //Assign target values for the bond lengths and angles
    CovalentGeometry cg(&lig, (lig.residues.front()));
    cg.AssignResidue();

    int i;
    for (i = 0; i < res->atomCount(); ++i)
    {
        res->atom(i)->copyChemicalData(*new_res->atom(i));
    }
    //Construct a coordinate generator
    CoordGenerator genXYZ(&lig, (lig.residues.front()));

    //Generate the coordinates with the new "stochastic distance geometry" method
    /* double sdgScore = */ genXYZ.sdgRefine();

    for (i = 0; i < res->atomCount(); ++i)
    {
        res->atom(i)->copyPosition(*new_res->atom(i));
    }

    //Now try getting coordinates the old way.
    //	double buildScore = genXYZ.GenMolStruct(log);

    //If we did better this time, replace the coordinates
    //	if (buildScore < sdgScore) {
    //		for(i=0; i<res->atomCount(); ++i) {
    //			res->atom(i)->copyPosition(*new_res->atom(i));
    //		}
    //	}

}

void GenerateCoordinates(Ligand *lig, std::string &log)
{


    //Assign some basic features of the molecule
    //	LigandPerceiver lp;
    //	lp.AssignAromaticBonds(&lig);
    //	lp.AssignHybridization(lig);
    //	lp.AssignAtomGeom(lig);
    //	lp.AssignChirality(lig);

    //Determine which atoms and bonds are cyclic and cluster
    //them into distinct ring systems
    //	lig->FindRingSystems();

    //Determine the constraints (planes, torsions, & chirals)
    //for the molecule
    LigDictionary dict(lig, (lig->residues.front()));
    dict.AnalyzeRings();
    dict.GenTorsions();
    dict.GenChirals();
    dict.GenDoubleBondImpropers();
    dict.FindAcycPlanes();
    dict.SetMolGeometry();

    //Assign target values for the bond lengths and angles
    CovalentGeometry cg(lig, (lig->residues.front()));
    cg.AssignResidue();

    //Construct a coordinate generator
    CoordGenerator genXYZ(lig, (lig->residues.front()));

    //Generate the coordinates
    genXYZ.GenMolStruct(log);
}

void sdgGenerateCoordinates(Ligand *lig)
{


    //Assign some basic features of the molecule
    //	LigandPerceiver lp;
    //	lp.AssignAromaticBonds(&lig);
    //	lp.AssignHybridization(lig);
    //	lp.AssignAtomGeom(lig);
    //	lp.AssignChirality(lig);

    //Determine which atoms and bonds are cyclic and cluster
    //them into distinct ring systems
    //	lig->FindRingSystems();

    //Determine the constraints (planes, torsions, & chirals)
    //for the molecule
    LigDictionary dict(lig, (lig->residues.front()));
    dict.AnalyzeRings();
    dict.GenTorsions();
    dict.GenChirals();
    dict.GenDoubleBondImpropers();
    dict.FindAcycPlanes();
    dict.SetMolGeometry();

    //Assign target values for the bond lengths and angles
    CovalentGeometry cg(lig, (lig->residues.front()));
    cg.AssignResidue();

    //Construct a coordinate generator
    CoordGenerator genXYZ(lig, (lig->residues.front()));

    //Generate the coordinates
    genXYZ.sdgGenMolStruct();
}

void sdgRefine(Ligand *lig)
{


    //Assign some basic features of the molecule
    //	LigandPerceiver lp;
    //	lp.AssignAromaticBonds(&lig);
    //	lp.AssignHybridization(lig);
    //	lp.AssignAtomGeom(lig);
    //	lp.AssignChirality(lig);

    //Determine which atoms and bonds are cyclic and cluster
    //them into distinct ring systems
    //	lig->FindRingSystems();

    //Determine the constraints (planes, torsions, & chirals)
    //for the molecule
    LigDictionary dict(lig, (lig->residues.front()));
    dict.AnalyzeRings();
    dict.GenTorsions();
    dict.GenChirals();
    dict.GenDoubleBondImpropers();
    dict.FindAcycPlanes();
    dict.SetMolGeometry();

    //Assign target values for the bond lengths and angles
    //	CovalentGeometry cg(lig, &(lig->residues.front()));
    //	cg.AssignResidue();

    //Construct a coordinate generator
    CoordGenerator genXYZ(lig, (lig->residues.front()));

    //Generate the coordinates
    genXYZ.sdgRefine();
}

} //namespace conflib


