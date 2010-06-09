#ifndef MI_MOLECULEBASE_H
#define MI_MOLECULEBASE_H

#include <QObject>
#include <map>

#include "model.h"
#include "Bond.h"
#include "RESIDUE_fwd.h"
#include "iterator.h"

namespace chemlib
{

    class MIMoleculeBase : public QObject
    {
        Q_OBJECT

        typedef std::map<MIMoleculeBase*, size_t> MoleculeRefCountMap;

        /**
         * Stores reference counts to objects of this class. Currently,
         * only increments and decrements in constructor and destructor.
         * Used to check if object has been deleted with isValid method.
         */
        static MoleculeRefCountMap refCounts;

        MIMoleculeBase(const MIMoleculeBase& /* mol */)
        {
        }
        MIMoleculeBase&operator=(const MIMoleculeBase& /* mol */)
        {
            return *((MIMoleculeBase*)0);
        }                                                                                      // NOTE: broken implementation just to avoid compiler warning, do not use!

    public:

        /**
         * Returns whether the given molecule is still valid (not been deleted).
         */
        static bool isValid(MIMoleculeBase *mol);

        MIMoleculeBase();
        MIMoleculeBase(Residue *reslist, const std::string &cmpd, Bond *conns, int nconns);
        virtual ~MIMoleculeBase();

#ifndef USE_ONLY_MIITER
// Use of these functions should be replaced by the MIIter<RESIDUE> versions.

        Residue *getResidues()
        {
            return residues;
        }

        Residue *getSymmResidues()
        {
            return SymmResidues;
        }

#endif
        MIIterBase<Residue> *GetResidues();

        MIIterBase<Residue> *GetSymmResidues();

        virtual MIAtom *GetAtom(int natom);

        bool IsCoordsChanged() const
        {
            return coords_changed;
        }

        void SetCoordsChanged(bool m = true)
        {
            coords_changed = m;
            if (m)
            {
                SetModified(m);
            }
        }

        bool GetModified() const
        {
            return modified;
        }

        void SetModified(bool m = true)
        {
            modified = m;
        }

        Residue *AddWater(float x, float y, float z, bool rebuild = true);

        //torsion handling
        void SetTatom(MIAtom *a1, MIAtom *a2)
        {
            Tatom1 = a1;
            Tatom2 = a2;
        }

        int SetupTorsion(MIAtom *a1, MIAtom *a2, MIAtomList *atoms);
        int SetupTorsion(MIAtom *a1, MIAtom *a2);
        void RotateTorsion(float deg);
        bool ClearTorsion();

        int Build(bool symmetryAtomsOnly = false);


        void Translate(float, float, float, MIAtomList *atoms);
        bool Contains(Residue *res);

        std::string compound;

        void ReplaceRes(Residue *oldres, Residue *dictres);
        Residue *InsertRes(Residue *atres, const Residue *dictres, int where, unsigned short chain_id = 0);
        void InsertResidues(Residue *atResidue, Residue *residues, int where, unsigned short chain_id = 0);
        void DeleteRes(Residue*);
        void DeleteResidues(std::vector<Residue*> residues);
        void DeleteAllResidues();
        void doDeleteRes(Residue *residue);
        void DeleteChain(Residue *chain);


        size_t SplitAtoms(MIAtomList &atoms, bool torsion_only);
        size_t LengthChain(unsigned short chain_id);
        MIAtom *GetTorsionAtom1() const
        {
            return Tatom1;
        }

        MIAtom *GetTorsionAtom2() const
        {
            return Tatom2;
        }

        std::vector<Bond>&getBonds()
        {
            return bonds;
        }

        std::vector<Bond>&getSymmetryBonds()
        {
            return symmetryBonds;
        }

        std::vector<Bond>::iterator BondsBegin()
        {
            return bonds.begin();
        }

        std::vector<Bond>::iterator BondsEnd()
        {
            return bonds.end();
        }

        bool AddBond(MIAtom*, MIAtom*);
        bool BreakBond(MIAtom*, MIAtom*);
        void FreeBonds();
        bool alreadybonded(MIAtom*, MIAtom*);

        void DeleteAtoms(MIAtomList atoms);
        void DeleteAtom(MIAtom *a);

        void ClearHbonds()
        {
            hbonds.clear();
        }

        void BuildHBonds();
        bool AddHBond(MIAtom*, MIAtom*);

        bool ReplaceMainChain(Residue *where, Residue *with, int nres);



        int getnresidues();

        int getnlinks()
        {
            return nlinks;
        }

        void BuildLinks();
        void ClearLinks();
        bool linked(MIAtom*, MIAtom*);
        bool linked(MIAtom*);

        bool BuildCB(Residue *res);
        bool Revert(const char *pathname);
        void FixAtomicNumbers();
        void SetSecStr(Residue *r1, Residue *r2, char sec_str);
        void Connect(Bond &connect);

        void SymmLink();
        void ClearSymmList();

        void FixChains();
        void SortChains();

        virtual void FixHeaders(std::vector<std::string>&)
        {
        }

        bool SavePDBFile(const char*);
        bool SavePDBFileVisible(const char *savepath);

        void SecStrFromAngles();
        void InitSeqPos();


        /**
         * Set the float to the torsion value if we find a match
         */
        bool GetTorsionValue(char*, Residue*);
        std::vector<std::string> *GetFileHead()
        {
            return &FileHead;
        }

        std::vector<std::string> *GetFileTail()
        {
            return &FileTail;
        }

        void setResidueNames(std::vector<Residue*> &residues, const std::string &name);
        void setChainId(Residue *chain, char c);
        void renumberChain(Residue *chain, int n);


        // General rules for deletion signalling:
        //
        //  All molecule, residue and atom deletion signalling should be handled
        //  by this class
        //
        //  Listeners to these signals should be prepared to get *ToBeDeleted
        //  signals for objects that they might not know about, since there is
        //  no creation signal
        //
        //  Signals are sent before and after a deletion event.
        //
        //  The the *Deleted (as opposed to the *ToBeDeleted) signals are called
        //  *after* the object is dead, dead, dead.  Unless you like segfaults,
        //  do not attempt to dereference the molecule pointer passed from
        //  moleculeDeleted.  Its only permissable use is as a raw address
        //  (e.g. for an index into a hash table.)
        //
        //  In general, only the highest-level appropriate signal is sent. For
        //  instance, when deleting an entire molecule, only the
        //  moleculeToBeDeleted signal is sent, not signals for each residue and
        //  or atom.
        //
        //  Therefore, if you care about atom deletion, you need to listen for
        //  residue and molecule deletion, too.  If you care about residue
        //  deletion, you need to listen to molecule deletion.
        //
        // If listeners need to know about each deleted residue or atom, they
        // may may wish to implement their signal handlers in terms of each
        // other.  For instance:
        //
        //     void moleculeToBeDeleted(MIMoleculeBase* mol) {
        //       std::vector<RESIDUE*> residues;
        //       for (MIIterBase<RESIDUE> *res=mol->GetResidues(); res; ++res) {
        //         residues.push_back(res);
        //       }
        //       for (MIIterBase<RESIDUE> *res=mol->GetSymmResidues(); res; ++res) {
        //         residues.push_back(res);
        //       }
        //       residuesToBeDeleted(mol,residues);
        //     }
        //     void residuesToBeDeleted(MIMoleculeBase *mol, std:vector<RESIDUE*> &residues) {
        //       MIAtomList atoms;
        //       for (size_t i=0; i < residues.size(); ++i) {
        //         atoms.insert(atoms.end(), residues[i].atoms.begin(), residues[i].atoms.end());
        //       }
        //       atomsToBeDeleted(mol,atoms);
        //     }
        //     void atomsToBeDeleted(MIMoleculeBase *mol, const MIAtomList &atoms) {
        //       // Do real work here.
        //     }

    signals:
        // sent when a [group of] atoms, but not an entire residue, is deleted
        void atomsToBeDeleted(chemlib::MIMoleculeBase*, const chemlib::MIAtomList&);
        void atomsDeleted(chemlib::MIMoleculeBase*);

        // sent when residue[s] deleted
        void residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&);
        void residuesDeleted(chemlib::MIMoleculeBase*);

        // sent when molecule is deleted
        void moleculeToBeDeleted(chemlib::MIMoleculeBase*);
        void moleculeDeleted(chemlib::MIMoleculeBase*);

        // sent when symmetry residues are about to be cleared
        void symmetryToBeCleared(chemlib::MIMoleculeBase*);

        // sent when atoms shown/hidden, bvalue/occ changed, or color changed
        void atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&);

        // sent when Build called, res inserted, id changed, renumber
        void moleculeChanged(chemlib::MIMoleculeBase*);

    public:
        std::vector<Bond> hbonds;

    protected:

        Residue *doInsertRes(Residue *atres, const Residue *dictres, int where, unsigned short chain_id);
        int searchbonds();
        int searchbonds(MIAtomList*);

        MIAtom *Tatom1;
        MIAtom *Tatom2;

        Residue *SymmResidues;
        Residue *residues;

        bool coords_changed;
        bool modified;
        char link_here[MAXNAME];
        char link_next[MAXNAME];
        int nlinks;
        int nresidues;

        std::vector<Bond> bonds;
        std::vector<Bond> connects;
        std::vector<Bond> symmetryBonds;
        std::vector<std::string> FileHead;
        std::vector<std::string> FileTail;

        virtual void PurgeSymmetryAtom(MIAtom*);
        virtual void PurgeAllAtoms();
        virtual void PurgeReferences(MIAtom*);
        void PurgeReferences(Residue*);
        virtual void PurgeAtom(MIAtom*);
        virtual void updateFixChainOptions(bool* /* breakByDiscontinuity */,
                                           bool* /* breakByNonpeptide */)
        {
        }

        void doDeleteAtom(MIAtom *a);

    private:
        void PurgeResidue(Residue*);
        void PurgeSymmetryResidues(Residue*);


    };

}


#endif // ifndef MI_MOLECULEBASE_H
