#ifndef MIFIT_MODEL_MOLECULE_H_
#define MIFIT_MODEL_MOLECULE_H_

#if defined __GNUC__
// only way to disable warning: 'class chemlib::MIMoleculeBase' should be explicitly initialized in the copy constructor
#pragma GCC system_header
#endif

#include <QObject>
#include <vector>
#include <deque>

#include <chemlib/MIMoleculeBase.h>
#include <chemlib/chemlib.h>

class Molecule;

#include "ViewPoint.h"
#include "Annotation.h"
#include "MoleculeType.h"
#include "CONTACT.h"

class CMapHeaderBase;
class SecondaryStructure;
class CMapHeaderBase;

class MolPrefsHandler
{
public:
    virtual ~MolPrefsHandler()
    {
    }

    virtual void operator()(bool *breakByDiscontinuity, bool *breakByNonpeptide) = 0;
};
void MISetMolPrefsHandler(MolPrefsHandler *mph);


class Molecule : public chemlib::MIMoleculeBase
{
    Q_OBJECT

public:
    typedef std::vector<Annotation*> AnnotationList;
    typedef std::vector<ATOMLABEL*> AtomLabelList;
    typedef std::vector<SURFDOT> SurfaceDots;

private:
    friend class MoleculeXmlHandler;

    std::deque<chemlib::MIAtom*> ribbonatoms;
    SurfaceDots dots;
    /**
     * 0= invisible else visible
     */
    int visible;
    int dots_visible;
    int labels_visible;
    int annots_visible;
    float srfdotsper;
    float srfboxsize;
    std::vector<short> savecolors;
    bool undoable;
    CMapHeaderBase *mapheader;

    bool drawBox;
    bool symmatoms_visible;
    float symm_radius;
    float symm_center[3];
    std::string alt_seq;

    SecondaryStructure *m_pSecondaryStructure;
    bool secStruc_ribbon;
    bool secStruc_schematic;

    AnnotationList annotations;
    AtomLabelList atomLabels;

    void doAtomLabelDelete(ATOMLABEL *label);
    void doAnnotationDelete(Annotation *annotation);

    Molecule(const Molecule& /* mol */)
    {
    }
    Molecule&operator=(const Molecule& /* mol */)
    {
        return *((Molecule*)0);
    }                                                                        // NOTE: broken implementation just to avoid compiler warning, do not use!

public:

    int nribbons;
    int modelnumber;
    /**
     * the full file name
     */
    std::string pathname;
    std::string author;
    std::string source;
    int ribbon_coloring;
    /**
     * to track current selections
     */
    int s_main;
    /**
     * to track current selections
     */
    int s_sides;
    /**
     * to track current selections
     */
    int s_nonprotein;
    /**
     * to track current selections
     */
    int s_link;
    /**
     * to track current selections
     */
    int s_waters;
    /**
     * to track current selections
     */
    int s_radiustype;

    int ModelType;
    bool HVisible;
    std::vector<CONTACT> hbondContacts;

    Molecule(chemlib::Residue *reslist, std::string cmpd, FILE *fp, chemlib::Bond *conns, int nconns, int type = MoleculeType::PDB);
    Molecule(int type = MoleculeType::Other);
    ~Molecule();

    void FixHeaders(std::vector<std::string> &headers);

    void TranslateResidueToCenter(chemlib::Residue *res, ViewPoint *viewpoint);
    void TranslateAtomsToCenter(std::vector<chemlib::MIAtom*>&, ViewPoint *viewpoint);

    AtomLabelList&getAtomLabels();

    /**
     * Adds an atom label to this molecule. The molecule takes ownership of the atom labels and
     * is responsible for deleting them.
     */
    void addAtomLabel(ATOMLABEL *label);
    void updateAtomLabels();
    void deleteAtomLabel(ATOMLABEL *label);
    void clearAtomLabels();
    void labelAtom(chemlib::MIAtom *atom, chemlib::Residue *res);
    void labelAtomStyle(chemlib::MIAtom *atom, int style);
    void unlabelAtom(chemlib::MIAtom *atom);
    void labelEveryNthResidue(int n);
    bool isAtomLabeled(chemlib::MIAtom *atom);
    ATOMLABEL *findLabelForAtom(chemlib::MIAtom *atom);
    void setAtomLabelVisible(ATOMLABEL *label, bool visible);
    void setAtomLabelColor(ATOMLABEL *label, unsigned char red, unsigned char green, unsigned char blue);
    void setAtomLabelText(ATOMLABEL *label, const char *text);

    bool symmAtomsVisible() const;
    void setSymmAtomsVisible(bool value);

signals:
    void atomLabelAdded(Molecule*, ATOMLABEL*);
    void atomLabelChanged(Molecule*, ATOMLABEL*);
    void atomLabelToBeDeleted(Molecule*, Molecule::AtomLabelList);
    void atomLabelDeleted(Molecule*);

public:
    AnnotationList&getAnnotations();

    /**
     * Adds an annotation to this molecule. The molecule takes ownership of the annotations and
     * is responsible for deleting them.
     */
    void addAnnotation(Annotation *annotation);

    /**
     * Deletes an annotation of this molecule.
     */
    void deleteAnnotation(Annotation *annotation);

    void addAnnotation(const std::string &s, float x = 0, float y = 0, float z = 0, int id = 0);

    /**
     * Deletes all geometry annotations.
     */
    void clearGeomAnnotations();

    /**
     * Deletes all annotations
     */
    void clearAnnotations();

signals:
    void annotationAdded(Molecule*, Annotation*);
    void annotationToBeDeleted(Molecule*, Annotation*);
    void annotationDeleted(Molecule*);

public:
    void ShowAnnotations();
    void HideAnnotations();
    int AnnotationsVisible();
    void DrawAnnotationBox(bool on);
    bool isDrawAnnotationBox();

    chemlib::Residue *MatchPentamer(std::string &pentdir, chemlib::Residue *start);

    chemlib::MIAtom *GetAtom(int natom);

    bool CheckCenter(float x, float y, float z);

    void GenSymmAtoms(ViewPoint*);


    CMapHeaderBase&GetMapHeader();
    void SetMapHeader(const CMapHeaderBase &mh);

    void PurgeAllAtoms();
    void PurgeAtom(chemlib::MIAtom*);
    void PurgeSymmetryAtom(chemlib::MIAtom*);
    void PurgeReferences(chemlib::MIAtom*);

    int SequenceIdentities();
    void DeleteLowerGap(chemlib::ResidueListIterator gap_point);
    void InsertLowerGap(chemlib::ResidueListIterator gap_point);
    void DeleteGap(chemlib::ResidueListIterator res);
    char GetSeq(int index);
    void InsertGap(chemlib::ResidueListIterator res);
    int WriteSequence(std::string path, int type);
    void ReadSequence(std::string path, int type, int skiplines);
    void SetSequence(std::string);
    std::string SeqString();
    unsigned int SeqLength();
    void SetSeq(char, int);
    void ShowHydrogens(bool);
    void ToggleHydrogens();
    bool UnDoable(Molecule *node);
    void Do();
    void UnDo();
    int getcolor(chemlib::ResidueListIterator, chemlib::MIAtom*, bool, int, int, std::string);
    long getndots();
    SURFDOT *GetDots();
    SurfaceDots&getDots();
    void setDotsColor(int color);
signals:
    void surfaceChanged(Molecule*);

public:
    long getnribboatomCount();
    void GetPDBInfo(FILE*);

    bool VisibleBounds(ViewPoint*, float&, float&, float&, float&, float&, float&);

    void Show();
    void Hide();
    void ShowAll();
    void HideAll();
    int Visible();
    void HideDots();
    void ShowDots();
    int DotsVisible();
    void ShowLabels();
    void HideLabels();
    int LabelsVisible();

    void Select(int, int, int, int, std::string, std::string, chemlib::Residue*, chemlib::Residue*,
                bool, int, int, int, int, int listtype = 0, int radiustype = 0);
    void FreeDots();

    void ClearRibbons();
    void BuildRibbons();
    chemlib::MIAtom *AllocRibbonAtom();
    long SurfaceCenter(ViewPoint*, float, float, bool ignore_hidden = true);
    long SolventSurface(ViewPoint*, float);
    long SurfaceAtom(chemlib::MIAtom*, float, bool ignore_hidden = true);
    long SurfaceAroundAtom(chemlib::MIAtom*, float, float);
    long SurfaceResidue(chemlib::Residue*, float, bool ignore_hidden = true);
    long SurfaceResidues(float, bool ignore_hidden = true);
    long Surface(chemlib::MIAtom*, bool ignore_hidden = true, bool send_signal = true);
    void Save(XMLArchive&);
    void Load(FILE *fp);

    void Translate(float, float, float, std::vector<chemlib::MIAtom*> *atoms); // need this b/c parent version is hidden by override below
    void Translate(float, float, float, std::vector<chemlib::MIAtom*> *atoms, SurfaceDots *dots);

    void Center(int &x, int &y, int &z);
    void Rotate(float rx, float ry, float rz, float cx, float cy, float cz, ViewPoint*, std::vector<chemlib::MIAtom*> *atoms, SurfaceDots *dots = NULL);

    int SaveSymmMolecule(chemlib::MIAtom *symatom, FILE *fp);
    int SeqMax();
    void MakeSecondaryStructure(bool bRibbon, bool bSchematic);
    void DeleteSecondaryStructure();
    SecondaryStructure *getSecondaryStructure();
    bool isSecondaryStructureRibbons();
    bool isSecondaryStructureSchematic();

    void toggleChainHidden(chemlib::Residue *chain);
    void toggleAtomHidden(chemlib::MIAtom *atom);
    void toggleResidueHidden(chemlib::Residue *residue);
    void setResidueColor(chemlib::Residue *residue, int color, int colorMethod);

    void toggleAtomsHidden(std::vector<chemlib::MIAtom*> &atoms);
    void toggleResiduesHidden(std::vector<chemlib::Residue*> &residues);
    void setResiduesColor(std::vector<chemlib::Residue*> &residues, int color, int colorMethod);

    void setAtomBValueAndOccupancy(chemlib::MIAtom *atom, float bvalue, float occ);
    void setAtomsBValueAndOccupancy(std::vector<chemlib::MIAtom*> atoms, float bvalue, float occ);
    void setAtomColor(chemlib::MIAtom *atom, int color);
    void setAtomsColor(std::vector<chemlib::MIAtom*> atoms, int color);
    void setChainColor(chemlib::Residue *chain, int color, int colorMethod);

    void setColor(int color, int colorMethod);

protected:

    void doAtomColor(chemlib::MIAtom *atom, int color);
    void doAtomBValueAndOccupancy(chemlib::MIAtom *atom, float bvalue, float occ);

    virtual void updateFixChainOptions(bool *breakByDiscontinuity, bool *breakByNonpeptide);

};

inline void Molecule::DrawAnnotationBox(bool on)
{
    drawBox = on;
}

inline bool Molecule::isDrawAnnotationBox()
{
    return drawBox;
}

inline CMapHeaderBase&Molecule::GetMapHeader()
{
    return *mapheader;
}

inline unsigned int Molecule::SeqLength()
{
    return alt_seq.length();
}

inline void Molecule::ToggleHydrogens()
{
    if (HVisible)
    {
        HVisible = false;
    }
    else
    {
        HVisible = true;
    }
    ShowHydrogens(HVisible);
}

inline bool Molecule::UnDoable(Molecule *node)
{
    return (undoable && node == this);
}

inline long Molecule::getndots()
{
    return dots.size();
}

inline SURFDOT*Molecule::GetDots()
{
    return &(dots[0]);
}

inline Molecule::SurfaceDots&Molecule::getDots()
{
    return dots;
}

inline long Molecule::getnribboatomCount()
{
    return ribbonatoms.size();
}

inline int Molecule::Visible()
{
    return visible;
}

inline void Molecule::HideDots()
{
    dots_visible = 0;
    surfaceChanged(this);
}

inline void Molecule::ShowDots()
{
    dots_visible = 1;
    surfaceChanged(this);
}

inline int Molecule::DotsVisible()
{
    return dots_visible;
}

inline void Molecule::ShowLabels()
{
    labels_visible = 1;
}

inline void Molecule::HideLabels()
{
    labels_visible = 0;
}

inline int Molecule::LabelsVisible()
{
    return labels_visible;
}

inline void Molecule::ShowAnnotations()
{
    annots_visible = 1;
}

inline void Molecule::HideAnnotations()
{
    annots_visible = 0;
}

inline int Molecule::AnnotationsVisible()
{
    return annots_visible;
}

inline bool Molecule::isSecondaryStructureRibbons()
{
    return secStruc_ribbon;
}

inline bool Molecule::isSecondaryStructureSchematic()
{
    return secStruc_schematic;
}

#endif /*MIFIT_MODEL_MOLECULE_H_*/
