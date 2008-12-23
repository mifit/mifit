#ifndef MI_SURF_H
#define MI_SURF_H

class Molecule;
class Stack;
class MISurface;

namespace MISurfaceType {
const unsigned int Molecular = 0;
const unsigned int Accessible = 1;
}

namespace MISurfaceSelectionMode {
const unsigned int AtomsOnly = 0;
const unsigned int SingleResidue = 1;
const unsigned int Residues = 2;
const unsigned int Peptide = 3;
const unsigned int EntireMolecule = 4;
}

// build a new surface for the given molecules and selection
unsigned int MIMakeSurface(const std::vector<Molecule*>& mols,
                           unsigned int surface_type,
                           unsigned int atom_count = 0, // # entries in selected array
                           const unsigned int* selected = 0, // contains one entry for each atom for each molecule in mols
                           unsigned int mask = 1);


// color entire surface a single color
bool MIColorSurface(unsigned int surfnum, short color);

// color surface subset
bool MIColorSurface(unsigned int surfnum,
                    short color,
                    const std::vector<Molecule*>& mols,
                    unsigned int atom_count = 0,
                    const unsigned int* selected = 0,
                    unsigned int selected_mask = 1);

// calculate/store distance to nearest atom for each vertex
bool MISurfaceCalcDistance(unsigned int surfnum,
                           const std::vector<Molecule*>& mols,
                           unsigned int atom_count = 0,
                           const unsigned int* selected = 0,
                           unsigned int selected_mask = 1);

// color by atom type
bool MIAtomColorSurface(unsigned int surfnum,
                        const std::vector<Molecule*>& mols,
                        unsigned int atom_count = 0,
                        const unsigned int* selected = 0,
                        unsigned int selected_mask = 1);

void MIDrawSurfaces();
void MIDeleteSurface(unsigned int surfnum);

unsigned int MISurfaceCount();

void MISurfaceVectorsFromStack(std::vector<Molecule*>& mols,
                               std::vector<unsigned int>& selected,
                               unsigned int selection_type,
                               Stack* AtomStack);


// cheezy interface to avoid global variables
unsigned int MIGetSurfaceType();
void MISetSurfaceType(unsigned int);

unsigned int MIGetSurfaceSelectionMode();
void MISetSurfaceSelectionMode(unsigned int);

bool MISetSurfaceQuality(unsigned int);
unsigned int MIGetSurfaceQuality();

bool MISetSurfaceSmoothingSteps(unsigned int);
unsigned int MIGetSurfaceSmoothingSteps();

float MIGetSurfaceAlpha();
void MISetSurfaceAlpha(float percent);

int MIGetSurfaceMinColor();
void MISetSurfaceMinColor(int);

int MIGetSurfaceMaxColor();
void MISetSurfaceMaxColor(int);

float MIGetSurfaceMinValue();
void MISetSurfaceMinValue(float);

float MIGetSurfaceMaxValue();
void MISetSurfaceMaxValue(float);

bool MIGetSurfaceUseGradColor();
void MISetSurfaceUseGradColor(bool state);


// each view maintains a different set of surfaces, so we need to know which view
// is active
void MISurfaceSetCurrentView(void* view);

#endif
