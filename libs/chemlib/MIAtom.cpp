#include <vector>
#include <string>

#include <nongui/nonguilib.h>
#include <util/utillib.h>

#include "model.h"
#include "MIAtom.h"
#include "mol_util.h"

namespace chemlib {

//#define MEMORY_CORRUPTION_DEBUG
#if defined(MEMORY_CORRUPTION_DEBUG)
static unsigned int magic_delete=UINT_MAX;
static unsigned int DELETION_COUNT=0;
static std::map<MIAtom*,unsigned int> DELETED_ATOMS;
#endif

static const unsigned int NTYPES = 9;

float MIAtom::MIAtomRadiusForType(int i) {
  const float radius[] = { 1.85F, 1.60F, 1.70F, 1.20F, 1.90F, 1.95F, 1.20F, 1.20F, 1.20F};
  if (i < 0 || i >= (int)NTYPES) {
    return radius[0];
  }

  return radius[i];
}

int MIAtom::MIGetAtomTypeFromName(const char* name) {
  const char atype[] = { 'C', 'O', 'N', 'H', 'S', 'F', '1', '2', '3'};

  for (unsigned int i = 0; i < NTYPES; i++) {
    if (name[0] == atype[i]) {
      return i;
    }
  }

  return 0;
}

MIAtom::AtomRefCountMap MIAtom::refCounts;

bool MIAtom::isValid(const MIAtom* atom) {
  bool result = false;
  if (atom != NULL && refCounts.find(atom) != refCounts.end()) {
    result = refCounts[atom] > 0;
  }
#ifdef MEMORY_CORRUPTION_DEBUG
  if (!result && atom) {
    printf("Invalid atom %p queried, was deletion #%d\n",atom,DELETED_ATOMS[atom]);
  }
#endif
  return result;
}

/**
 * Creates an atom. Note: Doesn't allocate space for anisotropic U's
 */
MIAtom::MIAtom()
  : type_(0),
  BValue_(15.0F),
  occ_(1.0F),   //Default to full occupancy
  color_(0),
  x_(0.0),
  y_(0.0),
  z_(0.0),
  dx_(0.0),
  dy_(0.0),
  dz_(0.0),
  weight_(0),
  radius_type_(0),
  symmop_(0),
  altloc_(' '),
  atomnumber_(0),
  atomicnumber_(6),
  isaromatic_(0),
  hybrid_(3),
  geom_(0),
  mass_(0),
  formal_charge_(0),
  charge_(0.0F),
  hcount_(0),
  iscyclic_(0),
  ring_system_(-1),
  smallest_aromatic_ring_(0),
  search_flag_(0),
  chiral_class_(0),
  chiral_order_(-1),
  U_(NULL)
 {
  ++refCounts[this];
  name_[0] = '\0';
  bondnumbers_.clear();
  nabors_.clear();
}

MIAtom::~MIAtom() {
  deleteAnisotropicity();
  --refCounts[this];
  if (refCounts[this] == 0) {
    refCounts.erase(this);
  }
#if defined(MEMORY_CORRUPTION_DEBUG)
  if (DELETION_COUNT==magic_delete)
    printf("Set magic_delete in debugger and set breakpoint here.\n");
  DELETED_ATOMS[this]=DELETION_COUNT++;
#endif
}


bool MIAtom::hasAnisotropicity() const {
  return U_ != NULL;
}

void MIAtom::newAnisotropicity() {
  if (hasAnisotropicity()) {
    deleteAnisotropicity();
  }
  U_ = (float*) malloc(6*sizeof(float));
  memset(U_, 0, 6*sizeof(float));
}

void MIAtom::deleteAnisotropicity() {
  if (U_) {
    free(U_);
    U_ = NULL;
  }
}

float MIAtom::U(size_t index) const {
  MI_ASSERT(hasAnisotropicity());
  float result = 0.0f;
  if (hasAnisotropicity() && index < 6) {
    result = U_[index];
  }
  return result;
}

float MIAtom::U(size_t index, float value) {
  MI_ASSERT(hasAnisotropicity());
  float result = 0.0f;
  if (hasAnisotropicity() && index < 6) {
    U_[index] = value;
    result = U_[index];
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    CountCyclicBonds
// Purpose:		Counts the number of ring bonds to the atom
// Input:       None, uses nabors data field
// Output:      The number of ring bonds
// Requires:	FindRingSystems() has been called
/////////////////////////////////////////////////////////////////////////////

//CHECK THIS...This is really the # of cyclic neighbors, since two cyclic
//atoms can be connected by an acyclic bond! -KWB 3/22/05

int MIAtom::CountCyclicBonds() const {

  if (!iscyclic_) {
    return 0;
  }

  MIAtom_const_iter nabor;
  int n_cyc_bonds = 0;

  for (nabor = nabors_.begin(); nabor != nabors_.end(); ++nabor) {
    if ((*nabor)->iscyclic_) {
      n_cyc_bonds++;
    }
  }

  return n_cyc_bonds;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ClearHybrid
// Purpose:		Sets the atom to a default bonding state, usually coincident
//				with setting all the bonds to single, and usually preceding
//				a call to GuessBondOrders, which uses 3D coordinates to ascertain
//				the bond orders and hybridization.
// Input:       None
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void MIAtom::ClearHybrid() {
  hybrid_ = 3;                               //sp3 is always the default, and
  isaromatic_ = 0;                           //we make no assumptions about
  hcount_ = 0;                               //aromaticity or hydrogens
}

float MIAtom::getRadius() const {
  return MIAtomRadiusForType(MIGetAtomTypeFromName(name_));
}

std::string chemlib::MIAtom::chiralClassToString(int chiralClass) {
  std::string chiralString;
  switch (chiralClass) {
    case CH_DEFAULT:
      chiralString = "default";
      break;
    case CH_NONE:
      chiralString = "none";
      break;
    case CH_TETRAHEDRAL:
      chiralString = "tetrahedral";
      break;
    case CH_ALLENE_LIKE:
      chiralString = "allene-like";
      break;
    case CH_TRIGONAL_BIPYRAMIDAL:
      chiralString = "trigonal bipyramidal";
      break;
    case CH_SQUARE_PLANAR:
      chiralString = "square planar";
      break;
    case CH_OCTAHEDRAL:
      chiralString = "octahedral";
      break;
    default:
      chiralString = "other";
      break;
  }
  return chiralString;
}

void MIAtom::getPosition(double a[3]) const {
  a[0] = x_;
  a[1] = y_;
  a[2] = z_;
}

void MIAtom::getPosition(float a[3]) const {
  a[0] = x_;
  a[1] = y_;
  a[2] = z_;
}

void MIAtom::copyPosition(const MIAtom& from) {
  x_ = from.x_;
  y_ = from.y_;
  z_ = from.z_;
}

void MIAtom::copyChemicalData(const MIAtom& source) {
  hybrid_ = source.hybrid_;
  geom_ = source.geom_;
  mass_ = source.mass_;
  charge_ = source.charge_;
  formal_charge_ = source.formal_charge_;
  hcount_ = source.hcount_;
  chiral_class(source.chiral_class());
  chiral_order(source.chiral_order());
  iscyclic_ = source.iscyclic_;
  smallest_aromatic_ring_ = source.smallest_aromatic_ring_;
}

void MIAtom::copyShallow(const MIAtom& atom) {
  setName(atom.name_);
  type_ = atom.type_;
  BValue_ = atom.BValue_;
  occ_ = atom.occ_;
  color_ = atom.color_;

  copyPosition(atom);

  dx_ = atom.dx_;
  dy_ = atom.dy_;
  dz_ = atom.dz_;
  weight_ = atom.weight_;
  radius_type_ = atom.radius_type_;
  symmop_ = atom.symmop_;
  altloc_ = atom.altloc_;
  atomnumber_ = atom.atomnumber_;
  atomicnumber_ = atom.atomicnumber_;
  isaromatic_ = atom.isaromatic_;

  copyChemicalData(atom);

  ring_system_ = atom.ring_system_;
  search_flag_ = atom.search_flag_;


  if (atom.hasAnisotropicity()) {
    newAnisotropicity();
    memcpy(U_, atom.U_, 6*sizeof(float));
  }
}


bool MIAtom::MIIsMainChainAtom(const MIAtom* atom) {
  if (!strcmp(atom->name_, "C")) {
    return true;
  }
  if (!strcmp(atom->name_, "O")) {
    return true;
  }

  if (!strcmp(atom->name_, "CA")) {
    return true;
  }
  if (!strcmp(atom->name_, "N")) {
    return true;
  }
  if (!strcmp(atom->name_, "H")) {
    return true;
  }
  return false;
}

bool MIAtom::MIIsSideChainAtom(const MIAtom* atom) {
  if (!strcmp(atom->name_, "C")) {
    return false;
  }
  if (!strcmp(atom->name_, "CA")) {
    return false;
  }
  if (!strcmp(atom->name_, "O")) {
    return false;
  }
  if (!strcmp(atom->name_, "N")) {
    return false;
  }
  if (!strcmp(atom->name_, "H")) {
    return false;
  }
  return true;
}

bool MIAtom::MIIsMainChainDNAAtom(const MIAtom* atom) {
  return (!MIIsSideChainDNAAtom(atom));
}

bool MIAtom::MIIsSideChainDNAAtom(const MIAtom* atom) {
  if (atom->name_[strlen(atom->name_)-1] == '*') {
    return 0;
  }
  if (!strcmp(atom->name_, "P")) {
    return 0;
  }
  if (!strcmp(atom->name_, "O1P")) {
    return 0;
  }
  if (!strcmp(atom->name_, "O2P")) {
    return 0;
  }
  return true;
}

bool MIAtom::MIIsHydrogen(const MIAtom* atom) {
  if (atom->name_[0] == 'H'
      || atom->name_[0] == '1'
      || atom->name_[0] == '2'
      || atom->name_[0] == '3'
      || atom->name_[0] == '4'
      || atom->name_[0] == '5') {
    return true;
  }

  return false;
}

//@{
// Given two atoms and the two residues containing them, return true if they are h-bonded.
// A quicky, non-rigorous method.
//@}
bool MIAtom::MIIsHBondable(const MIAtom* a1, const MIAtom* a2) {
  return ((a1->name_[0] == 'N' && a2->name_[0] == 'O')
          || (a1->name_[0] == 'O' && a2->name_[0] == 'N')
          || (a1->name_[0] == 'O' && a2->name_[0] == 'O') );
}

std::string MIAtom::liststring(const MIAtom* atom) {
  if (atom->altloc_ == ' ') {
    return format("%s %0.3f %0.3f %0.3f %0.2f %0.2f %d", atom->name_, atom->x_, atom->y_, atom->z_,
             atom->BValue_, atom->occ_, (int)atom->atomicnumber_);
  } else {
    return format("%s %c %0.3f %0.3f %0.3f %0.2f %0.2f %s", atom->name_, atom->altloc_, atom->x_, atom->y_, atom->z_,
             atom->BValue_, atom->occ_, Atomic_Name(atom->atomicnumber_));
  }
}

}


