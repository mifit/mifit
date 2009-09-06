#include "AtomTyper.h"
#include "LigandPerceiver.h"

using namespace std;
namespace chemlib {

AtomTyper::AtomTyper(const RESIDUE& res, const vector<Bond>& bonds) : m_mol(res, bonds) {
  LigandPerceiver lp;

  if (m_mol.GetNumAtoms() > 0) {
    lp.AssignHybridization(&m_mol);         //Probably should be set already!
    m_mol.FindRingSystems();
    lp.AssignImpHydrogens(&m_mol);
  }
}

}
