#ifndef mifit_model_SaveAtom_h
#define mifit_model_SaveAtom_h

namespace chemlib {

class MIAtom;

/**
 * Saves an atom for undo'ing later.
 */
class SaveAtom {
public:
  friend bool operator ==(SaveAtom a1, SaveAtom a2);
  float x, y, z;
  int color;
  MIAtom* from;

  SaveAtom();
  SaveAtom(MIAtom* a);
  int Restore() const;
  void RestoreColor(unsigned int mask) const;
};

}

#endif
