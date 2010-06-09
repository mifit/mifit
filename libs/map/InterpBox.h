#ifndef mifit_map_InterpBox_h
#define mifit_map_InterpBox_h

#include <vector>

namespace chemlib
{
    class MIAtom;
    class Residue;
}
class EMapBase;

//@{
// A fine grid preinterpolated from the map to speed up interpolations
// where the same area is repeatedly accessed.
//@}
class InterpBox
{
private:
    std::vector<float> grid_points;
    float min_x, min_y, min_z, max_x, max_y, max_z;
    int nx, ny, nz;
    float spacing;
    chemlib::MIAtom **atoms;
    int natoms;
    EMapBase *emap;
    void Init();
public:
    void ZeroModel(chemlib::Residue *model);
    InterpBox(std::vector<chemlib::MIAtom*> &atoms, EMapBase *from_map);
    float RDensity(std::vector<chemlib::MIAtom*> CurrentAtoms);
};


#endif // ifndef mifit_map_InterpBox_h
