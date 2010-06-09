#ifndef BOX_CLASS_H
#define BOX_CLASS_H

#include "MIAtom_fwd.h"
#include "Bond.h"
#include "Residue.h"

#include <vector>

namespace chemlib
{
    class Box
    {
    public:
#ifdef USE_VECTOR_RESIDUE
        Box(const std::vector<Residue>&, float);
#endif
        Box(const std::vector<Monomer*>&, float);
        bool Contains(const MIAtom&) const;
        void Expand(const MIAtom&, float margin);
        std::vector<float> Origin();
        std::vector<float> Dimensions();
        //		void Print();
    private:
        float bounds[3][2];             //Rows are X, Y, Z; Columns are min and max
    };

} //namespace chemlib
#endif //BOX_CLASS_H
