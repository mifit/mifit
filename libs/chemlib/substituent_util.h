#ifndef SUBSTITUENT_UTIL_H
#define SUBSTITUENT_UTIL_H

#include "Substituent.h"
#include <functional>

namespace chemlib
{
    void SortSubsCounterClockwise(std::vector<Substituent> &subs);

    void EqualizeSpacing(std::vector<Substituent> &subs);

    struct LeastTheta
        : public std::binary_function<const Substituent&,
                                      const Substituent&,
                                      bool>
    {
        bool operator()(const Substituent &c1, const Substituent &c2) const
        {
            return c1._theta < c2._theta;
        }

    };
} //namespace chemlib

#endif //SUBSTITUENT_UTIL_H
