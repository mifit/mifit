#ifndef SEQ_UTIL_H
#define SEQ_UTIL_H

#include <vector>

namespace chemlib
{
/////////////////////////////////////////////////////////////////////////////
// Function:    HammingDistance
// Purpose:		General utility to count the number of mismatches between two
//              vectors, ignoring the overhang if one is longer than the other
// Input:       Two vectors
// Output:      # of mismtaches
// Requires:
/////////////////////////////////////////////////////////////////////////////

    template<class T>
    int HammingDistance(std::vector<T> &v1, std::vector<T> &v2)
    {
        int dist = 0;
        for (unsigned int i = 0; i < v1.size() && i < v2.size(); ++i)
        {
            if (v1[i] != v2[i])
            {
                ++dist;
            }
        }

        return dist;
    }

    template<typename T>
    void InitializeArray(T *array, std::size_t size, T value)
    {
        for (unsigned int i = 0; i < size; ++i)
        {
            array[i] = value;
        }
    }

    template<class T>
    int GetIndex(const T query, const std::vector<T> &domain)
    {
        for (unsigned int i = 0; i < domain.size(); ++i)
        {
            if (query == domain[i])
            {
                return i;
            }
        }
        return -1;
    }

} //namespace chemlib
#endif //SEQ_UTIL
