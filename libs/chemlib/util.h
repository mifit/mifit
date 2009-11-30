#ifndef UTIL_H
#define UTIL_H

#include <functional>
#include <vector>
#include <map>
#include <stdlib.h>

//FIXME: clean up usage of fns in this file and make as private as possible


namespace chemlib
{
/////////////////////////////////////////////////////////////////////////////
// Function:    void_mem_fun_ref (with functor class void_mem_fun_ref_t)
// Purpose:		Produces a function object that can be used with void member
//				functions to pass a pointer to a generic algorithm in an idiom
//				such as: foreach(v.begin(), v.end(), void_mem_fun_ref(&Widget:method())
// Input:       Takes a pointer to a member function
// Output:      A functor class that holds the function pointer and offers an operator()
//				that invokes the pointed-to member function on an object.
//				(Meyers, Effective STL, pg 175)
// Note:		This functionality should be provided by the mem_fun_ref adapter
//				in the STL.  Alas, in Visual C++ this does not seem to work for
//				void functions.  Hence, this new version.
/////////////////////////////////////////////////////////////////////////////
    template<class T>
    class void_mem_fun_ref_t : public std::unary_function<T, void>
    {
    public:
        explicit void_mem_fun_ref_t(void(T::*funct)()) : ptr(funct)
        {
        }

        void operator()(T &input_obj) const
        {
            (input_obj.*ptr)();
        }

    private:
        void (T::*ptr)();
    };


    template<class T> inline
    void_mem_fun_ref_t<T> void_mem_fun_ref(void (T::*_funct)())
    {
        return (void_mem_fun_ref_t<T>(_funct));
    }

/////////////////////////////////////////////////////////////////////////////
// Function:    void_mem_fun (with functor class void_mem_fun_t)
// Purpose:		Produces a function object that can be used with void member
//				functions to pass a pointer to a generic algorithm in an idiom
//				such as: foreach(v.begin(), v.end(), void_mem_fun_ref(&Widget:method())
// Input:       Takes a pointer to a member function
// Output:      A functor class that holds the function pointer and offers an operator()
//				that invokes the pointed-to member function on an object.
//				(Meyers, Effective STL, pg 175)
// Note:		This functionality should be provided by the mem_fun adapter
//				in the STL.  Alas, in Visual C++ this does not seem to work for
//				void functions.  Hence, this new version.
/////////////////////////////////////////////////////////////////////////////
    template<class T>
    class void_mem_fun_t : public std::unary_function<T*, void>
    {
    public:
        explicit void_mem_fun_t(void(T::*funct)()) : ptr(funct)
        {
        }

        void operator()(T *input_obj_ptr) const
        {
            (input_obj_ptr->*ptr)();
        }

    private:
        void (T::*ptr)();
    };

    template<class T> inline
    void_mem_fun_t<T> void_mem_fun(void (T::*funct)())
    {
        return (void_mem_fun_t<T>(funct));
    }

//From Joe Gottman, posted to the Boost mailing list
//http://lists.boost.org/MailArchives/boost/msg21093.php
//Name changed to the suggestion posted by Jeff Garland
//All the comments are mine. -KWB

/////////////////////////////////////////////////////////////////////////////
// Function:    Has_Intersection
// Purpose:		Determine if two sorted ranges have at least one common element
// Input:       Beginning and ending iterators for each range
// Output:		True if the ranges intersect
// Requires:	1. That the ranges be sorted beforehand...otherwise it will often return
//				   false when the ranges do actually intersect
//				2. That the intersection is defined by equivalence, not equality.  I.e.
//				   the operator '==' is never applied, only inferred when the '<' operator
//				   returns false in both directions.  This is consistent with the STL sorting
//				   algorithms and with STL sets & maps.  A version based on equality could
//				   be constructed using the STL find algorithm, but would not be as efficient.
//				   See Scott Meyers "Effective STL", page 83 for a discussion of equivalence
//				   and equality.
/////////////////////////////////////////////////////////////////////////////
    template <typename Iterator1, typename Iterator2>
    bool Has_Intersection(Iterator1 start1, Iterator1 end1, Iterator2 start2, Iterator2 end2)
    {
        while ((start1 != end1) && (start2 != end2))
        {
            if (*start1 < *start2)                  //Walk the iterator that points to the lower
            {
                ++start1;                           //value up one notch
            }
            else if (*start2 < *start1)
            {
                ++start2;
            }
            else
            {
                return true;        //We've found 2 equivalent elements!
            }
        }
        return false;       //If we reach here, there are no equivalent elements
    }

/////////////////////////////////////////////////////////////////////////////
// Function:    AppendFromMap
// Purpose:		Copies all data from a map, and appends to a vector
// Input:       References to the vector and map
// Output:		None
// Requires:	That we don't need to copy the keys, just the values
//				(Also, there's probably some clever way to do this directly with
//				a combination of STL binders and generic algorithms)
/////////////////////////////////////////////////////////////////////////////
/*template <typename Key, typename T>
    void AppendFromMap(vector<T> &v, const map<Key, T> &m) {

    if (v.capacity() - v.size() < m.size()) {			//Expand the vector to accomodate
        v.reserve(v.size() + m.size());					//the new elements
    }

    map< Key, T >::const_iterator i, e = m.end();
    i = m.begin();
    while (i != e) {
        v.push_back(i->second);							//Copy the value (but not the key)
 ++i;
    }
   }*/

/////////////////////////////////////////////////////////////////////////////
// Function:    AppendFromPairs
// Purpose:		Copies the second value of each pair in a range and
//				appends those values to a vector
// Input:       A references to the vector, and iterators defining the range
// Output:		None
// Requires:
/////////////////////////////////////////////////////////////////////////////
    template <typename T, typename Iterator>
    void AppendFromPairs(std::vector<T> &target, Iterator start, Iterator end)
    {
        while (start != end)
        {
            target.push_back(start->second);
            ++start;
        }
    }

/////////////////////////////////////////////////////////////////////////////
// Function:    CopyKeys
// Purpose:		Copies all the keys in a map to a new map, initializing the values to T() (e.g 0)
// Input:       References to the old and new maps
// Output:		None
// Requires:
/////////////////////////////////////////////////////////////////////////////
/*template <typename Key, typename T>
    void CopyKeys(map<Key, T> &mNew, const map<Key, T> &mOld) {

    map<Key, T>::const_iterator i, e = mOld.end();
    i = mOld.begin();
    while (i != e) {
   //			mNew[i->first];							//Create an entry for each key with the value T()
        mNew.insert(i->first, T());
 ++i;
    }
   }*/

/////////////////////////////////////////////////////////////////////////////
// Function:    irand_approx
// Purpose:		(Quickly!) return a random, non-negative integer less than the
//				input value
// Input:       Size of the range from which to pick
// Output:		"Random" integer
// Requires:	That the user doesn't mind a little bias! If RAND_MAX isn't divisible
//				by the range, lower numbers will come up a bit more often.
/////////////////////////////////////////////////////////////////////////////
    inline size_t irand_approx(int range)
    {
        if (range < 1)
        {
            return 0;
        }
        else
        {
            return rand() % range;
        }
    }

} //namespace chemlib

#endif //UTIL_H
