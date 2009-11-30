#if !defined (AFX_CLUSTERLIST_H_)
#define AFX_CLUSTERLIST_H_

#include <chemlib/chemlib.h>

class Cluster
{
public:
    std::vector<chemlib::MIAtom*> atoms;
};

class ClusterList
{
    std::vector<Cluster> clusters;
    float dmax;
public:
    size_t size()
    {
        return clusters.size();
    }

    float GetMaxDistance()
    {
        return dmax;
    }

    void SetMaxDistance(float d)
    {
        dmax = d;
    }

    chemlib::RESIDUE *BuildClusters(chemlib::RESIDUE *first, chemlib::RESIDUE *last = NULL, int min_size = 5);
    ClusterList();
    virtual ~ClusterList();

};

#endif // !defined(AFX_CLUSTERLIST_H__DB3E55C3_663D_4F00_987B_96284A07A52F__INCLUDED_)
