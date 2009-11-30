#include <util/utillib.h>
#include <chemlib/RESIDUE_.h>

#include "clusterList.h"

using namespace chemlib;
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ClusterList::ClusterList()
{
    dmax = 2.0;
}

ClusterList::~ClusterList()
{

}

RESIDUE*ClusterList::BuildClusters(RESIDUE *first, RESIDUE *last, int min_size)
{
    vector<MIAtom*> atom_list;
    RESIDUE *res = first;
    MIAtom *a1, *a2;
    unsigned int i, j;
    int found, more;
    Cluster cluster;
    if (last != NULL)
    {
        last = last->next();
    }
    while (res != NULL && res != last)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            atom_list.push_back(res->atom(i));
        }
        res = res->next();
    }
    if (atom_list.size() == 0)
    {
        return 0;
    }

    do
    {
        a1 = atom_list[0];
        cluster.atoms.push_back(a1);
        atom_list.erase(atom_list.begin());
        do
        {
            i = 0;
            found = 0;
            while (i < atom_list.size())
            {
                a2 = atom_list[i];
                more = 0;
                for (j = 0; j < cluster.atoms.size(); j++)
                {
                    a1 = cluster.atoms[j];
                    if (AtomDist(*a1, *a2) < dmax)
                    {
                        cluster.atoms.push_back(a2);
                        atom_list.erase(atom_list.begin()+i);
                        found++;
                        more++;
                        break;
                    }
                }
                if (more > 0)
                {
                    i = 0;
                }
                else
                {
                    i++;
                }
            }
        } while (found > 0);
        if ((int)cluster.atoms.size() > min_size)
        {
            clusters.push_back(cluster);
        }
        cluster.atoms.clear();
    } while (atom_list.size() > 0);

    RESIDUE *reslist = NULL, *newres;
    int color = 1;
    for (i = 0; i < clusters.size(); i++)
    {
        newres = new RESIDUE();
        newres->setName(format("%d", color));
        newres->setType("CLUST");
        newres->set_chain_id('x');
        newres->setSecstr('X');
        for (j = 0; j < clusters[i].atoms.size(); j++)
        {
            MIAtom *atom = new MIAtom;
            atom->setPosition(clusters[i].atoms[j]->x(),
                              clusters[i].atoms[j]->y(),
                              clusters[i].atoms[j]->z());
            atom->setBValue(15.0);
            atom->setOcc(1.0);
            atom->setName(format("C%d", j+1).c_str());
            atom->setColor(color);
            atom->setAtomicnumber(6);
            atom->setType(MIAtom::MIGetAtomTypeFromName(atom->name()));
            newres->addAtom(atom);
        }
        if (!reslist)
        {
            reslist = newres;
            res = reslist;
        }
        else
        {
            res = res->insertResidue(newres);
        }
        color++;
    }
    return reslist;
}

