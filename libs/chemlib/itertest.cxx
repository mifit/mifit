#include "iterator.h"
#include "stdio.h"

// test harness for creation of MIIter classes
// named cxx to avoid being put into compilation of library


class RESIDUE
{
  public:
    RESIDUE() : resnum(0),next_res(0) {}
    int resnum;
    RESIDUE *next() { return next_res; }
    RESIDUE *next_res;
};

RESIDUE *BuildResList()
{
  RESIDUE *head=new RESIDUE;
  RESIDUE *res=head;

  for (unsigned int i=0;i<10; ++i)
  {
    res->resnum=i;
    if (i==9)
      break;
    res->next_res=new RESIDUE();
    res=res->next();
  }
  return head;
}

class MIMoleculeBase
{
  public:
    MIMoleculeBase(RESIDUE *res) : _res(res) {}
    MIIterBase<RESIDUE> *GetResidues();

  private:
    RESIDUE *_res;
};


MIIterBase<RESIDUE> *MIMoleculeBase::GetResidues()
{
  return new MISinglyLinkedListIter<RESIDUE>(_res);
}

void IterateTest(MIMoleculeBase *mol)
{
  unsigned int count=0;
  for (MIIter<RESIDUE> ri=mol->GetResidues(); ri; ++ri)
  {
    printf("Residue %d has number %d\n",count++,ri->resnum);
  }

  count=9;
  MIIter<RESIDUE> ri=mol->GetResidues();
  MIIter<RESIDUE> r2=mol->GetResidues();
  ri.Last();
  for (; ri; --ri)
  {
    printf("Backwards: residue %d has number %d\n",count--,ri->resnum);
    if (ri->resnum==5)
      r2=ri;
  }

  printf("Assignment test\n");
  for (; r2; ++r2)
  {
    printf("residue number %d\n",r2->resnum);
  }
}



int main(int argc, char **argv)
{
  RESIDUE *res=BuildResList();
  MIMoleculeBase *mol=new MIMoleculeBase(res);
  IterateTest(mol);


  mol=new MIMoleculeBase(0);
  MIIter<RESIDUE> ri=mol->GetResidues();
  ri.Last();
	RESIDUE *res2 = ri;
}
