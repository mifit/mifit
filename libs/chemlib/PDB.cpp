#include <stdio.h>

#include <nongui/nonguilib.h>
#include <math/mathlib.h>


#include <chemlib/Monomer.h>
#include "MIMolIOBase.h"
#include "PDB.h"
#include "mol_util.h"
#include "atom_util.h"

using namespace std;

namespace chemlib
{

PDB::PDB()
{
}

PDB::~PDB()
{
}

bool PDB::Write(FILE *fp, MIMolInfo &mol)
{
    //TODO: obviously glue in the MIFit calls for writing a pdb
    bool ret = false;

    if (!fprintf(fp, "COMPND    %3s\n", mol.res->type().c_str()))                   //Molecule/residue name
    {
        return false;
    }

    ret = SavePDB(fp, mol.res, ((mol.bonds.size() == 0) ? NULL : &mol.bonds[0]), mol.bonds.size(), true, &mol.header, &mol.tail);
    return ret;
}

//NOTE: this only returns a single residue from a PDB file
bool PDB::Read(FILE *fp, MIMolInfo &mol)
{
    vector<Bond> connects;

    // clear current mol
    MIMolInfo foo;
    mol = foo;
    mol.res = 0;

    Residue *tmp_res = LoadPDB(fp, &connects);
    if (!Monomer::isValid(tmp_res))
    {
        return false;
    }

    MIAtom *atom1, *atom2;
    vector<Bond> tmp_bonds;
    float dlimit, d;
    int i, j;
    for (i = 0; i < tmp_res->atomCount(); ++i)
    {
        tmp_res->atom(i)->setCharge(0.0F);
        for (j = i+1; j < tmp_res->atomCount(); ++j)
        {
            atom1 = tmp_res->atom(i);
            atom2 = tmp_res->atom(j);

            // check to see if distance bondable
            dlimit = (float)(BondLimit(atom1->name())+ BondLimit(atom2->name()));
            if ((d = (float)AtomDist(*atom1, *atom2)) < dlimit
                || chemlib::AlreadyBonded(atom1, atom2, connects))
            {
                CreateBond(atom1, atom2, tmp_bonds);
            }
        }
    }

    tmp_res->setPrefBonds(tmp_bonds);
    mol.res = new Residue(*tmp_res);
    mol.bonds = mol.res->prefBonds();
    mol.res->clearPrefBonds();
    mol.res->clearPrefAngles();
    FreeResidueList(tmp_res);
    return true;
}

void PDB::CreateBond(MIAtom *atom1, MIAtom *atom2, vector<Bond> &bonds)
{
    Bond bond;
    bond.setAtom1(atom1);
    bond.setAtom2(atom2);
    bond.getAtom1()->addType(AtomType::BONDED);
    bond.getAtom2()->addType(AtomType::BONDED);
    bond.type = B_NORMAL;
    bonds.push_back(bond);
}

#define strkeq(s1, s2, k) /* true if first k chars of string s1 and s2 match */ \
    ( k <= 0 ? 1 : s1[0] == s2[0] ? !strncmp(s1, s2, (k)) : 0)

static bool IGNORE_DUMMY = true;
void MISetIgnoreDummyAtomsOnLoad(bool ignore)
{
    IGNORE_DUMMY = ignore;
}

Residue *LoadPDB(FILE *f, std::vector<Bond> *connects)
{
#define column(n) (linebuf+(n)-1)
    char linebuf[200];
    double x, y, z;
    double Bvalue;
    double occupancy;
    char abuf[MAXATOMNAME];
    char rname[MAXNAME];
    char rtype[MAXNAME];
    char oldrtype[MAXNAME];
    char oldrname[MAXNAME];
    char chainid;
    int bonding_group, oldgroup = 0;
    int ter_count; /* used to form bonding_group value:
                      low byte is chain ID from PDB file,
                      next higher byte is ter_count */
    int is_atom, is_hetatm;
    long na = 0, n = 0, i, ii;
    int U11, U22, U33, U12, U13, U23;
    Residue *res, *res1;
    int j, ia[11], lb;
    static int modelno = -1;
    MIAtom *atom;
    MIAtomList atoms;
    na = 0;
    rewind(f);
    atoms.clear();
    connects->clear();

    if (NULL == (res = new Residue()))
    {
        return (NULL);
    }
    res1 = res;
    res->set_linkage_type(NTERMINUS);
    res->setSecstr('U');

    ter_count = 0;
    while (fgets(linebuf, sizeof linebuf, f) != NULL)
    {
        is_atom = strkeq("ATOM", linebuf, 4) || strkeq("atom", linebuf, 4);
        is_hetatm = is_atom ? 0 : strkeq("HETATM", linebuf, 6)
                    || strkeq("hetatm", linebuf, 6) ;
        if ( (is_atom || is_hetatm) )
        {
            atom = new MIAtom;
            lb = strlen(linebuf);
            if (lb < 53)
            {
                static int ntimes = 0;
                if (ntimes < 10)
                {
                    std::string m;
                    m = "ATOM record too short - not long enough to contain x,y,z info - ignored\n:";
                    m += linebuf;
                } //else if(ntimes == 10) Logger::message("Length message will not be repeated - more than 10 occurences");
                ntimes++;
                continue;
            }
            strncpy(abuf, column(13), 4);
            abuf[4] = '\0';
            i = 3;
            while (abuf[i] == ' ' && i > 0)
            {
                abuf[i] = '\0';
                i--;
            }
            ii = 0;
            while (i != 0 && abuf[ii] == ' ' && ii <= 4)
            {
                ii++;
            }
            atom->setName(abuf+ii);

            //Give the atom a color
            if (MIGetColorSetter())
            {
                (*MIGetColorSetter())(atom);
            }

            x  = atof(column(30));
            y  = atof(column(39));
            z  = atof(column(47));
            // catch nonsense atoms sometimes added by XPLOR
            if (x >= 9999.0 || y >= 9999.0 || z >= 9999.0)
            {
                if (IGNORE_DUMMY)
                {
                    continue;
                }
            }
            atom->setPosition((float)x, (float)y, (float)z);
            if (lb > 55)
            {
                occupancy = atof(column(55));
            }
            else
            {
                occupancy = 1.0;
            }
            atom->setOcc((float)occupancy);
            if (lb > 61)
            {
                Bvalue = atof(column(61));
            }
            else
            {
                Bvalue = 15.0;
            }
            if (Bvalue < 0.0)
            {
                Bvalue = 15.0;
            }
            atom->setBValue((float)Bvalue);
            atom->set_radius_type(0);
            atom->setType(MIAtom::MIGetAtomTypeFromName(atom->name()));
            atom->setMass(atoi(column(7))); // temporary storage for atomnumber
            // used for figuring CONECT info later on
            chainid = *column(22);
            atom->setAltloc(linebuf[16]);
            atom->setAtomicnumber(Atomic_Number(column(13)));
            atom->setAtomnumber(n + 1);
            atom->setSymmop(0);
            atom->setCharge(0);
            atom->set_formal_charge(0);
            atoms.push_back(atom);

            bonding_group = (ter_count<<8)|chainid;
            if (linebuf[17] == ' ' || linebuf[17] == '\''
                || (isalnum(linebuf[17]) && isalnum(linebuf[20])
                    && isalnum(linebuf[16]) && isalnum(linebuf[15])))
            {
                /* the ' is a special case from hyram leffert */
                sscanf(column(19), "%3s", rtype);
            }
            else
            {
                sscanf(column(18), "%3s", rtype);
            }
            sscanf(column(23), "%4s", rname);
            n++;
            na++;
            if (n == 1) // initialize firt residue.
            {
                strcpy(oldrname, rname);
                strcpy(oldrtype, rtype);
                oldgroup = bonding_group;
                res->setName(rname);
                res->setType(rtype);
                res->set_chain_id(bonding_group);
            }
            else if (bonding_group != oldgroup
                     || strcmp(oldrname, rname)
                     || strcmp(oldrtype, rtype) )
            {
                /* new residue */
                res->clearAtoms();
                for (j = 0; j < na-1; j++)
                {
                    res->addAtom(atoms[n-na+j]);
                }
                res->set_linkage_type(MIDDLE);
                res->setSecstr('U');
                res = res->insertResidue(new Residue());
                res->setName(rname);
                res->setType(rtype);
                res->set_chain_id(bonding_group);
                strcpy(oldrname, rname);
                strcpy(oldrtype, rtype);
                oldgroup = bonding_group;
                na = 1;
            }
        }
        else if (strkeq("TER", linebuf, 3) || strkeq("ter", linebuf, 3))
        {
            ter_count++;
        }
        else if (strkeq("ANISOU", linebuf, 6))
        {
            if (atom)
            {
                int an;
                sscanf(linebuf, "%*s%d", &an);
                if (an == atom->mass())
                {
                    atom->newAnisotropicity();
                    sscanf(column(29), "%7d%7d%7d%7d%7d%7d", &U11, &U22, &U33, &U12, &U13, &U23);
                    atom->U(0, U11/10000.0f);
                    atom->U(1, U22/10000.0f);
                    atom->U(2, U33/10000.0f);
                    atom->U(3, U12/10000.0f);
                    atom->U(4, U13/10000.0f);
                    atom->U(5, U23/10000.0f);
                }
                else
                {
                    //emess("ANISOU out of order");
                    //emess(linebuf);
                }
            }
        }
    }
    // if no residues, return NULL
    if (n == 0)
    {
        Logger::message("Error: No atoms found in file!");
        return NULL;
    }
    // last residue
    res->clearAtoms();
    for (j = 0; j < na; j++)
    {
        res->addAtom(atoms[n-na+j]);
    }
    res->setSecstr('U');
    res->set_chain_id(bonding_group);
    res->set_linkage_type(CTERMINUS);

    //  find the CONECT records and store them
    int nid, k, found;
    char fbuf[10];
    rewind(f);
    long ja;
    MIAtom *a1 = NULL;
    Bond connect;
    while (fgets(linebuf, sizeof linebuf, f) != NULL)
    {
        if (strkeq("CONECT", linebuf, 6))
        {
            nid = 1;
            fbuf[5] = '\0';
            for (j = 0; j < 5; j++)
            {
                fbuf[j] = linebuf[j+6];
            }
            if (sscanf(fbuf, "%d", &ia[0]) == 1)
            {
                found = 0;
                for (ja = 0; ja < n; ja++)
                {
                    if (atoms[ja]->mass() == ia[0])
                    {
                        found = 1;
                        a1 = atoms[ja];
                        break;
                    }
                }
                if (!found)
                {
                    continue; //go back to searching for CONECT's
                }
                for (i = 1; i < 11; i++)
                {
                    j = 5*i+6;
                    for (k = 0; k < 5; k++)
                    {
                        fbuf[k] = linebuf[k+j];
                    }
                    if (sscanf(fbuf, "%d", &ia[i]) == 1)
                    {
                        nid++;
                    }
                    else
                    {
                        break;
                    }
                    // search for these atoms numbers in atom list
                    for (ja = 0; ja < n; ja++)
                    {
                        if (atoms[ja]->mass() == ia[i])
                        {
                            // found it! make connect
                            connect.setAtom1(a1);
                            connect.setAtom2(atoms[ja]);
                            connects->push_back(connect);
                            break;
                        }
                    }
                }
            }
        }
    }
    modelno++;
    return (res1);
}

bool SavePDB(FILE *fp, Residue *res, Bond *Connects, int nConnects, bool mark_end, std::vector<std::string> *head, std::vector<std::string> *tail)
{
    /* if nconnects negative then save only visible atoms
     */
    char buf[200];
    char chainid, recid[7];
    int n = 0, nchain = 0;
    Residue *residue;
    //float occupancy = 1.0;
    //float Bvalue = 15.0;
    if (head != NULL)
    {
        if (head->size() > 0)
        {
            for (unsigned int i = 0; i < head->size(); i++)
            {
                fprintf(fp, "%s\n", (*head)[i].c_str());
            }
        }
    }
    /* now write out the atoms... */
    residue = res;
    char residueName[6];
    residueName[5] = '\0';
    while (residue != NULL)
    {
        strncpy(residueName, residue->name().c_str(), 4);
        if (strlen(residueName) == 0)
        {
            residueName[0] = '1';
            residueName[1] = '\0';
        }
        else if (strlen(residueName) < 5 && !isalpha(residue->name()[residue->name().size()-1]))
        {
            residueName[strlen(residueName)+1] = '\0';
            residueName[strlen(residueName)] = ' ';
        }
        //sprintf(buf,"Writing out %s_%s\n",residue->name,residue->type);
        //emess(buf);
        chainid = (char)(residue->chain_id() & 255);
        if (IsPeptide(*residue) || IsNucleic(residue))
        {
            strcpy(recid, "ATOM  ");
        }
        else
        {
            strcpy(recid, "HETATM");
        }
        for (int i = 0; i < residue->atomCount(); i++)
        {
            if (nConnects == (-1) && residue->atom(i)->isHidden())
            {
                continue;
            }
            n++;
            nchain++;
            const char *atomicname = Atomic_Name(residue->atom(i)->atomicnumber());
            if (strlen(residue->atom(i)->name()) < 4 && !(isdigit(residue->atom(i)->name()[0]))
                && atomicname[0] == ' ')
            {
                sprintf(buf,
                        "%6s%5d  %-3s%c%-4s%c%5s%11.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                        recid, n,
                        residue->atom(i)->name(), residue->atom(i)->altloc(), residue->type().c_str(),
                        chainid, residueName, (float)residue->atom(i)->x(),
                        (float)residue->atom(i)->y(), (float)residue->atom(i)->z(),
                        (float)residue->atom(i)->occ(), (float)residue->atom(i)->BValue(), atomicname);
            }
            else if (strlen(residue->atom(i)->name()) < 5)
            {
                sprintf(buf,
                        "%6s%5d %-4s%c%-4s%c%5s%11.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                        recid, n,
                        residue->atom(i)->name(), residue->atom(i)->altloc(), residue->type().c_str(),
                        chainid, residueName, (float)residue->atom(i)->x(),
                        (float)residue->atom(i)->y(), (float)residue->atom(i)->z(),
                        (float)residue->atom(i)->occ(), (float)residue->atom(i)->BValue(), atomicname);
            }
            else
            {
                sprintf(buf,
                        "%6s%5d%-5s%c%-4s%c%5s%11.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                        recid, n,
                        residue->atom(i)->name(), residue->atom(i)->altloc(), residue->type().c_str(),
                        chainid, residueName, (float)residue->atom(i)->x(),
                        (float)residue->atom(i)->y(), (float)residue->atom(i)->z(),
                        (float)residue->atom(i)->occ(), (float)residue->atom(i)->BValue(), atomicname);
            }
            fwrite(buf, sizeof(char), strlen(buf), fp);
            // temporarily save the atom number in mass for CONECTs
            residue->atom(i)->setMass(n);
            if (residue->atom(i)->hasAnisotropicity())
            {
                MIAtom *a = residue->atom(i);
                fprintf(fp, "ANISOU");
                fprintf(fp, "%5d", a->mass());
                if (strlen(residue->atom(i)->name()) < 4 && !(isdigit(residue->atom(i)->name()[0]))
                    && atomicname[0] == ' ')
                {
                    fprintf(fp, "  %-3s%c", a->name(), a->altloc());
                }
                else if (strlen(residue->atom(i)->name()) < 5)
                {
                    fprintf(fp, " %-4s%c", a->name(), a->altloc());
                }
                else
                {
                    fprintf(fp, "%-4s%c", a->name(), a->altloc());
                }
                fprintf(fp, "%-4s%c", residue->type().c_str(), chainid);
                fprintf(fp, "%4s", residueName);
                fprintf(fp, " %7d%7d%7d%7d%7d%7d       %2s\n",
                        ROUND(a->U(0)*10000.0),
                        ROUND(a->U(1)*10000.0),
                        ROUND(a->U(2)*10000.0),
                        ROUND(a->U(3)*10000.0),
                        ROUND(a->U(4)*10000.0),
                        ROUND(a->U(5)*10000.0),
                        atomicname);
            }
        }
        if (residue->linkage_type() &CTERMINUS)
        {
            n++;
            sprintf(buf, "TER   %5d      %-4s%c%5s\n", n, residue->type().c_str(), chainid, residueName);
            fwrite(buf, sizeof(char), strlen(buf), fp);
            //TER    1629      ARG A 106
        }

        residue = residue->next();
    }

    if (nConnects > 0 && Connects != NULL)
    {
        for (int i = 0; i < nConnects; i++)
        {
            sprintf(buf, "CONECT%5d%5d\n", Connects[i].getAtom1()->mass(),
                    Connects[i].getAtom2()->mass());
            fwrite(buf, sizeof(char), strlen(buf), fp);
        }
    }

    /* now copy everything else out of input pdb */
    if (tail != NULL)
    {
        if (tail->size() > 0)
        {
            for (unsigned int i = 0; i < tail->size(); i++)
            {
                fprintf(fp, "%s\n", (*tail)[i].c_str());
                if ((*tail)[i].compare("END") == 0)
                {
                    mark_end = false;
                }
            }
        }
    }
    if (mark_end)
    {
        sprintf(buf, "END\n");
        fwrite(buf, sizeof(char), strlen(buf), fp);
    }
    return true;
}

static MIColorSetter *COLOR_SETTER = 0;
void MIRegisterColorSetter(MIColorSetter *cs)
{
    COLOR_SETTER = cs;
}

MIColorSetter *MIGetColorSetter()
{
    return COLOR_SETTER;
}

} //namespace chemlib
