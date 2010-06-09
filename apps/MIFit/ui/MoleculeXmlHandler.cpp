#include "MoleculeXmlHandler.h"

#include <chemlib/chemlib.h>
#include <map/maplib.h>
#include "core/corelib.h"
#include <util/utillib.h>
#include <chemlib/Residue.h>

#include "EMap.h"

using namespace chemlib;
using namespace std;

MoleculeXmlHandler::MoleculeXmlHandler()
    : residueList(NULL),
      residue(NULL),
      prev_res(NULL),
      label(NULL)
{

    in_MapHeader = false;
    in_Annotation = false;
    in_Defaults = false;
    in_Molecule = false;
    in_EMap = false;
}

MoleculeXmlHandler::~MoleculeXmlHandler()
{
}

bool MoleculeXmlHandler::characters( const QString &chars)
{
    data += chars.toStdString();
    return true;
}

bool MoleculeXmlHandler::startElement(const QString &namespaceURI,
                                      const QString &name, const QString &qName, const QXmlAttributes &attributes)
{

    Q_UNUSED(namespaceURI);
    Q_UNUSED(qName);

    char buf[1024];
    data.clear();
    if (name == "Molecule")
    {
        molecule = new Molecule(MoleculeType::XML);
        in_Molecule = true;
        unsigned l = attributes.length();
        for (unsigned int i = 0; i < l; i++)
        {
            if (attributes.localName(i) == "name")
            {
                molecule->compound = attributes.value(i).toStdString();
            }
        }
    }

    if (name == "MapHeader")
    {
        in_MapHeader = true;
        if (in_Molecule && molecule->mapheader == NULL)
        {
            CMapHeader mapHeader;
            molecule->SetMapHeader(mapHeader);
        }
    }

    if (name == "Defaults")
    {
        in_Defaults = true;
    }

    if (name == "EMap")
    {
        in_EMap = true;
    }

    if (name == "ResidueList" || name == "SymmResidueList")
    {
        MI_ASSERT(residueList == NULL);
    }

    if (name == "Residue")
    {
        MI_ASSERT(residue == NULL);
        residue = new Residue();
        residue->setSecstr('U');
        if (residueList == NULL)
        {
            residue->set_linkage_type(NTERMINUS);
            residueList = residue;
        }
        else
        {
            residue->set_linkage_type(MIDDLE);
        }
        atoms.clear();
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "name")
            {
                residue->setName(attributes.value(i).toStdString());
            }
            if (attributes.localName(i) == "type")
            {
                residue->setType(attributes.value(i).toStdString());
            }
            if (attributes.localName(i) == "chain")
            {
                residue->set_chain_id((short unsigned int)atoi( attributes.value(i).toStdString().c_str() ));
            }
            if (attributes.localName(i) == "secstr")
            {
                strcpy(buf, attributes.value(i).toStdString().c_str());
                residue->setSecstr(buf[0]);
            }
            if (attributes.localName(i) == "linkage_type")
            {
                strcpy(buf, attributes.value(i).toStdString().c_str());
                residue->set_linkage_type((short unsigned int)atoi( attributes.value(i).toStdString().c_str() ));
            }
            if (attributes.localName(i) == "confomer")
            {
                residue->setConfomer((short unsigned int)atoi( attributes.value(i).toStdString().c_str() ));
            }
            if (attributes.localName(i) == "name1")
            {
                strcpy(buf, attributes.value(i).toStdString().c_str());
                residue->setName1(buf[0]);
            }
        }
    }

    if (name == "Atom")
    {
        atom = new MIAtom;
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "name")
            {
                atom->setName(attributes.value(i).toStdString().c_str());
            }

            if (attributes.localName(i) == "x")
            {
                atom->setX(atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "y")
            {
                atom->setY(atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "z")
            {
                atom->setZ(atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "color")
            {
                atom->setColor(atoi( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "altloc")
            {
                strcpy(buf, attributes.value(i).toStdString().c_str() );
                atom->setAltloc(buf[0]);
            }

            if (attributes.localName(i) == "occupancy")
            {
                atom->setOcc(atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "BValue")
            {
                atom->setBValue(atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "number")
            {
                atom->setAtomnumber(atoi( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "atomicnumber")
            {
                atom->setAtomicnumber(atoi( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "U11")
            {
                if (!atom->hasAnisotropicity())
                {
                    atom->newAnisotropicity();
                }
                atom->U(0, atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "U22")
            {
                if (!atom->hasAnisotropicity())
                {
                    atom->newAnisotropicity();
                }
                atom->U(1, atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "U33")
            {
                if (!atom->hasAnisotropicity())
                {
                    atom->newAnisotropicity();
                }
                atom->U(2, atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "U12")
            {
                if (!atom->hasAnisotropicity())
                {
                    atom->newAnisotropicity();
                }
                atom->U(3, atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "U13")
            {
                if (!atom->hasAnisotropicity())
                {
                    atom->newAnisotropicity();
                }
                atom->U(4, atof( attributes.value(i).toStdString().c_str() ));
            }

            if (attributes.localName(i) == "U23")
            {
                if (!atom->hasAnisotropicity())
                {
                    atom->newAnisotropicity();
                }
                atom->U(5, atof( attributes.value(i).toStdString().c_str() ));
            }

        }
    }

    if (name == "Bond")
    {
        bond.setAtom1(NULL);
        bond.setAtom2(NULL);
        bond.type = 0;
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "atom1")
            {
                int natom = atoi( attributes.value(i).toStdString().c_str() );
                bond.setAtom1(molecule->GetAtom(natom));
            }
            if (attributes.localName(i) == "atom2")
            {
                int natom = atoi( attributes.value(i).toStdString().c_str() );
                bond.setAtom2(molecule->GetAtom(natom));
            }
            if (attributes.localName(i) == "type")
            {
                bond.type = atoi(attributes.value(i).toStdString().c_str()) ;
            }
        }
    }

    if (name == "H_Bond")
    {
        bond.setAtom1(NULL);
        bond.setAtom2(NULL);
        bond.type = 0;
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "atom1")
            {
                int natom = atoi( attributes.value(i).toStdString().c_str() );
                bond.setAtom1(molecule->GetAtom(natom));
            }
            if (attributes.localName(i) == "atom2")
            {
                int natom = atoi( attributes.value(i).toStdString().c_str() );
                bond.setAtom2(molecule->GetAtom(natom));
            }
        }
    }

    //<Connect atom1="665" atom2="1635" />
    if (name == "Connect")
    {
        bond.setAtom1(NULL);
        bond.setAtom2(NULL);
        bond.type = 0;
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "atom1")
            {
                int natom = atoi( attributes.value(i).toStdString().c_str() );
                bond.setAtom1(molecule->GetAtom(natom));
            }
            if (attributes.localName(i) == "atom2")
            {
                int natom = atoi( attributes.value(i).toStdString().c_str() );
                bond.setAtom2(molecule->GetAtom(natom));
            }
        }
    }
    //<Label atom="34" label="VAL 3 CA" red="255" green="255" blue="255" n="8" xo="2" yo="-2" />
    if (name == "Label")
    {
        MIAtom *atom = NULL;
        char labelString[128];
        int xo, yo;
        int visible;
        unsigned char red, green, blue;
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "atom")
            {
                int natom = atoi( attributes.value(i).toStdString().c_str() );
                atom = molecule->GetAtom(natom);
            }
            if (attributes.localName(i) == "label")
            {
                strncpy(labelString, attributes.value(i).toStdString().c_str(), 128);
            }
            if (attributes.localName(i) == "red")
            {
                red = (unsigned char)atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "green")
            {
                green = (unsigned char)atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "blue")
            {
                blue = (unsigned char)atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "visible")
            {
                visible = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "xo")
            {
                yo = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "yo")
            {
                xo = atoi( attributes.value(i).toStdString().c_str() );
            }
        }
        if (residue == NULL)
        {
            residue = residue_from_atom(molecule->getResidues(), atom);
        }
        if (residue != NULL && atom != NULL)
        {
            label = new ATOMLABEL(residue, atom);
            label->label(labelString);
            label->visible(visible != 0);
            label->xOffset(xo);
            label->yOffset(yo);
            label->red(red);
            label->green(green);
            label->blue(blue);
        }
    }
    if (name == "SymmOps" && in_MapHeader)
    {
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "spgpno")
            {
                mapheader.spgpno = atoi( attributes.value(i).toStdString().c_str() );
                mapheader.SetSymmOps();
            }
        }
    }
    if (name == "Limits" && in_MapHeader)
    {
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "hmin")
            {
                mapheader.hmin = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "hmax")
            {
                mapheader.hmax = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "kmin")
            {
                mapheader.kmin = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "kmax")
            {
                mapheader.kmax = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "lmin")
            {
                mapheader.lmin = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "lmax")
            {
                mapheader.lmax = atoi( attributes.value(i).toStdString().c_str() );
            }
        }
    }
    if (name == "MapHeaderOptions" && in_MapHeader)
    {
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "M_Coefficient")
            {
                mapheader.M_Coefficient = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "N_Coefficient")
            {
                mapheader.N_Coefficient = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "fc_is_fom")
            {
                mapheader.fc_is_fom = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "Bsolvent")
            {
                mapheader.Bsolvent = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "Ksolvent")
            {
                mapheader.Ksolvent = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "use_aniso")
            {
                mapheader.use_aniso = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "use_bulksolvent")
            {
                mapheader.use_bulksolvent = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "title")
            {
                mapheader.title = attributes.value(i).toStdString().c_str();
            }
            if (attributes.localName(i) == "crystal_name")
            {
                mapheader.crystal_name = attributes.value(i).toStdString().c_str();
            }
        }
    }
    if (name == "UnitCell" && in_MapHeader)
    {
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "a")
            {
                mapheader.a = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "b")
            {
                mapheader.b = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "c")
            {
                mapheader.c = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "alpha")
            {
                mapheader.alpha = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "beta")
            {
                mapheader.beta = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "gamma")
            {
                mapheader.gamma = atof( attributes.value(i).toStdString().c_str() );
            }
        }
    }
    if (name == "Annotation")
    {
        in_Annotation = true;
        ann = new Annotation();
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "Id")
            {
                ann->m_id = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "x")
            {
                ann->m_x = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "y")
            {
                ann->m_y = atof( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "z")
            {
                ann->m_z = atof( attributes.value(i).toStdString().c_str() );
            }
        }
    }
    if (name == "Color" && in_Annotation == true)
    {
        unsigned char Red = ann->m_color.red;
        unsigned char Green = ann->m_color.green;
        unsigned char Blue = ann->m_color.blue;
        for (int i = 0; i < attributes.length(); i++)
        {
            if (attributes.localName(i) == "red")
            {
                Red = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "green")
            {
                Green = atoi( attributes.value(i).toStdString().c_str() );
            }
            if (attributes.localName(i) == "blue")
            {
                Blue = atoi( attributes.value(i).toStdString().c_str() );
            }
        }
        ann->m_color = PaletteColor(Red, Green, Blue);
    }
    return true;
}

bool MoleculeXmlHandler::endElement(const QString &namespaceURI, const QString &name, const QString &qName)
{
    Q_UNUSED(namespaceURI);
    Q_UNUSED(qName);

    // No escapes are legal here
    if ( name == "Atom" )
    {
        atoms.push_back(atom);
    }
    if ( name == "Residue" )
    {
        // finish off residue
        residue->setAtoms(atoms);
        atoms.clear();
        // add new res to the end of the list
        if (prev_res == NULL)
        {
            prev_res = residue;
            residue = NULL;
        }
        else
        {
            prev_res = prev_res->insertResidue(residue);
            residue = NULL;
        }
    }
    if ( name == "RibbonAtoms" )
    {
        vector<MIAtom*>::iterator i, e;
        i = atoms.begin();
        e = atoms.end();
        for (; i != e; i++)
        {
            (*i)->setType(AtomType::RIBBONATOM);
            molecule->ribbonatoms.push_back(*i);
            molecule->nribbons++;
            delete *i;
        }
        atoms.clear();
    }
    if ( name == "MapHeader" )
    {
        in_MapHeader = false;
        if (in_Molecule)
        {
            molecule->SetMapHeader(mapheader);
        }
        else if (in_EMap)
        {
            mapHeaders.push_back(mapheader);
        }
    }
    if ( name == "Molecule" )
    {
        in_Molecule = false;
        molecule->InitSeqPos();
        molecules.push_back(molecule);
    }

    if ( name == "Compound" )
    {
        molecule->compound = data;
    }

    if ( name == "Defaults" )
    {
        in_Defaults = false;
    }

    if ( name == "EMap" )
    {
        in_EMap = false;
    }

    if ( name == "Annotation" )
    {
        in_Annotation = false;
        // add annotation to Molecule
        if (molecule)
        {
            molecule->addAnnotation(ann);
        }
    }

    if ( name == "Text" && in_Annotation)
    {
        ann->SetText(data.c_str());
    }

    if ( name == "Author" )
    {
        molecule->author = data;
    }

    if ( name == "Pathname" )
    {
        molecule->pathname = data;
    }

    if ( name == "Source" )
    {
        molecule->source = data;
    }

    if ( name == "dots_visible" )
    {
        molecule->dots_visible = atoi(data.c_str()) != 0;
    }

    if ( name == "HVisible" )
    {
        molecule->HVisible = atoi(data.c_str()) != 0;
    }

    if ( name == "labels_visible" )
    {
        molecule->labels_visible = atoi(data.c_str()) != 0;
    }

    if ( name == "modelnumber" )
    {
        molecule->modelnumber = atoi(data.c_str()) != 0;
    }

    if ( name == "link_here" )
    {
        strncpy(molecule->link_here, data.c_str(), MAXNAME);
    }

    if ( name == "link_next" )
    {
        strncpy(molecule->link_next, data.c_str(), MAXNAME);
    }

    if ( name == "Bond" )
    {
        if (bond.getAtom1() && bond.getAtom2())
        {
            molecule->bonds.push_back(bond);
        }
    }

    if ( name == "H_Bond" )
    {
        if (bond.getAtom1() && bond.getAtom2())
        {
            molecule->hbonds.push_back(bond);
        }
    }

    if ( name == "Connect" )
    {
        if (bond.getAtom1() && bond.getAtom2())
        {
            molecule->connects.push_back(bond);
            molecule->Connect(bond);
        }
    }

    if ( name == "Label" )
    {
        if (label != NULL)
        {
            if (label->atom() != NULL)
            {
                molecule->addAtomLabel(label);
            }
            else
            {
                delete label;
            }
        }
    }

    if (name == "ResidueList")
    {
        if (prev_res)
        {
            prev_res->set_linkage_type(CTERMINUS);
        }
        molecule->residues = residueList;
        residueList = NULL;
        prev_res = NULL;
    }
    if (name == "SymmResidueList")
    {
        if (prev_res)
        {
            prev_res->set_linkage_type(CTERMINUS);
        }
        molecule->SymmResidues = residueList;
        residueList = NULL;
        prev_res = NULL;
    }
    return true;
}

bool MoleculeXmlHandler::error(const QXmlParseException &exception)
{
    errorString_ = QObject::tr("Error in file %1 at line %2, char %3: %4")
                   .arg(exception.systemId())
                   .arg(exception.lineNumber())
                   .arg(exception.columnNumber())
                   .arg(exception.message());
    return false;
}

bool MoleculeXmlHandler::fatalError(const QXmlParseException &exception)
{
    errorString_ = QObject::tr("Fatal error in file %1 at line %2, char %3: %4")
                   .arg(exception.systemId())
                   .arg(exception.lineNumber())
                   .arg(exception.columnNumber())
                   .arg(exception.message());
    return false;
}

bool MoleculeXmlHandler::warning(const QXmlParseException &exception)
{
    errorString_ = QObject::tr("nWarning in file %1 at line %2, char %3: %4")
                   .arg(exception.systemId())
                   .arg(exception.lineNumber())
                   .arg(exception.columnNumber())
                   .arg(exception.message());
    return true;
}
