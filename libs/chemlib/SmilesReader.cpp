#include "SmilesReader.h"
#include <chemlib/Monomer.h>

//#include <ctype.h>
//#include <algorithm>


namespace chemlib
{

Monomer *SmilesToMol(const std::string &smiles,
                     std::vector<Bond> &bonds,
                     std::string &error_message)
{

    if (!ValidateSmiles(smiles, error_message))
    {
        return 0;
    }

    Ligand lig(Ligand::Smiles);
    LigandPerceiver lp;
    SmilesReader sr(smiles.c_str(), &lig);

    sr.TraverseString();

    if (sr.CountAtoms() > 0)
    {
        lp.AssignHybridization(&lig);
        lp.AssignAtomGeom(&lig);
        lp.AssignChirality(&lig);

        lig.FindRingSystems();
    }
    else
    {
        error_message = "No atoms read from Smiles string!";
        return 0;
    }

    CovalentGeometry cg(&lig, lig.residues.front());
    cg.AssignResidue();


    std::vector<Monomer*> residues;
    lig.Export(residues, bonds);
    return residues.front();
}

bool SmilesToLig(const std::string &smiles,
                 Ligand **lig,
                 std::string &error_message)
{
    if (!ValidateSmiles(smiles, error_message))
    {
        return false;
    }

    *lig = new Ligand(Ligand::Smiles);
    LigandPerceiver lp;

    SmilesReader sr(smiles.c_str(), *lig);
    sr.TraverseString();

    if (sr.CountAtoms() > 0)
    {
        lp.AssignHybridization(*lig);
        lp.AssignAtomGeom(*lig);
        lp.AssignChirality(*lig);

        (*lig)->FindRingSystems();

        CovalentGeometry cg(*lig, (*lig)->residues.front());
        cg.AssignResidue();
        return true;
    }
    else
    {
        error_message = "No atoms read from Smiles string!";
        return false;
    }
}

bool ValidateSmiles(const std::string &smi, std::string &smi_err)
{
    int i;
    int bracket_count;
    int rc_count[100];
    int rc;
    bool at_bond = false;

    bracket_count = 0;
    for (i = 0; i < 100; ++i)
    {
        rc_count[i] = 0;
    }

    std::string::const_iterator cur_pos;

    for (cur_pos = smi.begin(); cur_pos != smi.end(); ++cur_pos)
    {

        if (isdigit(*cur_pos) && bracket_count == 0)
        {
            rc = *cur_pos - '0';
            rc_count[rc]++;
        }

        if (*cur_pos == '%' && bracket_count == 0)
        {
            ++cur_pos;
            rc = (*cur_pos - '0') * 10;
            ++cur_pos;
            rc += *cur_pos - '0';
            rc_count[rc]++;
        }


        if (*cur_pos == '[')
        {
            if (bracket_count > 0)
            {
                smi_err = "Unexpected \'[\' character. Brackets cannot be nested.\n";
                return false;
            }
            ++bracket_count;
        }

        if (*cur_pos == ']')
        {
            if (bracket_count < 1)
            {
                smi_err = "Unexpected \']\' character\n";
                return false;
            }
            --bracket_count;
        }

        if (*cur_pos == '.' && at_bond)
        {
            smi_err = "Disconnect symbol cannot follow bond symbol.\n";
            return false;
        }

        if (*cur_pos == '(' && at_bond)
        {
            smi_err = "Parenthesis cannot follow bond symbol.\n";
            return false;
        }

        if (*cur_pos == ')' && at_bond)
        {
            smi_err = "Parenthesis cannot follow bond symbol.\n";
            return false;
        }

        if ((bracket_count == 0
             && *cur_pos == '-')
            || *cur_pos == '='
            || *cur_pos == '#'
            || *cur_pos == ':'
            || *cur_pos == '/'
            || *cur_pos == '\\')
        {
            if (at_bond)
            {
                smi_err = "Two consecutive bond symbols\n";
                return false;
            }
            at_bond = true;
        }
        else
        {
            at_bond = false;
        }
    }


    if (at_bond)
    {
        smi_err = "Smiles strings cannot end with a bond symbol.\n";
        return false;
    }

    if (bracket_count != 0)
    {
        smi_err = "Unmatched square bracket.\n";
        return false;
    }

    for (i = 0; i < 100; ++i)
    {
        if (rc_count[i] % 2 != 0)
        {
            char buf[1024];
            sprintf(buf, "Ring number %d is opened but not closed", i);
            smi_err = buf;
            return false;
        }
    }

    return true;
}

//SMILES string and target mol can be passed to the Constructor
SmilesReader::SmilesReader(const char *smistr, Ligand *lig)
{

    _smistr = new char[strlen(smistr)+1];
    strcpy(_smistr, smistr);

    _lig = lig;


    // Initialize variables
    _cur_pos = _smistr;
    _bond_order = SINGLEBOND;
    _bond_stereo = 0;

    _natoms = _lig->GetNumAtoms();

    char resnum[MAXNAME];
    sprintf(resnum, "%4d", (unsigned int)(_lig->residues.size()+1));
    _res = _lig->AddRes("LIG", resnum);

    _ring_closer._reader = this;

}

//SMILES string and target mol can be passed to the Constructor
SmilesReader::SmilesReader(const char *smistr, Ligand *lig, Residue &res, std::vector<Bond> &bonds)
{

    _smistr = new char[strlen(smistr)+1];
    strcpy(_smistr, smistr);

    _lig = lig;


    _res = _lig->AddRes(res, bonds);

    // Initialize variables
    _cur_pos = _smistr;
    _bond_order = SINGLEBOND;
    _bond_stereo = 0;

    _natoms = _lig->GetNumAtoms();

    _ring_closer._reader = this;
}

//The Destructor
SmilesReader::~SmilesReader()
{
    delete[] _smistr;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    CountAtoms
// Input:       A Smiles string
// Output:      The number of _atoms in the Ligand
// Requires:
/////////////////////////////////////////////////////////////////////////////

int SmilesReader::CountAtoms()
{
    char *ptr;
    int _natoms = 0;
    ptr = _smistr;
    while (*ptr)
    {
        switch (*ptr)
        {
        case 'C': _natoms++;
            ptr++;
            break;
        case 'N': _natoms++;
            ptr++;
            break;
        case 'O': _natoms++;
            ptr++;
            break;
        case 'P': _natoms++;
            ptr++;
            break;
        case 'S': _natoms++;
            ptr++;
            break;
        case 'F': _natoms++;
            ptr++;
            break;
        case 'I': _natoms++;
            ptr++;
            break;
        case 'B': _natoms++;
            ptr++;
            break;
        case 'H': _natoms++;
            ptr++;
            break;
        case 'D': _natoms++;
            ptr++;
            break;                        //Deuterium
        case 'T': _natoms++;
            ptr++;
            break;                        //Tritium

        case 'c': _natoms++;
            ptr++;
            break;
        case 'n': _natoms++;
            ptr++;
            break;
        case 'o': _natoms++;
            ptr++;
            break;
        case 'p': _natoms++;
            ptr++;
            break;
        case 's': _natoms++;
            ptr++;
            break;

        case '[': _natoms++;
            while (*ptr != ']')
            {
                ptr++;
            }
            ptr++;
            break;
        default: ptr++;
        }
    }
    return _natoms;
}

void SmilesReader::TraverseString()
{
    _res->reserveAtoms(CountAtoms());

    while (*_cur_pos)
    {
        if (isdigit(*_cur_pos))             //Digits are ring opens or closes
        {
            _ring_closer.ProcessSimple(*_cur_pos);
        }
        else
        {
            switch (*_cur_pos)
            {
            //If we get a "disconnect" symbol, check the stack for uncapped branches
            case '.':
                _bond_order = 0; //A gimmick in SMILES to code disconnects
                break;
            case '(':
                //Create a new branch
                if (_roots.size() > 0)
                {
                    _roots.push(_roots.top());
                }
                else
                {
                    //Throw error
                }
                break;
            case ')':
                //Close the current branch
                if (_roots.size() > 0)
                {
                    _roots.pop();
                }
                else
                {
                    //Throw error
                }
                break;
            case '%':
                _ring_closer.ProcessComplex(_cur_pos);
                _cur_pos += 2;
            case '-':
                //Set the bond order to single
                _bond_order = SINGLEBOND;
                break;
            case '=':
                //Set the bond order to double
                _bond_order = DOUBLEBOND;
                break;
            case '#':
                //Set the bond order to triple
                _bond_order = TRIPLEBOND;
                break;
            case ':':
                //Set the bond order to aromatic
                _bond_order = PARTIALDOUBLEBOND;
                break;
            case '/':
                _bond_stereo = STEREO_UP;
                break;
            case '\\':
                _bond_stereo = STEREO_DOWN;
                break;

            //These cases handle unbracketed _atoms, which in SMILES must
            //of these _atoms: C,N,O,S,P,F,Cl,Br,I,B,H,D,T,c,n,o,s,p
            case 'C':
                _cur_pos++;
                if (*_cur_pos == 'l' || *_cur_pos == 'L')
                {
                    // Add a chlorine atom
                    AddBasicAtom(17);
                }
                else
                {
                    _cur_pos--;
                    // Add a carbon atom
                    AddBasicAtom(6);
                }
                break;
            case 'N':
                //Add a nitrogen atom
                AddBasicAtom(7);
                break;
            case 'O':
                //Add an oxygen atom
                AddBasicAtom(8);
                break;
            case 'P':
                //Add a phosphorus atom
                AddBasicAtom(15);
                break;
            case 'S':
                //Add a sulfur atom
                AddBasicAtom(16);
                break;
            case 'F':
                //Add a fluorine atom
                AddBasicAtom(9);
                break;
            case 'B':
                _cur_pos++;
                if (*_cur_pos == 'r' || *_cur_pos == 'R')
                {
                    // Add a bromine atom
                    AddBasicAtom(35);
                }
                else
                {
                    _cur_pos--;
                    // Add a boron atom
                    AddBasicAtom(5);
                }
                break;
            case 'I':
                //Add an iodine atom
                AddBasicAtom(53);
                break;
            case 'H':
                //Add a hydrogen atom
                AddBasicAtom(1);
                break;
            case 'D':
                //Add a deuterium (2H) atom
                AddBasicAtom(1);
                break;
            case 'T':
                //Add a tritium atom
                AddBasicAtom(1);
                break;
            case 'c':
                //Add an aromatic carbon atom
                AddBasicAtom(6, true);
                break;
            case 'n':
                //Add an aromatic nitrogen atom
                AddBasicAtom(7, true);
                break;
            case 'o':
                //Add an aromatic oxygen atom
                AddBasicAtom(8, true);
                break;
            case 'p':
                //Add an aromatic phosphorus atom
                AddBasicAtom(15, true);
                break;
            case 's':
                //Add an aromatic sulfur atom
                AddBasicAtom(16, true);
                break;
            case '[':
                ProcessBracketAtom();
                break;
                // Unusual _atoms, or _atoms that need unusual specs (isotope,chiral,
                // charge, # of hydrogens) must be placed in brackets
                //		default:
                //Throw an error--this symbol is unrecognized
            } //End switch
        }  //End else

        _cur_pos++;
    } //End while

    //	_lig->atoms.assign(_atoms, _atoms + natoms);

    //Copy _atoms into the new "LIG" residue
    // make an array of MIAtom pointers and set them
    //	 if(_res->atoms() = (MIAtom * *)calloc(_res->natoms, sizeof(MIAtom *))){
    //		 int j;
    //		 for(j=0;j<_res->natoms;j++){
    //			 _res->atom(j)=_atoms[j];
    //		 }
    //	 }
    //	 else {
    //		//Throw memory allocation error
    //	 }
}

void SmilesReader::AddBasicAtom(int atomic_number, bool is_aromatic)
{
    AddAtom(atomic_number,
            is_aromatic,
            0,                                      //0 represents unspecified mass
            0,                                      //neutral charge assumed on "basic" atoms
            -1,                                     //-1 represents undetermined hydrogen count
            0,                                      //chiral class of 0 means no chirality specified
            0);                                     //chiral order immaterial when chiral_class=0
}

void SmilesReader::AddAtom(int atomic_number,
                           bool is_aromatic,
                           int mass,
                           int charge,
                           int hcount,
                           int chiral_class,
                           int chiral_order)
{

    MIAtom atom;
    //Create an atom. Init the properties using the parameters passed.
    //And push it onto the vector of _atoms.

    atom.setAtomicnumber(atomic_number);
    atom.setAtomnumber(_natoms + 1);
    atom.setIsaromatic(is_aromatic);
    atom.setMass(mass);
    atom.set_formal_charge(charge);
    atom.setHcount(hcount);
    atom.chiral_class(chiral_class);
    atom.chiral_order(chiral_order);

    atom.setName(format("%s%d", Left_Atomic_Name(atom.atomicnumber()), atom.atomnumber()).c_str());


    //Add this atom to the Ligand, capturing the address into which it is stored
    //(Should improve this...this address can be invalidated if the
    // STL vector is reallocated)
    MIAtom *atm_ptr = new MIAtom;
    atm_ptr->copyShallow(atom);
    _res->addAtom(atm_ptr);

    //Form a bond to the previous atom (if there is one)
    if (_roots.size() > 0 && _bond_order)
    {
        _lig->AddBond(_roots.top(), atm_ptr, _bond_order, _bond_stereo);
    }

    //Remove the last atom from the root/branch stack
    if (_roots.size() > 0)
    {
        _roots.pop();
    }

    //Reset bond parameters to default
    _bond_order = SINGLEBOND;
    _bond_stereo = 0;

    //Push this atom onto the root/branch stack, so a bond will be formed to
    //the next atom.
    _roots.push(atm_ptr);
    _natoms++;
}

//*	Function to parse the information from the "sublanguage of a bracketed
//* SMILES atom from the following syntax:
//*
//*   atom : '[' <mass> symbol <chiral> <hcount> <sign<charge>> ']'
//*       ;
void SmilesReader::ProcessBracketAtom()
{
    bool arom;                                  //Flag for aromaticity
    char elem_symbol[3];                        //One or two-letter symbol for the element type

    _cur_pos++;                                 //Advance beyond the initial '[' char

    //**Read atomic mass (optional)**//
    int mass = GetNumber();                     //Read leading digits (isotope information)


    //**Read atom symbol (required)**//
    if (islower(*_cur_pos))                     //MIAtom names in lower case imply aromaticity
    {
        arom = true;
        elem_symbol[0] = toupper(*_cur_pos++);
    }
    else
    {
        arom = false;
        elem_symbol[0] = *_cur_pos++;
    }


    if (islower(*_cur_pos))                     //If this is followed by a lower-case letter
    {
        elem_symbol[1] = toupper(*_cur_pos);    //that is the second letter of the symbol
        _cur_pos++;
    }
    else
    {
        elem_symbol[1] = elem_symbol[0];
        elem_symbol[0] = ' ';                  //Fill in with a space, to match our lookup tbl
    }
    elem_symbol[2] = '\0';

    int atomic_number = Atomic_Number(elem_symbol);     //Convert symbol to the atomic num


    //**Read chirality info (optional)**//

    int chiral_class = CH_NONE;                     //Chirality information is stored
    int chiral_order = 0;                           //with two bits, one for class, one for order


    if (*_cur_pos == '@')                           //The syntax for chirality info is
    {
        _cur_pos++;                                 // chival : '@' <chiclass> <chiorder>
        if (isdigit(*_cur_pos))                     //        | chiral '@'
        {
            chiral_class = CH_DEFAULT;              //        ;
            chiral_order = GetNumber();
        }
        else if (*_cur_pos == '@')             //"@@" is a shortcut for "@2", the second
        {
            chiral_class = CH_DEFAULT;              // order of the default class
            chiral_order = 2;
            _cur_pos++;
        }
        else if (*_cur_pos == 'T'              //These specs are rarely needed, since
                 && *(_cur_pos+1) == 'H')          //the default class is usually applied.
        {
            _cur_pos += 2;                          //Currently supports five classes.
            chiral_class = CH_TETRAHEDRAL;
            chiral_order = GetNumber();
        }
        else if (*_cur_pos == 'A'
                 && *(_cur_pos+1) == 'L')
        {
            _cur_pos += 2;
            chiral_class = CH_ALLENE_LIKE;
            chiral_order = GetNumber();
        }
        else if (*_cur_pos == 'S'
                 && *(_cur_pos+1) == 'P')
        {
            _cur_pos += 2;
            chiral_class = CH_SQUARE_PLANAR;
            chiral_order = GetNumber();
        }
        else if (*_cur_pos == 'T'
                 && *(_cur_pos+1) == 'B')
        {
            _cur_pos += 2;
            chiral_class = CH_TRIGONAL_BIPYRAMIDAL;
            chiral_order = GetNumber();
        }
        else if (*_cur_pos == 'O'
                 && *(_cur_pos+1) == 'H')
        {
            _cur_pos += 2;
            chiral_class = CH_OCTAHEDRAL;
            chiral_order = GetNumber();
        }
        else
        {
            chiral_class = CH_DEFAULT;              //"@" is a shortcut for "@1", the first
            chiral_order = 1;                       //order of the default class.
        }
    }

    //**Read attached hydrogen info (optional)**//
    int hcount = 0;                                     //Number of hydrogens on this atom

    if (*_cur_pos == 'H')
    {
        _cur_pos++;                                     //Hydrogen count is spec'd
        if (isdigit(*_cur_pos))                         //with an "H" followed by zero or more
        {
            hcount = GetNumber();                       //digits
        }
        else
        {
            hcount = 1;
        }
    }


    if (hcount == 1                                          //If this is a tetrahedral center
        && (_roots.size() == 0 || _bond_order == 0))        //with a hydrogen as the first substituent
    {
        chiral_order = (chiral_order == 1) ? 2 : 1;         //reverse the chiral order to reflect
    }                                                       //the change in reference axis

    //**Read formal charge (optional)**//
    int charge = 0;                                     //Charge of this atom

    if (*_cur_pos == '+')                           //Read positive charges
    {
        _cur_pos++;

        if (isdigit(*_cur_pos))                     //Charge is specified with a '+' followed
        {
            charge = GetNumber();                   //by zero or more digits
        }
        else if (*_cur_pos == '+')                 //"++" is a shortcut for "+2", "+++" is a
        {
            charge = 1 + CountRepetitions('+');    //'shortcut' for "+3", etc.
        }
        else
        {
            charge = 1;
        }
    }
    else if (*_cur_pos == '-')                    //Read negative charges
    {
        _cur_pos++;

        if (isdigit(*_cur_pos))                     //As above, charge is specified with a '-'
        {
            charge = -GetNumber();                  //followed by zero or more digits
        }
        else if (*_cur_pos == '-')
        {
            charge = -1 - CountRepetitions('-');
        }
        else
        {
            charge = -1;
        }
    }

    if (*_cur_pos == '\0')
    {
        //		cout << "WARNING, unexpectedly reached end of SMILES string"
        //			 << " while parsing a bracketed atom of type "  << elem_symbol << endl
        //			 << "Missing a ] character?" << endl;
    }
    if (*_cur_pos != ']')
    {
        //		cout << "WARNING, reached unexpected character \"" << *_cur_pos << "\""
        //			 << "when a \"]\" character was expected while parsing a bracketed atom of type "
        //			 << elem_symbol << endl
        //			 << "Missing a ] character?" << endl;
    }

    AddAtom(atomic_number, arom, mass, charge, hcount, chiral_class, chiral_order);
}

//Converts a string of digits to a (non-negative) integer value.
//The string of digits starts at the the current position (SmilesReader::_cur_pos),
//and ends with the first non-digit character.
int SmilesReader::GetNumber()
{
    int n = 0;
    while (isdigit(*_cur_pos))
    {
        n *= 10;
        n += *_cur_pos - '0';
        _cur_pos++;
    }
    return n;
}

//Counts the number of consecutive instances of a given character, starting
//at the current position.
int SmilesReader::CountRepetitions(char repeated)
{
    int n = 0;
    while (*_cur_pos == repeated)
    {
        n++;
        _cur_pos++;
    }
    return n;
}

//The Constructor
SmiRingClosures::SmiRingClosures()
{

    //Set all the half-bonds to be unoccupied
    int i;
    for (i = 0; i < 100; i++)
    {
        _open_half_bonds[i] = false;
        _short_half_bonds[i] = -1;
    }

}

//The Destructor
SmiRingClosures::~SmiRingClosures()
{
}

//Takes a character '0' through '9'
void SmiRingClosures::ProcessSimple(char key)
{
    _ringnum = (int) key - '0';
    if (_open_half_bonds[_ringnum])             //Check if this bond is already
    {
        Complete(_short_half_bonds[_ringnum]);  //half-done
    }
    else
    {
        Create(&_short_half_bonds[_ringnum]);
    }
}

void SmiRingClosures::ProcessComplex(char *key)
{
    _ringnum  = (int) 10 * (key[1] - '0');   //Convert the next two ASCI chars to
    _ringnum += (int) key[2] - '0';          //an integer

    if (_open_half_bonds[_ringnum])             //Check if this bond is already
    {
        Complete(_short_half_bonds[_ringnum]);  //half-done
    }
    else
    {
        Create(&_short_half_bonds[_ringnum]);
    }
}

void SmiRingClosures::Create(int *halfbond)
{
    Bond bond;


    if (_reader->_roots.size() > 0)
    {
        *halfbond = _reader->_lig->bonds.size();

        bond.setAtom1(_reader->_roots.top());  //Store a ptr to the atom,
        bond.setOrder(_reader->_bond_order);     //and two bond flags
        bond.stereo = _reader->_bond_stereo;
        bond.type = B_CONNECT;

        (_reader->_roots.top())->addBondnumber(_reader->_lig->bonds.size()); //Record this bond

        _reader->_bond_order = SINGLEBOND;      //Reset the bond flags
        _reader->_bond_stereo = 0;              //to default

        _open_half_bonds[_ringnum] = true;  //Store that we have an open
                                            //ring-closure here

        _reader->_lig->bonds.push_back(bond);   //Store this bond struct in the Ligand obj

    }
    else
    {
        //Throw error
    }
}

void SmiRingClosures::Complete(int halfbond)
{
    if (halfbond < 0)
    {
        exit(1);
    }

    Bond *bond;
    bond = &(_reader->_lig->bonds[halfbond]);

    //**Determine bond order.
    //There are two places where the bond order may be set--
    //creation or completion.  Set bond order to whichever is greater
    //than 1, but both should not be set differently.
    if (_reader->_bond_order == SINGLEBOND && bond->getOrder() != SINGLEBOND)
    {
        _reader->_bond_order = bond->getOrder();
    }
    else if (_reader->_bond_order != SINGLEBOND && bond->getOrder() != SINGLEBOND
             && _reader->_bond_order != bond->getOrder())
    {
        //Throw error
    }

    //**Determine bond stereo flags.
    //Same as above, there are two spots a flag could be set --
    //creation or completion.  Set stereo flag to whichever is
    //non-zero, but throw error if both are non-zero and different.
    if (_reader->_bond_stereo == 0 && bond->stereo > 0)
    {
        _reader->_bond_stereo = bond->stereo;
    }
    else if (_reader->_bond_stereo > 0 && bond->stereo > 0
             && _reader->_bond_stereo != bond->stereo)
    {
        //Throw error
    }

    bond->setAtom2(_reader->_roots.top());              //Complete the bond
    bond->getAtom1()->addNabor(bond->getAtom2());
    bond->getAtom2()->addNabor(bond->getAtom1());
    bond->getAtom2()->addBondnumber(halfbond);
    bond->setOrder(_reader->_bond_order);
    bond->stereo = _reader->_bond_stereo;

    _reader->_bond_order = SINGLEBOND;      //Reset the bond flags
    _reader->_bond_stereo = 0;              //to default

    _open_half_bonds[_ringnum] = false; //Done with this
                                        //ring-closure
    _short_half_bonds[_ringnum] = -1;

}

}
