#include "valence.h"

namespace chemlib
{
std::vector<int> GetValenceStates(int atomicnumber)
{

    unsigned int vt;                //The classification of the valence type
    vt = GetValenceType(atomicnumber);

    std::vector<int> valences;
    switch (vt)
    {
    case ValenceType::vtH:                                  //The first four types have only one
        valences.push_back(1);                              //valence state.  The client can just
        return valences;                                    //add hydrogen atoms until this valence is
    case ValenceType::vtO:                                  //reached
        valences.push_back(2);
        return valences;
    case ValenceType::vtB:
        valences.push_back(3);
        return valences;
    case ValenceType::vtC:
        valences.push_back(4);
        return valences;
    case ValenceType::vtN:                                  //Nitrogen and sulfur can take on different
        valences.push_back(3);                              //valence states.  The client will assume
        valences.push_back(5);                              //the valence to be the *lowest* of these that
        return valences;                                    //is consistent with the current bonding pattern
    case ValenceType::vtS:
        valences.push_back(2);
        valences.push_back(4);
        valences.push_back(6);
        return valences;
    case ValenceType::vtI:
        valences.push_back(1);
        valences.push_back(3);
        valences.push_back(5);
        valences.push_back(7);
        return valences;
    default:
        return valences;
    }
}

unsigned int GetValenceType(int atomicnumber)
{

    switch (atomicnumber)
    {
    //		case 1: return ValenceType::vtH;
    //		case 3: return ValenceType::vtH;
    //		case 4: return ValenceType::vtO;
    case 5: return ValenceType::vtB;
    case 6: return ValenceType::vtC;            //carbon
    case 7: return ValenceType::vtN;            //nitrogen
    case 8: return ValenceType::vtO;            //oxygen
    //		case 9: return ValenceType::vtH;		//fluorine
    //		case 11: return ValenceType::vtH;
    //		case 12: return ValenceType::vtO;
    //		case 13: return ValenceType::vtB;
    case 14: return ValenceType::vtC;           //silicon
    case 15: return ValenceType::vtN;           //phosphorus
    case 16: return ValenceType::vtS;           //sulfur
    //		case 17: return ValenceType::vtI;		//chlorine
    //		case 19: return ValenceType::vtH;
    //		case 20: return ValenceType::vtO;
    //		case 31: return ValenceType::vtB;
    case 32: return ValenceType::vtC;
    case 33: return ValenceType::vtN;           //arsenic
    case 34: return ValenceType::vtS;           //selenium
    //		case 35: return ValenceType::vtI;		//bromine
    //		case 37: return ValenceType::vtH;		//rubidium
    //		case 38: return ValenceType::vtO;		//strontium
    //		case 49: return ValenceType::vtB;		//indium
    //		case 51: return ValenceType::vtN;
    //		case 52: return ValenceType::vtS;		//tellurium
    //		case 53: return ValenceType::vtI;		//iodine
    //		case 55: return ValenceType::vtH;		//cesium
    //		case 56: return ValenceType::vtO;		//barium
    //		case 83: return ValenceType::vtN;
    //		case 84: return ValenceType::vtS;
    //		case 85: return ValenceType::vtI;		//astatine
    //		case 87: return ValenceType::vtH;		//francium
    //		case 88: return ValenceType::vtO;		//radium
    default: return 0;
    }
}

} //namespace chemlib
