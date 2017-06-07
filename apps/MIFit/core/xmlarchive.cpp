#include <nongui/nonguilib.h>
#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>
#include <map/maplib.h>

#include "xmlarchive.h"
#include "Annotation.h"

using namespace chemlib;

//////////////////////////////////////////////////////////////////////
// XMLArchive Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

static std::string encode(const char *str)
{
    std::string out;
    int length = strlen(str);
    for (int i = 0; i < length; ++i)
    {
        switch (str[i])
        {
        case '"':
            out += "&quot;";
            break;
        case '&':
            out += "&amp;";
            break;
        case '\'':
            out += "&apos;";
            break;
        case '<':
            out += "&lt;";
            break;
        case '>':
            out += "&gt;";
            break;
        default:
            out += str[i];
            break;
        }
    }
    return out;
}

bool IsXML(FILE *fp)
{
    bool found_xml = false;
    int count = 0;
    std::string buf;
    std::unique_ptr<io> ioObj(io::defaultIo());
    io &file = *ioObj;
    file.attach(fp);
    file.rewind();
    while (file.readLine(buf) != 0 && count < 20)
    {
        if (strstr(buf.c_str(), "<?xml ") )
        {
            found_xml = true;
        }
        count++;
    }
    file.rewind();
    return found_xml;
}

XMLArchive::XMLArchive(const char *pathname, int m)
    : CArchive(pathname, m)
{
    inTag = false;
    f = format("<?xml version = \"1.0\" encoding=\"UTF-8\"?>\n\n");
    Write(f);
    //f=format("<!-- Created by MIFit XMLArchive. User: %s on %s -->\n\n", ::wxGetUserId().c_str(), ::wxGetHostName().c_str() );
    //Write(f);
}

XMLArchive::~XMLArchive()
{

}

void XMLArchive::WriteField(const char *tag, int value)
{
    calc_indent();
    f = format("%s<%s>%d</%s>\n", indent, tag, value, tag);
    Write(f);
}

void XMLArchive::WriteField(const char *tag, float value)
{
    calc_indent();
    f = format("%s<%s>%f</%s>\n", indent, tag, value, tag);
    Write(f);
}

void XMLArchive::WriteField(const char *tag, const char *value)
{
    calc_indent();
    f = format("%s<%s>%s</%s>\n", indent, tag, encode(value).c_str(), tag);
    Write(f);
}

void XMLArchive::WriteField(const char *tag)
{
    /* write an empty field */
    calc_indent();
    f = format("%s<%s/>\n", indent, tag);
    Write(f);
}

void XMLArchive::WriteHBond(Bond &hbond)
{
    /* write an hbond */
    calc_indent();
    f = format("%s<H_Bond atom1=\"%d\" atom2=\"%d\" />\n", indent, (int)hbond.getAtom1()->atomnumber(), (int)hbond.getAtom2()->atomnumber());
    Write(f);
}

void XMLArchive::WriteConnect(Bond &bond)
{
    /* write a connect */
    calc_indent();
    f = format("%s<Connect atom1=\"%d\" atom2=\"%d\" />\n", indent, (int)bond.getAtom1()->atomnumber(), (int)bond.getAtom2()->atomnumber());
    Write(f);
}

void XMLArchive::WriteColor(PaletteColor &c)
{
    /* write a connect */
    calc_indent();
    f = format("%s<Color red=\"%d\" green=\"%d\" blue=\"%d\" />\n", indent, (int)c.red, (int)c.green, (int)c.blue);
    Write(f);
}

void XMLArchive::WriteLabel(const ATOMLABEL &label)
{
    /* write a label */
    calc_indent();
    f = format("%s<Label atom=\"%d\" label=\"%s\" red=\"%d\" green=\"%d\" blue=\"%d\" visible=\"%d\" xo=\"%d\" yo=\"%d\" />\n",
               indent, (int)label.atom()->atomnumber(), label.label().c_str(), (int)label.red(),
               (int)label.green(), (int)label.blue(), (int)label.isVisible(), label.xOffset(), label.yOffset());
    Write(f);
}

void XMLArchive::WriteDot(SURFDOT &dot)
{
    /* write a dot */
    calc_indent();
    f = format("%s<SurfDot color=\"%d\" x=\"%f\" y=\"%f\" z=\"%f\" width=\"%d\"/>\n", indent, (int)dot.color, dot.x, dot.y, dot.z, (int)dot.w);
    Write(f);
}

void XMLArchive::WriteSymmCenter(float center[3])
{
    /* write symm_center */
    calc_indent();
    f = format("%s<SymmCenter x=\"%f\" y=\"%f\" z=\"%f\" />\n", indent, center[0], center[1], center[2]);
    Write(f);
}

void XMLArchive::WriteAnnotation(Annotation *ann)
{
    /* write an annotation */
    ann->WriteXML(this);
}

void XMLArchive::WriteMapHeader(CMapHeaderBase &mh)
{
    /* write a mapheader */
    if (!MapHeaderOK(&mh))
    {
        return;
    }
    BeginTag("MapHeader");
    calc_indent();
    f = format("%s<UnitCell a=\"%f\" b=\"%f\" c=\"%f\" alpha=\"%f\" beta=\"%f\" gamma=\"%f\" />\n", indent,
               mh.a, mh.b, mh.c, mh.alpha, mh.beta, mh.gamma);
    Write(f);
    f = format("%s<Limits hmin=\"%d\" hmax=\"%d\" kmin=\"%d\" kmax=\"%d\" lmin=\"%d\" lmax=\"%d\" resmin=\"%f\" resmax=\"%f\" />\n", indent,
               mh.hmin, mh.hmax, mh.kmin, mh.kmax, mh.lmin, mh.lmax, mh.resmin, mh.resmax);
    Write(f);
    f = format("%s<SymmOps spgpno=\"%d\" />\n", indent, mh.spgpno);
    Write(f);
    f = format("%s<MapHeaderOptions crystal_name=\"%s\" M_Coefficient=\"%f\" N_Coefficient=\"%f\" fc_is_fom=\"%d\" Bsolvent=\"%f\" Ksolvent=\"%f\" use_aniso=\"%d\" use_bulksolvent=\"%d\" title=\"%s\" />\n",
               indent, mh.crystal_name.c_str(), mh.M_Coefficient, mh.N_Coefficient,
               (int)mh.fc_is_fom, mh.Bsolvent, mh.Ksolvent, (int)mh.use_aniso, (int)mh.use_bulksolvent, mh.title.c_str());
    Write(f);
    if (mh.nNCRSymmops > 0)
    {
        for (int i = 0; i < mh.nNCRSymmops; i++)
        {
            f = format("%s<NCRSymmop n=\"%d\" m1=\"%f\" m2=\"%f\" m3=\"%f\" m4=\"%f\" m5=\"%f\" m6=\"%f\" m7=\"%f\" m8=\"%f\" m9=\"%f\" m10=\"%f\" m11=\"%f\" m12=\"%f\" />\n", indent,
                       i+1, mh.NCRSymmops[i][0], mh.NCRSymmops[i][1], mh.NCRSymmops[i][2], mh.NCRSymmops[i][3], mh.NCRSymmops[i][4],
                       mh.NCRSymmops[i][5], mh.NCRSymmops[i][6], mh.NCRSymmops[i][7], mh.NCRSymmops[i][8], mh.NCRSymmops[i][9],
                       mh.NCRSymmops[i][10], mh.NCRSymmops[i][11]);
            Write(f);
        }
    }
    EndTag();
}

void XMLArchive::WriteBond(Bond &hbond)
{
    /* write an hbond */
    calc_indent();
    f = format("%s<Bond atom1=\"%d\" atom2=\"%d\" type=\"%d\" />\n", indent, (int)hbond.getAtom1()->atomnumber(), (int)hbond.getAtom2()->atomnumber(), (int)hbond.type);
    Write(f);
}

void XMLArchive::WriteField(const MIAtom *a)
{
    /* write an ATOM field */
    calc_indent();
    char altloc = a->altloc();
    if (!isprint(altloc))
    {
        altloc = ' ';
    }
    f = format("%s<Atom name=\"%s\" x=\"%f\" y=\"%f\" z=\"%f\" color=\"%d\" altloc=\"%c\" occupancy=\"%f\" BValue=\"%f\" number=\"%d\" atomicnumber=\"%d\"",
               indent, a->name(), a->x(), a->y(), a->z(), (int)a->color(), altloc, a->occ(), a->BValue(),
               (int)a->atomnumber(), (int)a->atomicnumber());
    if (a->hasAnisotropicity())
    {
        std::string s;
        s = format(" U11=\"%f\" U22=\"%f\" U33=\"%f\" U12=\"%f\" U13=\"%f\" U23=\"%f\"",
                   a->U(0), a->U(1), a->U(2), a->U(3), a->U(4), a->U(5));
        f += s;
    }
    f += "/>\n";
    Write(f);
}

void XMLArchive::WriteField(const Residue &res)
{
    /* write a RESIDUE field */
    tags.push_back("Residue");
    calc_indent();
    char secstr = res.secstr();
    // this fixes a wierd bug - if these happen to be 0 - disaster!
    if (!isalpha(secstr))
    {
        secstr = 'U';
    }
    if (!isalpha(res.name1()))
    {
        secstr = 'X';
    }
    f = format("%s<Residue name=\"%s\" type=\"%s\" chain=\"%d\" secstr=\"%c\" confomer=\"%d\" flags=\"%d\" linkage_type=\"%d\" name1=\"%c\" seqpos=\"%d\" >\n",
               indent, res.name().c_str(), res.type().c_str(), (int)res.chain_id(), secstr, (int)res.confomer(), (int)res.flags(), (int)res.linkage_type(), res.name1(), res.seqpos());
    Write(f);
    inTag++;
    for (int i = 0; i < res.atomCount(); i++)
    {
        WriteField(res.atom(i));
    }
    EndTag();
}

void XMLArchive::BeginTag(const char *t)
{
    tags.push_back(t);
    calc_indent();
    f = format("%s<%s>\n", indent, t);
    Write(f);
    inTag++;
}

void XMLArchive::BeginTag(const char *t, const char *attr)
{
    tags.push_back(t);
    calc_indent();
    f = format("%s<%s %s>\n", indent, t, attr);
    Write(f);
    inTag++;
}

void XMLArchive::EndTag()
{
    inTag--;
    if (inTag < 0)
    {
#ifdef DEBUG
        Logger::message("Warning: parentheses level in XMLArchive less than 0: programmer error");
#endif
        inTag = 0;
    }
    calc_indent();
    f = format("%s</%s>\n", indent, tags[inTag].c_str());
    Write(f);
    tags.erase(tags.begin() + inTag);
}

void XMLArchive::calc_indent()
{
    int i;
    int ns = (int) inTag * 8;
    if ( (unsigned int) ns >= sizeof(indent))
    {
        ns = sizeof(indent -1);
    }
    if (ns < 0)
    {
        ns = 0;
    }
    if (ns == 0)
    {
        indent[0] = '\0';
    }
    else
    {
        for (i = 0; i < ns; i++)
        {
            indent[i] = ' ';
        }
        indent[ns] = '\0';
    }
}

