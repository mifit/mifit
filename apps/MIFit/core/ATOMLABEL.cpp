#include "ATOMLABEL.h"
#include "RESIDUE.h"
#include <chemlib/Monomer.h>

using namespace chemlib;


int ATOMLABEL::defaultStyle_ = 0;
int ATOMLABEL::defaultSize_ = 12;
unsigned char ATOMLABEL::defaultRed_ = 255;
unsigned char ATOMLABEL::defaultGreen_ = 255;
unsigned char ATOMLABEL::defaultBlue_ = 255;

ATOMLABEL::ATOMLABEL()
    : atom_(NULL),
      residue_(NULL),
      useDefaultStyle_(true),
      style_(defaultStyle_),
      customLabel(false),
      visible_(false),
      xo(2),
      yo(-2),
      useDefaultColor_(true),
      red_(defaultRed_),
      green_(defaultGreen_),
      blue_(defaultBlue_),
      useDefaultSize_(true),
      size_(defaultSize_)
{

    label_ = labelString(residue_, atom_, style_);
}

ATOMLABEL::ATOMLABEL(const Monomer *residue, const MIAtom *atom)
    : atom_(atom),
      residue_(residue),
      useDefaultStyle_(true),
      style_(defaultStyle_),
      customLabel(false),
      visible_(true),
      xo(2),
      yo(-2),
      useDefaultColor_(true),
      red_(defaultRed_),
      green_(defaultGreen_),
      blue_(defaultBlue_),
      useDefaultSize_(true),
      size_(defaultSize_)
{

    label_ = labelString(residue_, atom_, style_);
}

const Monomer *ATOMLABEL::residue() const
{
    return residue_;
}

const MIAtom*ATOMLABEL::atom() const
{
    return atom_;
}

bool ATOMLABEL::useDefaultStyle() const
{
    return useDefaultStyle_;
}

void ATOMLABEL::useDefaultStyle(bool value)
{
    useDefaultStyle_ = value;
    style(defaultStyle_);
}

int ATOMLABEL::style() const
{
    if (useDefaultStyle_)
    {
        return defaultStyle_;
    }
    return style_;
}

void ATOMLABEL::style(int value)
{
    style_ = value;
    useDefaultStyle_ = false;
    label_ = labelString(residue_, atom_, style_);
}

bool ATOMLABEL::isCustomLabel()
{
    return customLabel;
}

bool ATOMLABEL::isVisible() const
{
    return visible_;
}

void ATOMLABEL::visible(bool on)
{
    visible_ = on;
}

const std::string&ATOMLABEL::label() const
{
    if (!customLabel && useDefaultStyle_ && defaultStyle_ != style_)
    {
        style_ = defaultStyle_;
        label_ = labelString(residue_, atom_, style_);
    }
    return label_;
}

void ATOMLABEL::label(const char *text)
{
    customLabel = true;
    label_ = text;
}

int ATOMLABEL::xOffset() const
{
    return xo;
}

void ATOMLABEL::xOffset(int offset)
{
    xo = offset;
}

int ATOMLABEL::yOffset() const
{
    return yo;
}

void ATOMLABEL::yOffset(int offset)
{
    yo = offset;
}

bool ATOMLABEL::useDefaultColor() const
{
    return useDefaultColor_;
}

void ATOMLABEL::useDefaultColor(bool value)
{
    useDefaultColor_ = value;
}

unsigned char ATOMLABEL::red() const
{
    if (useDefaultColor_)
    {
        return defaultRed_;
    }
    return red_;
}

void ATOMLABEL::red(unsigned char value)
{
    red_ = value;
    useDefaultColor_ = false;
}

unsigned char ATOMLABEL::green() const
{
    if (useDefaultColor_)
    {
        return defaultGreen_;
    }
    return green_;
}

void ATOMLABEL::green(unsigned char value)
{
    green_ = value;
    useDefaultColor_ = false;
}

unsigned char ATOMLABEL::blue() const
{
    if (useDefaultColor_)
    {
        return defaultBlue_;
    }
    return blue_;
}

void ATOMLABEL::blue(unsigned char value)
{
    blue_ = value;
    useDefaultColor_ = false;
}

bool ATOMLABEL::useDefaultSize() const
{
    return useDefaultSize_;
}

void ATOMLABEL::useDefaultSize(bool value)
{
    useDefaultSize_ = value;
}

int ATOMLABEL::size() const
{
    if (useDefaultSize_)
    {
        return defaultSize_;
    }
    return size_;
}

void ATOMLABEL::size(int value)
{
    size_ = value;
}

int ATOMLABEL::defaultStyle()
{
    return defaultStyle_;
}

void ATOMLABEL::defaultStyle(int value)
{
    defaultStyle_ = value;
}

void ATOMLABEL::defaultColor(unsigned char red, unsigned char green, unsigned char blue)
{
    defaultRed_ = red;
    defaultGreen_ = green;
    defaultBlue_ = blue;
}

unsigned char ATOMLABEL::defaultRed()
{
    return defaultRed_;
}

unsigned char ATOMLABEL::defaultGreen()
{
    return defaultGreen_;
}

unsigned char ATOMLABEL::defaultBlue()
{
    return defaultBlue_;
}

int ATOMLABEL::defaultSize()
{
    return defaultSize_;
}

void ATOMLABEL::defaultSize(int value)
{
    defaultSize_ = value;
}

std::string ATOMLABEL::labelString(const Monomer *res, const MIAtom *atom, int style)
{

    static char id[128];
    const std::string &t = res->type();
    const std::string &r = res->name();
    const std::string &a = atom->name();
    char l = atom->altloc();
    char c = (char)(res->chain_id()&255);
    char s;

    if (t.length() == 0 && r.length() == 0 && a.length() == 0)
    {
        strcpy(id, "???");
        return (id);
    }

    switch (style)
    {
    default:
    case 0:
        if (c != ' ')
        {
            if (l == ' ')
            {
                sprintf(id, "%s %s %c %s", t.c_str(), r.c_str(), c, a.c_str());
            }
            else
            {
                sprintf(id, "%s %s %c %s_%c", t.c_str(), r.c_str(), c, a.c_str(), l);
            }
        }
        else
        if (l == ' ')
        {
            sprintf(id, "%s %s %s", t.c_str(), r.c_str(), a.c_str());
        }
        else
        {
            sprintf(id, "%s %s %s_%c", t.c_str(), r.c_str(), a.c_str(), l);
        }
        break;
    case 1:
        if (c != ' ')
        {
            sprintf(id, "%s %s %c", t.c_str(), r.c_str(), c);
        }
        else
        {
            sprintf(id, "%s %s", t.c_str(), r.c_str());
        }
        break;
    case 2:
        sprintf(id, "%s %s", t.c_str(), r.c_str());
        break;
    case 3:
        sprintf(id, "%s", t.c_str());
        break;
    case 4:
        sprintf(id, "%s", r.c_str());
        break;
    case 5:
        if (l == ' ')
        {
            sprintf(id, "%s", a.c_str());
        }
        else
        {
            sprintf(id, "%s_%c", a.c_str(), l);
        }
        break;
    case 6:
        s = singleletter(t.c_str());
        if (c != ' ')
        {
            sprintf(id, "%c%s %c", s, r.c_str(), c);
        }
        else
        {
            sprintf(id, "%c%s", s, r.c_str());
        }
        break;
    case 7:
        s = singleletter(t.c_str());
        sprintf(id, "%c%s", s, r.c_str());
        break;
    case 8:
        s = singleletter(t.c_str());
        sprintf(id, "%c", s);
        break;
    }
    return std::string(id);
}

