#ifndef mifit_model_ATOMLABEL_H
#define mifit_model_ATOMLABEL_H

#include <string>
namespace chemlib
{
    class MIAtom;
    class Monomer;
}

/**
 * An atom label structure.
 */
class ATOMLABEL
{

    /**
     * Pointer to atom labeled. If atom x,y,z is changed label follows.
     */
    const chemlib::MIAtom *atom_;

    const chemlib::Monomer *residue_;

    bool useDefaultStyle_;

    mutable int style_;

    bool customLabel;

    mutable std::string label_;

    bool visible_;

    /**
     * pixel offset from computed x y position
     */
    int xo;

    /**
     * pixel offset from computed x y position
     */
    int yo;

    bool useDefaultColor_;
    unsigned char red_;
    unsigned char green_;
    unsigned char blue_;

    bool useDefaultSize_;
    int size_;

    static int defaultStyle_;
    static unsigned char defaultRed_;
    static unsigned char defaultGreen_;
    static unsigned char defaultBlue_;
    static int defaultSize_;

public:

    ATOMLABEL();

    ATOMLABEL(const chemlib::Monomer *residue, const chemlib::MIAtom *atom);

    const chemlib::Monomer *residue() const;

    const chemlib::MIAtom *atom() const;

    bool useDefaultStyle() const;

    void useDefaultStyle(bool value);

    int style() const;

    void style(int value);

    bool isCustomLabel();

    bool isVisible() const;

    void visible(bool on);

    const std::string &label() const;
    void label(const char *text);

    int xOffset() const;
    void xOffset(int offset);

    int yOffset() const;
    void yOffset(int offset);

    bool useDefaultColor() const;
    void useDefaultColor(bool value);

    unsigned char red() const;
    void red(unsigned char value);

    unsigned char green() const;
    void green(unsigned char value);

    unsigned char blue() const;
    void blue(unsigned char value);

    bool useDefaultSize() const;
    void useDefaultSize(bool value);

    int size() const;
    void size(int value);

    static void defaultColor(unsigned char red, unsigned char green, unsigned char blue);
    static unsigned char defaultRed();
    static unsigned char defaultGreen();
    static unsigned char defaultBlue();
    static int defaultStyle();
    static void defaultStyle(int value);
    static int defaultSize();
    static void defaultSize(int value);

    static std::string labelString(const chemlib::Monomer *res, const chemlib::MIAtom *atom, int style);
};

#endif // ifndef mifit_model_ATOMLABEL_H
