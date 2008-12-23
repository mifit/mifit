#ifndef mifit_figurelib_Drawing_h
#define mifit_figurelib_Drawing_h

#include <vector>
#include <string>

#include "corelib.h"

#ifndef RAD2DEG
#define RAD2DEG (57.295779513082320876798154814105)     /* 180 / PI */
#endif
#ifndef DEG2RAD
#define DEG2RAD (0.017453292519943295769236907684886)   /* PI / 180 */
#endif

namespace moldraw {

//NOTE: this class now encapsulates everything to be drawn in a simple
//      vector of draw_item objects.  It can then be trivially subclassed
//      to render the item list in a QCanvas or OpenGL canvas, etc.

class Drawing {
public:
  virtual ~Drawing() {}

  static const float ANNOTATION_LABEL_SIZE;
  static const float RESIDUE_LABEL_SIZE;
  static const float HBOND_LABEL_SIZE;
  static const float SUN_LABEL_SIZE;
  static const float ATOM_LABEL_SIZE;
  static const float ATOM_LABEL_OFFSET;
  static const float ATOM_RADIUS;
  static const float UNIRES_RADIUS;

  static const float ATOM_SPOKE_WIDTH;
  static const float SUN_BORDER_WIDTH;
  static const float STD_BOND_WIDTH;
  static const float LIG_BOND_WIDTH;
  static const float SPOKE_EXTENT;

  static const int MAX_SUN_CYCLES;
  static const int MAX_RESLABEL_CYCLES;

  virtual void FitToPage(std::vector<float> origin, std::vector<float> dimensions) = 0;
  virtual void Finish() = 0;

  virtual void DrawAtom(float x, float y, float r, const PaletteColor& color);
  virtual void DrawSingleBond(float x1, float y1, float x2, float y2, float width);
  virtual void DrawDoubleBond(float x1, float y1, float x2, float y2, float width);
  virtual void DrawHydrogenBond(float x1, float y1, float x2, float y2, float width, float font_size, float dist);
  virtual void DrawSun(float x, float y, float r, float width, float font_size, const std::string& name);
  virtual void DrawSpokes(float x1, float y1, float r1, float x2, float y2, float r2, float, float width1, float width2);
  virtual void DrawLabel(float x, float y, float off, float off_x, float off_y, float font_size, const std::string& content);


  protected:
    class draw_item {
      public:
        enum ItemType { Atom, SingleBond, DoubleBond, HydrogenBond, Sun, Spokes, Label };

        draw_item() : x1(0.0f), y1(0.0f), x2(0.0f), y2(0.0f), width(0.0f), width2(0.0f),
                      r(0.0f), r2(0.0f), font_size(0.0f), dist(0.0f), off(0.0f), 
                      off_x(0.0f), off_y(0.0f), text(""), font(""), color(PaletteColor(0,0,0)) {}

        ItemType type;
        float x1, y1, x2, y2, width, width2, r, r2, font_size, dist, off, off_x, off_y;
        std::string text;
        std::string font;
        PaletteColor color;
    };

    std::vector<draw_item> _items;

};

}

#endif
