#ifndef mifit_figurelib_Spokes_h
#define mifit_figurelib_Spokes_h

#include "Shape.h"

namespace moldraw
{

    class Spokes : Shape
    {
    public:
        Spokes(Drawing *dp,
               float x1,
               float y1,
               float r1,
               float x2,
               float y2,
               float r2,
               float extent,
               float width1,
               float width2);
        virtual ~Spokes()
        {
        }

        virtual void Draw();
    private:
        float _x1;
        float _y1;
        float _r1;
        float _x2;
        float _y2;
        float _r2;
        float _extent;
        float _width1;
        float _width2;
    };

}

#endif // ifndef mifit_figurelib_Spokes_h
