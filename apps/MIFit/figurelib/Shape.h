#ifndef mifit_figurelib_Shape_h
#define mifit_figurelib_Shape_h

#include <vector>

namespace moldraw
{
    class Drawing;

    class Shape
    {
    public:
        Shape(Drawing *dp);
        virtual ~Shape();

        virtual void Draw() = 0;

    protected:
        Drawing *_dp;
        float _line_width;
    };


}

#endif // ifndef mifit_figurelib_Shape_h
