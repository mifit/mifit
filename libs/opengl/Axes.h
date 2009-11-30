#ifndef mi_opengl_Axes_h
#define mi_opengl_Axes_h

#include <memory>

namespace mi
{
    namespace opengl
    {

        class GLFont;

        class Axes
        {

            static int labelsFontSize;

            std::auto_ptr<GLFont> labelsFont;

            int fontSize;

            int length;

            float glUnitsPerPixel;

            void drawText(const char *text, float x, float y, float z, float offset, float scale, float inverseRotation[16]);

        public:

            Axes(std::auto_ptr<GLFont> font);

            void setFontSize(int fontSize);

            void setLength(int length);

            void setGlUnitsPerPixel(float glUnitsPerPixel);

            void render();

        };

    }
}

#endif // ifndef mi_opengl_Axes_h
