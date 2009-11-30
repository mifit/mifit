#ifndef mi_opengl_GLFont_h
#define mi_opengl_GLFont_h

#include <QFont>
#include <QFontMetrics>
#include <QGLWidget>
#include <QImage>
#include <QPainter>
#include <QPen>
#include <QString>

namespace mi
{
    namespace opengl
    {

        class GLFont
        {

            QGLWidget *glWidget_;
            QPen pen_;
            QFont font_;
            QFontMetrics fontMetrics_;
            QPainter painter_;
            QImage image_;

        public:

            GLFont(QGLWidget *widget);
            virtual ~GLFont();

            GLFont(const GLFont &font);
            GLFont&operator=(const GLFont &font);

            void setGlWidget(QGLWidget *widget)
            {
                glWidget_ = widget;
            }

            void render(const QString &text);

            const QPen&pen() const
            {
                return pen_;
            }

            void setPen(QPen pen)
            {
                pen_ = pen;
            }

            const QFont&font() const
            {
                return font_;
            }

            void setFont(const QFont &font)
            {
                font_ = font;
            }

            const QFontMetrics&fontMetrics() const
            {
                return fontMetrics_;
            }

        };

    }
}

#endif // ifndef mi_opengl_GLFont_h
