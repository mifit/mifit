#ifndef mi_opengl_Text_h
#define mi_opengl_Text_h

#include <QtGlobal>

class QChar;
class QFont;
class QFontMetrics;
class QString;

namespace mi
{
namespace opengl
{

class TextPrivate;

class Text
{
public:
    Text(const QFont &f);
    virtual ~Text();

    QFont font() const;
    QFontMetrics fontMetrics() const;

    void renderText(float x, float y, const QString &text);


private:
    Q_DISABLE_COPY(Text)

    TextPrivate *const d;
};

} // namespace opengl
} // namespace mi

#endif // mi_opengl_Text_h
