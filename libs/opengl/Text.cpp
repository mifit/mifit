#include "Text.h"

#include <cmath>
#include <iostream>
#include <QtCore/QHash>
#include <QtCore/QSysInfo>
#include <QtGui/QPainter>
#include <QtGui/QPixmap>
#include <QtOpenGL/QGLFormat>
#include <QtOpenGL/QGLFramebufferObject>

namespace
{
    const int TEXTURE_SIZE = 256;

    struct CharData
    {
        GLuint textureId;
        uint width;
        uint height;
        GLfloat s[2];
        GLfloat t[2];
    };

} // anonymous namespace


namespace mi
{
namespace opengl
{

struct TextPrivate
{
    TextPrivate(const QFont &f);
    ~TextPrivate();

    void allocateTexture();
    CharData &createCharacter(QChar c);

    QFont font;
    QFontMetrics fontMetrics;

    QHash<ushort, CharData> characters;
    QList<GLuint> textures;

    GLint xOffset;
    GLint yOffset;
};

TextPrivate::TextPrivate(const QFont &f)
    : font(f), fontMetrics(f), xOffset(0), yOffset(0)
{
}

TextPrivate::~TextPrivate()
{
    foreach (GLuint texture, textures)
        glDeleteTextures(1, &texture);
}

void TextPrivate::allocateTexture()
{
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    QImage image(TEXTURE_SIZE, TEXTURE_SIZE, QImage::Format_ARGB32);
    image.fill(Qt::transparent);
    image = QGLWidget::convertToGLFormat(image);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, TEXTURE_SIZE, TEXTURE_SIZE,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, image.bits());

    textures += texture;
}

CharData &TextPrivate::createCharacter(QChar c)
{
    ushort unicodeC = c.unicode();
    if (characters.contains(unicodeC))
        return characters[unicodeC];

    if (textures.empty())
        allocateTexture();

    GLuint texture = textures.last();

    GLsizei width = fontMetrics.width(c);
    GLsizei height = fontMetrics.height();

    QPixmap pixmap(width, height);
    pixmap.fill(Qt::transparent);
    QPainter painter;
    painter.begin(&pixmap);
    painter.setRenderHints(QPainter::HighQualityAntialiasing
                           | QPainter::TextAntialiasing);
    painter.setFont(font);
    painter.setPen(Qt::white);

    painter.drawText(0, fontMetrics.ascent(), c);
    painter.end();


    QImage image = QGLWidget::convertToGLFormat(pixmap.toImage());
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, xOffset, yOffset, width, height, GL_RGBA,
                    GL_UNSIGNED_BYTE, image.bits());

    CharData& character = characters[unicodeC];
    character.textureId = texture;
    character.width = width;
    character.height = height;
    character.s[0] = static_cast<GLfloat>(xOffset) / TEXTURE_SIZE;
    character.t[0] = static_cast<GLfloat>(yOffset) / TEXTURE_SIZE;
    character.s[1] = static_cast<GLfloat>(xOffset + width) / TEXTURE_SIZE;
    character.t[1] = static_cast<GLfloat>(yOffset + height) / TEXTURE_SIZE;

    xOffset += width;
    if (xOffset + fontMetrics.maxWidth() >= TEXTURE_SIZE)
    {
        xOffset = 1;
        yOffset += height;
    }
    if (yOffset + fontMetrics.height() >= TEXTURE_SIZE)
    {
        allocateTexture();
        yOffset = 1;
    }
    return character;
}

Text::Text(const QFont &f) : d(new TextPrivate(f))
{
}

Text::~Text()
{
    delete d;
}

QFont Text::font() const
{
    return d->font;
}

QFontMetrics Text::fontMetrics() const
{
    return d->fontMetrics;
}

//! Renders text at given x, y.
void Text::renderText(float x, float y, const QString &text)
{
    // If the current context's device is not active for painting, the
    // texture generation does not work. This may be specific to the way
    // MIFit is setup.
    if (!QGLContext::currentContext()->device()->paintingActive())
        return;

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_TEXTURE_BIT);
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    GLuint texture = 0;
    glTranslatef(x, y, 0);
    for (int i = 0; i < text.length(); ++i)
    {
        CharData &c = d->createCharacter(text[i]);
        if (texture != c.textureId)
        {
            texture = c.textureId;
            glBindTexture(GL_TEXTURE_2D, texture);
        }


        glBegin(GL_QUADS);
        glTexCoord2f(c.s[0], c.t[0]);
        glVertex2f(0, 0);

        glTexCoord2f(c.s[1], c.t[0]);
        glVertex2f(c.width, 0);

        glTexCoord2f(c.s[1], c.t[1]);
        glVertex2f(c.width, c.height);

        glTexCoord2f(c.s[0], c.t[1]);
        glVertex2f(0, c.height);
        glEnd();

        glTranslatef(c.width, 0, 0);
    }

    glPopMatrix();
    glPopAttrib();
}

} // namespace opengl
} // namespace mi
