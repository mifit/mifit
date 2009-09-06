#include "GLFont.h"

namespace mi {
namespace opengl {

GLFont::GLFont(QGLWidget* widget)
  : glWidget_(widget),
    pen_(Qt::white, 2),
    font_("Helvetica", 12, QFont::Bold),
    fontMetrics_(font_),
    painter_(),
    image_(10, 10, QImage::Format_ARGB32) {

  // Determine size for rendered text
  painter_.begin(&image_);
  painter_.setRenderHint(QPainter::Antialiasing, true);
  painter_.setRenderHint(QPainter::TextAntialiasing, true);
  painter_.setFont(font_);
  painter_.setPen(pen_);
  fontMetrics_ = painter_.fontMetrics();
  painter_.end();

}

GLFont::~GLFont() {
}

GLFont::GLFont(const GLFont& font)
: glWidget_(font.glWidget_),
  pen_(font.pen_),
  font_(font.font_),
  fontMetrics_(font.fontMetrics_),
  painter_(),
  image_(font.image_.size(), QImage::Format_ARGB32) {
}

GLFont& GLFont::operator=(const GLFont& font) {
  if (&font != this) {
    glWidget_ = font.glWidget_;
    pen_ = font.pen_;
    font_ = font.font_;
    fontMetrics_ = font.fontMetrics_;
  }
  return *this;
}

void GLFont::render(const QString& text) {

  if (!glWidget_) {
    return;
  }
  QRect fontRect = fontMetrics_.boundingRect(text);

  // Create proper size image and render text
  image_ = QImage(fontRect.size(), QImage::Format_ARGB32);
  image_.fill(0);
  painter_.begin(&image_);
  painter_.setRenderHint(QPainter::Antialiasing, true);
  painter_.setRenderHint(QPainter::TextAntialiasing, true);
  painter_.setFont(font_);
  painter_.setPen(pen_);
  painter_.drawText(0, image_.height(), text);
  painter_.end();

  glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT);

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Transfer text image to OpenGL
  QImage finalImage = image_.mirrored(false, true);
  GLuint id = glWidget_->bindTexture(finalImage);
  glWidget_->drawTexture(QPointF(0.0, 0.0), id);
  glWidget_->deleteTexture(id);

  glPopAttrib();
}

}
}
