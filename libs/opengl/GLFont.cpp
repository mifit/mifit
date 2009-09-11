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
  glWidget_->renderText(0.0, 0.0, 0.0, text, font_);
}

}
}
