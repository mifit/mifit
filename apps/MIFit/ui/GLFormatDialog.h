#ifndef GLFORMATDIALOG_H
#define GLFORMATDIALOG_H

#include <QtGui/QDialog>
#include <QtOpenGL/QGLFormat>

namespace Ui {
  class GLFormatDialog;
}

/**
 * Displays information about the OpenGL format for a given current format,
 * the default format, and the default overlay format. When the default format,
 * is displayed, a button provides access to change the default (using GLFormatEdit).
 *
 * The information is display as formated text to easily allow users to provide the
 * information by copy/paste.
 *
 */
class GLFormatDialog : public QDialog {
  Q_OBJECT
public:
  GLFormatDialog(QWidget *parent = 0);
  ~GLFormatDialog();

public slots:
  void setCurrentFormat(const QGLFormat& format);

private:
  Ui::GLFormatDialog *m_ui;
  QGLFormat currentFormat;

private slots:
  void formatTypeChanged(const QString& type);
  void changeDefault();

};

#endif // GLFORMATDIALOG_H
