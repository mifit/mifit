#ifndef GLFORMATEDIT_H
#define GLFORMATEDIT_H

#include <QtGui/QDialog>
#include <QtOpenGL/QGLFormat>

namespace Ui {
    class GLFormatEdit;
}
class QAbstractButton;

/**
 * Editor dialog for changing an OpenGL format.
 *
 * Used to allow user to change the default OpenGL format before
 * a main window is created. The default selected by Qt may not
 * work properly on some computers.
 */
class GLFormatEdit : public QDialog {
    Q_OBJECT
public:
    GLFormatEdit(QWidget *parent = 0);
    ~GLFormatEdit();

    QGLFormat format() const;

public slots:
    void setFormat(const QGLFormat& format);
    void reset();

private:
    Ui::GLFormatEdit *m_ui;
    QGLFormat format_;

private slots:
    void buttonClicked(QAbstractButton* button);
};

#endif // GLFORMATEDIT_H
