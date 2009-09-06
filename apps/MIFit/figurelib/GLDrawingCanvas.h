#include <vector>
#include "Drawing.h"


#include "ui/MIEventHandler.h"
#include <QGLWidget>
#include <QMouseEvent>

class MIMenuBar;

namespace moldraw {


class GLDrawingCanvas : public QGLWidget, public Drawing, public MIEventHandler
{
  Q_OBJECT
  public:
    GLDrawingCanvas(QWidget *parent);
    virtual ~GLDrawingCanvas();

    virtual void FitToPage(std::vector<float> origin, std::vector<float> dimensions);
    virtual void Finish();

    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);

    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();

  private slots:
    void OnEditCopy();
    void OnPrint();
    void OnExportImage();
    void OnSelectFont();
    void OnAddAnnotation();
    void OnChangeFonts(const MIActionEvent &evt);

  private:
    void myRenderText(int x, int y, const QString & str,
                      const QFont & fnt = QFont(), int listBase = 2000);
    void myRenderText(double x, double y, double z, const QString & str,
                      const QFont & fnt = QFont(), int listBase = 2000);

    void toWorldCoords(const QPoint &pos, float &x, float &y);
    float orig_x, orig_y, dim_x, dim_y;
    unsigned int _selected_item;
    MIMenuBar *menuBar;
    std::string _current_font;
};

}
