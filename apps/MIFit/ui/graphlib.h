#ifndef graph_HEADER
#define graph_HEADER

#include <string>
#include <QFont>
#include <QGLWidget>
#include <vector>


class GraphColor
{
public:
    GraphColor() : _r(0), _g(0), _b(0)
    {
    }

    GraphColor(unsigned char r, unsigned char g, unsigned char b)
        : _r(r), _g(g), _b(b)
    {
    }

    unsigned char _r, _g, _b;
};


class CURVELABEL
{
public:
    std::string label;
    int flag;
    GraphColor color;
};

class GR_POINT
{
public:
    GR_POINT() : data(0), color(GraphColor(255, 255, 255)), width(1)
    {
    }

    float x, y;
    int flag;
    std::string label;
    void *data;
    GraphColor color;
    int width;
};

class GraphWindow : public QGLWidget
{
    Q_OBJECT

public:
    GraphWindow(QWidget *parent);
    ~GraphWindow()
    {
    }

    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void redraw();

    void keyPressEvent(QKeyEvent *e);

    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);


    const std::vector<GR_POINT>&GetData()
    {
        return data;
    }

    void graph_push_scale(double x1, double x2, double y1, double y2);
    void graph_pop_scale();

signals:
    void keyPress(int keycode, bool shift);
    void pick(const GR_POINT &point);
    void mouseOver(int id);

private:
    void myRenderText(int x, int y, const QString &str,
                      const QFont &fnt = QFont(), int listBase = 2000);
    void myRenderText(double x, double y, double z, const QString &str,
                      const QFont &fnt = QFont(), int listBase = 2000);

    bool mouseHelper(QMouseEvent *e, float &fx, float &fy, float &xrange, float &yrange);
    void finishMouseEvent(QMouseEvent *e, float fx, float fy);
    void DrawMultiLineString(float x, float y, const std::string &c);

    std::vector<GR_POINT> data;
    std::string xlabel;
    int xlabelsize;     /* not currently implemented */
    std::string ylabel;
    int ylabelsize;     /* not currently implemented */
    std::string title;
    int titlesize;     /* not currently implemented */
    int graphstyle;
    std::vector<CURVELABEL> curvelabel;
    double maxx;     /*axis scaling - use GR_AUTO for autoscaling*/
    double maxy;
    double minx;
    double miny;
    int xdivisions, ydivisions;
    int xminor, yminor;
    QFont bigFont;

    float rubber_x1, rubber_x2, rubber_y1, rubber_y2;
    std::vector<double> scale_stack;
    int _last_close_id;

public:
    void graph_scale(double x1, double x2, double y1, double y2);
    int graph_point(float x, float y, int flag, GraphColor color = GraphColor(0, 0, 0), int width = 1);
    int graph_point(float x, float y, int flag, void *data, GraphColor color = GraphColor(0, 0, 0), int width = 1);
    bool graph_replacepoint(int i, float x, float y, int flag);
    bool graph_labelpoint(int i, const char *string);
    void graph_removepointsafter(int i);
    void graph_clear()
    {
        std::vector<GR_POINT>().swap(data); // was data.clear();
        std::vector<CURVELABEL>().swap(curvelabel); // was curvelabel.clear();
        _last_close_id = -1;
    }

    void graph_xlabel(const char *s)
    {
        xlabel = s;
    }

    void graph_ylabel(const char *s)
    {
        ylabel = s;
    }

    void graph_title(const char *s)
    {
        title = s;
    }

    void graph_xdivisions(int div)
    {
        if (div < 0)
        {
            div = 0;
        }
        xdivisions = div;
    }

    void graph_ydivisions(int div)
    {
        if (div < 0)
        {
            div = 0;
        }
        ydivisions = div;
    }

    void graph_xminor(int div)
    {
        if (div < 0)
        {
            div = 0;
        }
        xminor = div;
    }

    void graph_yminor(int div)
    {
        if (div < 0)
        {
            div = 0;
        }
        yminor = div;
    }

    void graph_set_style(int style)
    {
        graphstyle = style;
    }

    void graph_curvelabel(int i, const char *label, int flag, GraphColor color = GraphColor(0, 0, 0));
private:
};


#endif // ifndef graph_HEADER
