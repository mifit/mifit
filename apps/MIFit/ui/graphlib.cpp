#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "graphlib.h"
#include "asplib.h"
#include "Application.h"

#include <QFont>
#include <QFontMetrics>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QSettings>

#define LABEL_ROOM (4*fheight)
#define TITLE_ROOM (3*fheight)

GraphWindow::GraphWindow(QWidget *parent)
    : QGLWidget(parent)
{
    /* initialize graph data */
    xlabel = "X";
    ylabel = "Y";
    title = "Graph";
    xlabelsize = 10;
    ylabelsize = 10;
    titlesize = 14;
    xdivisions = 5;
    ydivisions = 5;
    xminor = 5;
    yminor = 5;
    minx = asplib::GR_AUTO;
    maxx = asplib::GR_AUTO;
    miny = asplib::GR_AUTO;
    maxy = asplib::GR_AUTO;
    graphstyle = 0;
    bigFont = QFont();
    bigFont.setWeight(QFont::Light);
    rubber_x1 = rubber_x2 = rubber_y1 = rubber_y2 = 0.0f;
    _last_close_id = -1;

    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);
}

int Twidth(QFont &font, const std::string &c)
{
    QFontMetrics fm(font);
    return fm.width(c.c_str());
}

int Theight(QFont &font)
{
    QFontMetrics fm(font);
    return fm.height();
}



void GraphWindow::keyPressEvent(QKeyEvent *key_evt)
{
    emit keyPress(key_evt->key(),
                  (key_evt->modifiers() & (Qt::ControlModifier | Qt::ShiftModifier | Qt::AltModifier)));
}


bool GraphWindow::mouseHelper(QMouseEvent *e, float &fx, float &fy, float &xrange, float &yrange)
{

    QPoint pt = e->pos();
    int x = width(), y = height();
    int xorigin, yorigin, xend, yend;
    int fheight = Theight(bigFont);

    /*  add border */
    xend = x-20;
    /* this leaves LABEL_ROOM for labels below graph */
    yorigin = y-LABEL_ROOM-2*fheight;
    xorigin = 20 + Twidth(bigFont, ylabel);
    yend = TITLE_ROOM;

    xrange = (maxx-minx);
    yrange = (maxy-miny);

    fx = ((pt.x() - xorigin)*xrange/(xend-xorigin)) + minx;
    fy = ((pt.y() - yorigin)*yrange/(yend-yorigin)) + miny;

    return true;
}


void GraphWindow::finishMouseEvent(QMouseEvent *evt, float fx, float fy)
{
    std::vector<GR_POINT> close_pts;
    std::vector<unsigned int> close_ids;
    for (unsigned int i = 0; i < data.size(); ++i)
    {
        if (data[i].flag & (asplib::GR_SMALLBOX | asplib::GR_SMALLCIRCLE))
        {
            if (fabs(fx-data[i].x) < 5 && fabs(fy-data[i].y) < 5)
            {
                close_ids.push_back(i);
                close_pts.push_back(data[i]);
            }
        }
    }

    int closest = 0;
    float mindist = 999.99f;
    for (unsigned int i = 0; i < close_pts.size(); ++i)
    {
        float xsq = (fx-close_pts[i].x);
        xsq *= xsq;
        float ysq = (fy-close_pts[i].y);
        ysq *= ysq;
        float distsq = xsq+ysq;
        if (distsq < mindist)
        {
            closest = i;
            mindist = distsq;
        }
    }

    if (close_pts.size())
    {
        emit mouseOver(close_ids[closest]);
        if (evt->button() & Qt::LeftButton)
        {
            if (_last_close_id != -1 && _last_close_id < (int)data.size())
            {
                data[_last_close_id].flag &= ~asplib::GR_BOX;
            }
            _last_close_id = close_ids[closest];
            data[_last_close_id].flag |= asplib::GR_BOX;
                emit pick(close_pts[closest]);
        }
    }
}

void GraphWindow::mousePressEvent(QMouseEvent *evt)
{
    float xrange, yrange, fx, fy;
    if (!mouseHelper(evt, fx, fy, xrange, yrange))
    {
        return;
    }

    if (evt->button() & (Qt::LeftButton | Qt::RightButton))
    {
        if ((fx > maxx || fx < minx || fy > maxy || fy < miny))
        {
            graph_pop_scale();
            rubber_x1 = rubber_x2 = fx;
            rubber_y1 = rubber_y2 = fy;
            redraw();
            return;
        }
    }

    if (evt->button() & Qt::RightButton)
    {
        rubber_x1 = fx-fabs(fx*0.05);
        rubber_y1 = fy+fabs(fy*0.05);
        rubber_x2 = fx;
        rubber_y2 = fy;
        redraw();
        return;
    }

    finishMouseEvent(evt, fx, fy);
}

void GraphWindow::mouseReleaseEvent(QMouseEvent *evt)
{
    float xrange, yrange, fx, fy;
    if (!mouseHelper(evt, fx, fy, xrange, yrange))
    {
        return;
    }

    if (evt->button() & Qt::RightButton)
    {
        if (rubber_x2 < rubber_x1)
        {
            float tmp = rubber_x1;
            rubber_x1 = rubber_x2;
            rubber_x2 = tmp;
        }
        if (rubber_y2 < rubber_y1)
        {
            float tmp = rubber_y1;
            rubber_y1 = rubber_y2;
            rubber_y2 = tmp;
        }
        if (fabs(rubber_x1-rubber_x2) > xrange/20.0f  && fabs(rubber_y1-rubber_y2) > yrange/20.0f)
        {
            graph_push_scale(rubber_x1, rubber_x2, rubber_y1, rubber_y2);
        }
        rubber_x1 = rubber_x2 = rubber_y1 = rubber_y2 = 0.0f;
        redraw();
        return;
    }
    finishMouseEvent(evt, fx, fy);
}

void GraphWindow::mouseMoveEvent(QMouseEvent *evt)
{
    float xrange, yrange, fx, fy;
    if (!mouseHelper(evt, fx, fy, xrange, yrange))
    {
        return;
    }

    if (evt->buttons() & Qt::RightButton)
    {
        rubber_x2 = fx;
        rubber_y2 = fy;
        redraw();
        return;
    }
    finishMouseEvent(evt, fx, fy);
}

void GraphWindow::initializeGL()
{
    glClearColor(1.0, 1.0, 1.0, 1.0);

    // line rendering with the depth test enabled produces "interesting"
    // artifacts.  when two line segments start and end at the same
    // coordinates, the actual point of joining is not rendered, leading to a
    // "dotted" appearance.  You can see the same effect in the main GL
    // window, and I've seen this effect multiple times in other applications
    // glEnable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
}

void GraphWindow::resizeGL(int w, int h)
{
    glViewport(0, 0, (GLint)w, (GLint)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, (GLint)w, (GLint)h, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
}

void GraphWindow::redraw()
{
    updateGL();
}

void GraphWindow::graph_scale(double x1, double x2, double y1, double y2)
{
    minx = x1;
    miny = y1;
    maxx = x2;
    maxy = y2;
}

void GraphWindow::graph_push_scale(double x1, double x2, double y1, double y2)
{
    scale_stack.push_back(minx);
    scale_stack.push_back(miny);
    scale_stack.push_back(maxx);
    scale_stack.push_back(maxy);
    minx = x1;
    miny = y1;
    maxx = x2;
    maxy = y2;
}

void GraphWindow::graph_pop_scale()
{
    if (scale_stack.size() < 4)
    {
        return;
    }
    maxy = scale_stack.back();
    scale_stack.pop_back();
    maxx = scale_stack.back();
    scale_stack.pop_back();
    miny = scale_stack.back();
    scale_stack.pop_back();
    minx = scale_stack.back();
    scale_stack.pop_back();
}

int GraphWindow::graph_point(float x, float y, int flag, GraphColor color, int width)
{
    GR_POINT gr;
    gr.x = x;
    gr.y = y;
    gr.flag = flag;
    gr.color = color;
    gr.width = width;
    data.push_back(gr);
    return data.size()-1;
}

int GraphWindow::graph_point(float x, float y, int flag, void *user_data, GraphColor color, int width)
{
    GR_POINT gr;
    gr.x = x;
    gr.y = y;
    gr.flag = flag;
    gr.data = user_data;
    gr.color = color;
    gr.width = width;
    data.push_back(gr);
    return data.size()-1;
}

bool GraphWindow::graph_replacepoint(int i, float x, float y, int flag)
{
    if (data.size() <= (unsigned int)i)
    {
        return false;
    }
    data[i].x = x;
    data[i].y = y;
    data[i].flag = flag;
    return true;
}

bool GraphWindow::graph_labelpoint(int i, const char *string)
{
    if (data.size() <= (unsigned int)i)
    {
        return false;
    }
    if (!string)
    {
        data[i].label = "";
        data[i].flag &= ~asplib::GR_LABEL;
    }
    else
    {
        data[i].label = string;
        data[i].flag |= asplib::GR_LABEL;
    }
    return (true);
}

void GraphWindow::graph_removepointsafter(int i)
{
    data.resize(i);
}

void GraphWindow::graph_curvelabel(int /* i */, const char *label, int flag, GraphColor color)
{
    CURVELABEL l;
    //if(i >= curvelabel.size() ) curvelabel.reserve(i+1);
    l.label = label;
    l.flag = flag;
    l.color = color;
    curvelabel.push_back(l);
}


static void stringSplit(std::string str, std::string delim, std::vector<std::string> &results)
{
    unsigned int cutAt;
    while ( (cutAt = str.find_first_of(delim)) != str.npos)
    {
        if (cutAt > 0)
        {
            results.push_back(str.substr(0, cutAt));
        }
        str = str.substr(cutAt+1);
    }
    if (str.length() > 0)
    {
        results.push_back(str);
    }
}


//taken from Qt 4.2.3 to avoid assertion failure with faulty glx driver(s) from xorg 6.9/7.0 on linux.
void GraphWindow::myRenderText(int x, int y, const QString &str, const QFont &font, int listBase)
{
#ifdef __linux__
    static bool use_alternate_rendering = false;
    static bool firsttime = true;
    if (firsttime)
    {
        firsttime = false;
        use_alternate_rendering =
            QSettings().value("View Parameters/AlternateTextRendering", false).toBool();
    }

    if (use_alternate_rendering)
    {
        makeCurrent();
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_TEXTURE_2D);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, width(), height(), 0, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glAlphaFunc(GL_GREATER, 0.0);
        glEnable(GL_ALPHA_TEST);
        glRasterPos2i(0, 0);
        glBitmap(0, 0, 0, 0, x, -y, NULL);
        glListBase(fontDisplayListBase(font, listBase));
        QByteArray cstr(str.toLatin1());
        glCallLists(cstr.size(), GL_UNSIGNED_BYTE, cstr.constData());

        // can't support this, but not needed
        //if (font.underline() || font.strikeOut() || font.overline())
        // qt_drawFontLining(x, y, str, font);

        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopAttrib();
        return;
    }
#endif // ifdef __linux__
    renderText(x, y, str, font, listBase);
}

//taken from Qt 4.2.3 to avoid assertion failure with faulty glx driver(s) from xorg 6.9/7.0 on linux.
void GraphWindow::myRenderText(double x, double y, double z, const QString &str, const QFont &font,
                               int listBase)
{
#ifdef __linux__
    static bool use_alternate_rendering = false;
    static bool firsttime = true;
    if (firsttime)
    {
        firsttime = false;
        use_alternate_rendering =
            QSettings().value("View Parameters/AlternateTextRendering", false).toBool();
    }

    if (use_alternate_rendering)
    {
        makeCurrent();
        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glDisable(GL_TEXTURE_1D);
        glDisable(GL_TEXTURE_2D);
        glDisable(GL_CULL_FACE);

        glRasterPos3d(x, y, z);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glAlphaFunc(GL_GREATER, 0.0);
        glEnable(GL_ALPHA_TEST);
        glListBase(fontDisplayListBase(font, listBase));
        QByteArray cstr(str.toLatin1());
        glCallLists(cstr.size(), GL_UNSIGNED_BYTE, cstr.constData());

        glPopAttrib();
        return;
    }
#endif // ifdef __linux__
    renderText(x, y, z, str, font, listBase);
}



void GraphWindow::DrawMultiLineString(float x, float y, const std::string &c)
{
    std::vector<std::string> results;
    std::string delim("\n");
    std::string str(c.c_str());

    stringSplit(str, delim, results);

    int lineSpacing = Theight(bigFont);
    for (unsigned int i = 0; i < results.size(); ++i)
    {
        myRenderText(x, y, 0.0, results[i].c_str());
        y += lineSpacing;
    }
}

void drawCircle (int x, int y, int r)
{
    float cx, cy;
    glBegin(GL_POLYGON);
    for (int i = 0; i <= 360; i += 18)
    {
        cx = (float)r * cos(i*DEG2RAD);
        cy = (float)r * sin(i*DEG2RAD);
        glVertex3f(cx + x, cy + y, 0.0f);
    }
    glEnd();
}



#define String(a, b, c) { DrawMultiLineString((a), (b), (c)); }

#define Point(a, b) {glBegin(GL_POINTS); glVertex2f((a), (b)); glEnd(); }

#define Line(x1, y1, x2, y2) {glBegin(GL_LINES); glVertex2f((x1), (y1)); glVertex2f((x2), (y2)); glEnd(); }

#define DashLine(a, b, c, d) {glLineStipple(4, 0xAAAA); glEnable(GL_LINE_STIPPLE); Line((a), (b), (c), (d)); glDisable(GL_LINE_STIPPLE); }
#define DotLine(a, b, c, d)  {glLineStipple(1, 0xAAAA); glEnable(GL_LINE_STIPPLE); Line((a), (b), (c), (d)); glDisable(GL_LINE_STIPPLE); }

#define Label(a, b, c) {Line((a), (b)+8, (a), (b)+14); String((a), (b)+14, (c));}

#define Cross(a, b)      {Line((a)-5, (b), (a)+5, (b)); Line((a), (b)-5, (a), (b)+5);}
#define SmallCross(a, b) {Line((a)-3, (b), (a)+3, (b)); Line((a), (b)-3, (a), (b)+3);}

#define Box(a, b)      {glRectf((a)-4, (b)-4, (a)+4, (b)+4);}
#define SmallBox(a, b) {glRectf((a)-2, (b)-2, (a)+2, (b)+2); }

#define Circle(a, b)      {drawCircle((a), (b), 5);}
#define SmallCircle(a, b) {drawCircle((a), (b), 3);}


void GraphWindow::paintGL()
{
    if (!isVisible())
        return;

    char buff[200];
    int xorigin, yorigin, xend, yend;
    float xint, yint, x1, y1, x2, y2;
    float xmin = FLT_MAX, xmax = -FLT_MAX, ymin = FLT_MAX, ymax = -FLT_MAX;
    int i;
    float res;
    int ndata = data.size();
    int fheight = Theight(bigFont);
    float ydiv, xdiv;
    int x = width(), y = height();

    resizeGL(x, y); // should be unnecessary, but isn't.  something's screwing up the matrix and this resets it.

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glColor3f(0.0f, 0.0f, 0.0f);
    glLineWidth(1.0);

    /*  add border */
    xend = x-20;
    /* this leaves LABEL_ROOM for labels below graph */
    yorigin = y-LABEL_ROOM-2*fheight;
    xorigin = 20 + Twidth(bigFont, ylabel);
    yend = TITLE_ROOM;

    String(xorigin, yend-TITLE_ROOM/2, title);

    /* draw axes */
    Line(xorigin, yorigin, xorigin, yend);
    Line(xorigin, yorigin, xend, yorigin);
    Line(xend, yorigin, xend, yend);
    Line(xorigin, yend, xend, yend);
    /* label axes */
    x1 = xorigin-10-Twidth(bigFont, ylabel);
    /* if even ydivisions then label is moved to avoid collision */
    y1 = yorigin/2+yend/2 - (ydivisions+1)%2*fheight;
    String(x1, y1, ylabel);
    x1 = xorigin/2.+xend/2.-Twidth(bigFont, xlabel)/2.;
    y1 = yorigin+2*fheight;
    String(x1, y1, xlabel);
    if (ndata > 0)
    {
        for (i = 0; i < ndata; i++)
        {
            if (xmax < data[i].x)
            {
                xmax = data[i].x;
            }
            if (ymax < data[i].y)
            {
                ymax = data[i].y;
            }
            if (xmin > data[i].x)
            {
                xmin = data[i].x;
            }
            if (ymin > data[i].y)
            {
                ymin = data[i].y;
            }
        }
        if (maxx != asplib::GR_AUTO)
        {
            xmax = maxx;
        }
        if (maxy != asplib::GR_AUTO)
        {
            ymax = maxy;
        }
        if (minx != asplib::GR_AUTO)
        {
            xmin = minx;
        }
        if (miny != asplib::GR_AUTO)
        {
            ymin = miny;
        }
        if (xmax != xmin) /* i.e. if there is only 1 data point */
        {
            xint = (xend-xorigin)/(xmax-xmin); /* interval for 1 in x */
        }
        else
        {
            xint = 1.0;
        }
        if (ymax != ymin)
        {
            yint = (yend-yorigin)/(ymax-ymin);
        }
        else
        {
            yint = 1.0;
        }
        ydiv = ydivisions;
        xdiv = xdivisions;
        for (i = 0; i <= xdivisions; i++)
        {
            if (!(graphstyle&asplib::GR_CUBEDTORES))
            {
                sprintf(buff, "%0.3g", xmin+(float)i*(xmax-xmin)/xdiv);
            }
            else
            {
                res = xmin+(float)i*(xmax-xmin)/xdiv;
                res = 0.5/ pow(res, 1.0f/3.0f);
                sprintf(buff, "%0.3g", res);
            }
            x1 = xorigin+i*(xend-xorigin)/xdiv;
            y1 = yorigin + fheight;
            String(x1-Twidth(bigFont, buff)/2.0, y1, buff);
            x2 = x1;
            y1 = yorigin;
            y2 = yend;
            DashLine(x1, y1, x2, y2);
        }
        for (i = 0; i <= ydivisions; i++)
        {
            sprintf(buff, "%0.3g", ymin+(float)i*(ymax-ymin)/ydiv);
            y1 = yorigin+i*(yend-yorigin)/ydiv;
            x1 = xorigin-2-Twidth(bigFont, buff);
            String(x1, y1  + fheight/2.0, buff);
            y2 = y1;
            x1 = xorigin;
            x2 = xend;
            DashLine(x1, y1, x2, y2);
        }
        /* minor x marks */
        ydiv = ydivisions*yminor;
        xdiv = xdivisions*xminor;
        for (i = 0; i < xdiv; i++)
        {
            x1 = x2 = xorigin+i*(xend-xorigin)/xdiv;
            y1 = yorigin;
            y2 = yorigin-4;
            Line((int)x1, (int)y1, (int)x2, (int)y2);
            y1 = yend;
            y2 = yend+4;
            Line((int)x1, (int)y1, (int)x2, (int)y2);
        }
        for (i = 0; i < ydiv; i++)
        {
            y1 = y2 = yorigin+i*(yend-yorigin)/ydiv;
            x1 = xorigin;
            x2 = xorigin+4;
            Line((int)x1, (int)y1, (int)x2, (int)y2);
            x1 = xend;
            x2 = xend-4;
            Line((int)x1, (int)y1, (int)x2, (int)y2);
        }

        if (ndata > 1)
        {
            for (i = 0; i < ndata-1; i++)
            {
                glColor3ub(data[i].color._r, data[i].color._g, data[i].color._b);
                glLineWidth((GLfloat)data[i].width);

                // implement clipping
                if (data[i].x < minx || data[i].x > maxx
                    || data[i].y < miny || data[i].y > maxy)
                {
                    continue;
                }


                x1 = xorigin + (data[i].x - xmin)*xint;
                y1 = yorigin + (data[i].y - ymin)*yint;
                x2 = xorigin + (data[i+1].x - xmin)*xint;
                y2 = yorigin + (data[i+1].y - ymin)*yint;
                if (!(data[i].flag & asplib::GR_BREAK))
                {
                    if (data[i].flag&asplib::GR_DASH)
                    {
                        DashLine(x1, y1, x2, y2);
                    }
                    else
                    {
                        Line((int)x1, (int)y1, (int)x2, (int)y2);
                    }
                }
                if (data[i].flag & asplib::GR_CROSS)
                {
                    Cross(x1, y1);
                }
                if (data[i].flag & asplib::GR_BOX)
                {
                    Box(x1, y1);
                }
                if (data[i].flag & asplib::GR_CIRCLE)
                {
                    Circle(x1, y1);
                }
                if (data[i].flag & asplib::GR_SMALLCROSS)
                {
                    SmallCross(x1, y1);
                }
                if (data[i].flag & asplib::GR_SMALLBOX)
                {
                    SmallBox(x1, y1);
                }
                if (data[i].flag & asplib::GR_SMALLCIRCLE)
                {
                    SmallCircle(x1, y1);
                }
                if (data[i].flag & asplib::GR_POINT)
                {
                    Point(x1, y1);
                }
                if (data[i].flag & asplib::GR_LABEL)
                {
                    Label(x1, y1, data[i].label);
                }
            }
        }
        /* do last point's cross */
        i = ndata -1;
        x1 = xorigin + (data[i].x - xmin)*xint;
        y1 = yorigin + (data[i].y - ymin)*yint;
        if (data[i].flag & asplib::GR_CROSS)
        {
            Cross(x1, y1);
        }
        if (data[i].flag & asplib::GR_BOX)
        {
            Box(x1, y1);
        }
        if (data[i].flag & asplib::GR_CIRCLE)
        {
            Circle(x1, y1);
        }
        if (data[i].flag & asplib::GR_SMALLCIRCLE)
        {
            SmallCircle(x1, y1);
        }
        if (data[i].flag & asplib::GR_SMALLCROSS)
        {
            SmallCross(x1, y1);
        }
        if (data[i].flag & asplib::GR_SMALLBOX)
        {
            SmallBox(x1, y1);
        }
        if (data[i].flag & asplib::GR_POINT)
        {
            Point(x1, y1);
        }
        if (data[i].flag & asplib::GR_LABEL)
        {
            Label(x1, y1, data[i].label);
        }

        glColor3f(0.0f, 0.0f, 0.0f);

        /* draw any curve labels */
        x1 = xorigin;
        y1 = yorigin+2.5*fheight;

        for (i = 0; (unsigned int)i < curvelabel.size(); i++)
        {
            glColor3ub(curvelabel[i].color._r, curvelabel[i].color._g, curvelabel[i].color._b);

            float x3, y3;
            if (curvelabel[i].label.c_str()[0] != 0)
            {
                //		if(!(curvelabel[i].flag&asplib::GR_BREAK)){
                //			if(curvelabel[i].flag&asplib::GR_DASH){
                //				DashLine(x1,y1,x2,y2);
                //			}else {Line((int)x1,(int)y1,(int)x2,(int)y2);}}
                if (curvelabel[i].flag & asplib::GR_CROSS)
                {
                    Cross(x1, y1);
                }
                if (curvelabel[i].flag & asplib::GR_BOX)
                {
                    Box(x1, y1);
                }
                if (curvelabel[i].flag & asplib::GR_CIRCLE)
                {
                    Circle(x1, y1);
                }
                if (curvelabel[i].flag & asplib::GR_SMALLCROSS)
                {
                    SmallCross(x1, y1);
                }
                if (curvelabel[i].flag & asplib::GR_SMALLBOX)
                {
                    SmallBox(x1, y1);
                }
                if (curvelabel[i].flag & asplib::GR_SMALLCIRCLE)
                {
                    SmallCircle(x1, y1);
                }
                if (curvelabel[i].flag & asplib::GR_POINT)
                {
                    Point(x1, y1);
                }
                x3 = x1+8;
                y3 = y1+fheight/4.0;
                String(x3, y3, curvelabel[i].label);
                y1 += fheight*3/2;
                if (y1 > y-fheight)
                {
                    x1 += (xend-xorigin)/2.0;
                    y1 = yorigin+2.5*fheight;
                }
            }
        }


        // rubber band
        if (rubber_x1 != rubber_x2 && rubber_y1 != rubber_y2)
        {
            glColor3i(128, 128, 128);

            x1 = xorigin + (rubber_x1 - xmin)*xint;
            y1 = yorigin + (rubber_y1 - ymin)*yint;
            x2 = xorigin + (rubber_x2 - xmin)*xint;
            y2 = yorigin + (rubber_y2 - ymin)*yint;
            glLineStipple(4, 0xAAAA);
            glEnable(GL_LINE_STIPPLE);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glRectf(x1, y1, x2, y2);
            glDisable(GL_LINE_STIPPLE);
        }
    }

    glPopAttrib();
}
