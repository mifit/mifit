#include <cmath>
#include <math/Point3.h>
#include <opengl/Camera.h>
#include <opengl/Frustum.h>
#include <opengl/OpenGL.h>
#include <opengl/Viewport.h>

#include "core/corelib.h"
#include <chemlib/RESIDUE_.h>

#include "GLOverviewCanvas.h"
#include "Displaylist.h"
#include "MIMainWindow.h"
#include "MIGLWidget.h"

#include <QMouseEvent>
#include <QKeyEvent>

using namespace chemlib;
using namespace mi::math;
using namespace mi::opengl;
using namespace std;

class MIGLWidget;
static MIGLWidget *GetView()
{
    return MIMainWindow::instance()->currentMIGLWidget();
}


static void findCenter(std::vector<Bond> &bonds, Vector3<float> &min, Vector3<float> &max, Vector3<float> &center)
{
    min.set((float) INT_MAX, (float) INT_MAX, (float) INT_MAX);
    max.set((float) INT_MIN, (float) INT_MIN, (float) INT_MIN);
    std::vector<Bond>::iterator bondIter = bonds.begin();
    while (bondIter != bonds.end())
    {
        MIAtom *a1 = (*bondIter).getAtom1();
        MIAtom *a2 = (*bondIter).getAtom2();
        bondIter++;
        if (a1 != NULL)
        {
            min.x = std::min(min.x, a1->x());
            min.y = std::min(min.y, a1->y());
            min.z = std::min(min.z, a1->z());
            max.x = std::max(max.x, a1->x());
            max.y = std::max(max.y, a1->y());
            max.z = std::max(max.z, a1->z());
        }
        if (a2 != NULL)
        {
            min.x = std::min(min.x, a2->x());
            min.y = std::min(min.y, a2->y());
            min.z = std::min(min.z, a2->z());
            max.x = std::max(max.x, a2->x());
            max.y = std::max(max.y, a2->y());
            max.z = std::max(max.z, a2->z());
        }
    }
    center.set(max);
    center.add(min);
    center.scale(0.5f);
}

GLOverviewCanvas::GLOverviewCanvas(QWidget *parent)
    : QGLWidget(parent)
{

    viewport = new Viewport();
    frustum = new Frustum();
    camera = new Camera();

    viewpointCenterMarkerRadius = 3.0f;
    selectedResidueMarkerRadius = 2.0f;
    selectedResidueMarkerResolution = 8.0f;
    //initializeForRendering();
}

GLOverviewCanvas::~GLOverviewCanvas()
{
}


void GLOverviewCanvas::initializeGL()
{
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void GLOverviewCanvas::resizeGL(int width, int height)
{
    viewport->set(0, 0, width, height);
    //render();
}

void GLOverviewCanvas::setIdentity(float rotation[])
{
    rotation[0] = 1.0f;
    rotation[1] = 0.0f;
    rotation[2] = 0.0f;
    rotation[3] = 0.0f;

    rotation[4] = 0.0f;
    rotation[5] = 1.0f;
    rotation[6] = 0.0f;
    rotation[7] = 0.0f;

    rotation[8] = 0.0f;
    rotation[9] = 0.0f;
    rotation[10] = 1.0f;
    rotation[11] = 0.0f;

    rotation[12] = 0.0f;
    rotation[13] = 0.0f;
    rotation[14] = 0.0f;
    rotation[15] = 1.0f;
}

void GLOverviewCanvas::createAlphaCarbonTrace()
{
    if (!GetView())
        return;

    std::vector<chemlib::Bond>().swap(bonds); // was bonds.clear();
    Displaylist *displaylist = GetView()->GetDisplaylist();
    if (displaylist == NULL)
    {
        return;
    }
    list<Molecule*>::iterator molIter = displaylist->begin();
    list<Molecule*>::iterator molIterEnd = displaylist->end();
    while (molIter != molIterEnd)
    {
        BuildCALinks(bonds, (*molIter)->getResidues());
        molIter++;
    }

    findCenter(bonds, min, max, center);

    Point3<float> minPos(min);
    Point3<float> maxPos(max);
    float maxDistance = minPos.distance(maxPos) / 2.0f;

    float fieldOfView = 40.0f;
    double cameraDistance = maxDistance * 2.0 / (2.0 *  tan(toRadians(fieldOfView) * 0.5));
    cameraOffset.set(0.0f, 0.0f, -cameraDistance);
    frustum->setFieldOfView(fieldOfView);
    frustum->setPerspective(false);
    frustum->setFocalLength(cameraDistance);
    frustum->setNearClipping(0.01f);
    frustum->setFarClipping(3.0f * cameraDistance);
}

/**
 * Adjusts color for depth.
 */
void GLOverviewCanvas::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (GetView() == NULL)
    {
        return;
    }

    createAlphaCarbonTrace();

    ViewPoint *viewpoint = GetView()->viewpoint;

    float viewmat[3][3];
    viewpoint->copymatrix(viewmat);

    Vector3<float> target(viewpoint->getcenter(0), viewpoint->getcenter(1), viewpoint->getcenter(2));
    Vector3<float> eye(target);
    eye.add(cameraOffset);
    camera->setEye(eye);
    camera->lookAt(target);

    float fogColor[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    glEnable(GL_FOG);
    glFogfv(GL_FOG_COLOR, fogColor);
    glFogf(GL_FOG_START, 0.5 * frustum->getFocalLength());
    glFogf(GL_FOG_END, 1.2 * frustum->getFocalLength());
    glFogi(GL_FOG_MODE, GL_LINEAR);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    viewport->render();
    frustum->render(*viewport);

    camera->render();

    glTranslatef(viewpoint->getcenter(0), viewpoint->getcenter(1), viewpoint->getcenter(2));

    float rotation[16];
    setIdentity(rotation);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            rotation[i + j*4] = viewmat[i][j];
        }
    }
    glMultMatrixf(rotation);

    glPushMatrix();

    glTranslatef(-viewpoint->getcenter(0), -viewpoint->getcenter(1), -viewpoint->getcenter(2));

    const RESIDUE *res = GetView()->getFocusResidue();
    MIAtom *focusedResidueAlphaCarbon = NULL;
    if (res != NULL)
    {
        focusedResidueAlphaCarbon = atom_from_name("CA", *res);
    }

    vector<Bond>::iterator bondIter = bonds.begin();
    while (bondIter != bonds.end())
    {
        MIAtom *a1 = (*bondIter).getAtom1();
        MIAtom *a2 = (*bondIter).getAtom2();
        int c1 = ::getColor(a1);
        glColor3fv(getColor(c1));
        glBegin(GL_LINES);
        glVertex3f(a1->x(), a1->y(), a1->z());
        glVertex3f(a2->x(), a2->y(), a2->z());
        glEnd();
        // mark the focusres CA
        if (a1 == focusedResidueAlphaCarbon)
        {
            glColor3f(1.0f, 1.0f, 1.0f);
            drawSelectedResidueMarker(a1->x(), a1->y(), a1->z());
        }
        else if (a2 == focusedResidueAlphaCarbon)
        {
            glColor3f(1.0f, 1.0f, 1.0f);
            drawSelectedResidueMarker(a2->x(), a2->y(), a2->z());
        }
        bondIter++;
    }

    // Draw bounding box of view
    glTranslatef(viewpoint->getcenter(0), viewpoint->getcenter(1), viewpoint->getcenter(2));

    // Apply inverse rotation to keep elements aligned with camera.
    setIdentity(rotation);
    for (int i2 = 0; i2 < 3; ++i2)
    {
        for (int j2 = 0; j2 < 3; ++j2)
        {
            rotation[i2 + j2*4] = viewmat[j2][i2];
        }
    }
    glMultMatrixf(rotation);

    // Disable depth test to render on top of other elements.
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_FOG);

    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    glVertex3f(0.0f, -viewpointCenterMarkerRadius, 0.0f);
    glVertex3f(0.0f, viewpointCenterMarkerRadius, 0.0f);
    glVertex3f(-viewpointCenterMarkerRadius, 0.0f, 0.0f);
    glVertex3f(viewpointCenterMarkerRadius, 0.0f, 0.0f);
    glEnd();

    float viewWidth = viewpoint->getwidth()/2.0f;
    float viewHeight = viewpoint->getheight()/2.0f;
    glBegin(GL_LINES);
    glVertex2f(-viewWidth, -viewHeight);
    glVertex2f(viewWidth, -viewHeight);
    glVertex2f(-viewWidth, viewHeight);
    glVertex2f(viewWidth, viewHeight);
    glVertex2f(-viewWidth, -viewHeight);
    glVertex2f(-viewWidth, viewHeight);
    glVertex2f(viewWidth, -viewHeight);
    glVertex2f(viewWidth, viewHeight);
    glEnd();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_FOG);

    glPopMatrix();

    glFlush();
}

float*GLOverviewCanvas::getColor(int colorIndex)
{
    static float color[3];
    int c = PaletteIndex(colorIndex);
    color[0] = Colors::RPallette[c] / 255.0f;
    color[1] = Colors::GPallette[c] / 255.0f;
    color[2] = Colors::BPallette[c] / 255.0f;
    return color;
}

void GLOverviewCanvas::drawSelectedResidueMarker(float x, float y, float z)
{
    static GLUquadricObj *quadObj;
    if (!quadObj)
    {
        quadObj = gluNewQuadric();
    }
    glPushMatrix();
    glTranslatef(x, y, z);
    gluQuadricDrawStyle(quadObj, GLU_FILL);
    gluQuadricNormals(quadObj, GLU_SMOOTH);
    gluSphere(quadObj, selectedResidueMarkerRadius, selectedResidueMarkerResolution, selectedResidueMarkerResolution);
    glPopMatrix();
}

void GLOverviewCanvas::mousePressEvent(QMouseEvent *e)
{
    e->accept();
    if (!GetView())
        return;
    GetView()->mousePressEvent(e);
    if (GetView()->ViewChanged())
        updateGL();
}
void GLOverviewCanvas::mouseReleaseEvent(QMouseEvent *e)
{
    e->accept();
    if (!GetView())
        return;
    GetView()->mouseReleaseEvent(e);
    if (GetView()->ViewChanged())
        updateGL();
}
void GLOverviewCanvas::mouseMoveEvent(QMouseEvent *e)
{
    e->accept();
    if (!GetView())
        return;
    GetView()->mouseMoveEvent(e);
    if (GetView()->ViewChanged())
        updateGL();
}
void GLOverviewCanvas::mouseDoubleClickEvent(QMouseEvent *e)
{
    e->accept();
    if (!GetView())
        return;
    GetView()->mouseDoubleClickEvent(e);
    if (GetView()->ViewChanged())
        updateGL();
}

void GLOverviewCanvas::keyPressEvent(QKeyEvent *e)
{
    e->accept();

    unsigned short nRepCnt = 1;
    unsigned short nFlags = 0;

    unsigned short nChar;
    QString s = e->text();
    if (s.size()==0)
    {
        return;
    }
    nChar = s[0].toAscii();

    bool shift = (e->modifiers() & Qt::ShiftModifier);
    if (nChar <= 256 && isupper(nChar) && !shift)
    {
        nChar = tolower(nChar);
    }

    if (nChar == 186 && !shift)
    {
        nChar = ';';
    }
    if (nChar == 186 && shift)
    {
        nChar = ':';
    }
    if (nChar == 188 && !shift)
    {
        nChar = ',';
    }
    if (nChar == 188 && shift)
    {
        nChar = '<';
    }
    if (nChar == 190 && !shift)
    {
        nChar = '.';
    }
    if (nChar == 190 && shift)
    {
        nChar = '>';
    }
    if (nChar == 219 && !shift)
    {
        nChar = '[';
    }
    if (nChar == 219 && shift)
    {
        nChar = '{';
    }
    if (nChar == 221 && !shift)
    {
        nChar = ']';
    }
    if (nChar == 221 && shift)
    {
        nChar = '}';
    }

    GetView()->OnKeyDown(nChar, nRepCnt, nFlags);

    if (GetView()->ViewChanged())
    {
        updateGL();
    }
}

QSize GLOverviewCanvas::sizeHint() const
{
    return QSize(200, 200);
}
