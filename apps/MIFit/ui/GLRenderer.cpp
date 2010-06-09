#include <math/Matrix4.h>
#include <math/Vector3.h>
#include <math/Vector4.h>
#include <opengl/Camera.h>
#include <opengl/Frustum.h>
#include <opengl/Light.h>
#include <opengl/OpenGL.h>
#include <opengl/GLFont.h>

#include <nongui/nonguilib.h>
#include <map/maplib.h>
#include "core/corelib.h"
#include <chemlib/Monomer.h>
#include "GLRenderer.h"
#include "TargaImage.h"
#include "EMap.h"
#include <climits>
#include <vector>
#undef X
#undef Y
#undef Z

#include "Application.h"
#include "Displaylist.h"

#include <QGLWidget>
#include <QPainter>


using namespace chemlib;
using namespace mi::math;
using namespace mi::opengl;
using namespace std;

static const bool useClockwiseFrontFaces = false;
static float white[4] = { 1.0f, 1.0f, 1.0f, 1.0f };

static void logOpenGLErrors(const char *file, int line)
{
    int error = glGetError();
    if (error != 0)
    {
        Logger::log("OpenGL error %s (%d) at line %d in %s", gluErrorString(error), error, line, file);
    }
}

GLRenderer::GLRenderer()
    :  pickingEnabled(false),
      qcanvas_win(NULL),
      popStackImage(NULL),
      minimizeStackImage(NULL),
      unminimizeStackImage(NULL),
      closeStackImage(NULL),
      hideHydrogens(false),
      light(NULL),
      light2(NULL),
      frustum(NULL),
      camera(NULL),
      currentModel(true),
      fogEnabled(false),
      antialiasLines(true),
      joinBondsOfSameColor(true),
      showBondOrders(false),
      fontSize(ATOMLABEL::defaultSize()),
      viewVectorSet(false)
{

    popStackImage = new TargaImage(Application::instance()->GetMolimageHome() + "/data/images/popStack.tga");
    minimizeStackImage = new TargaImage(Application::instance()->GetMolimageHome() + "/data/images/minimizeStack.tga");
    unminimizeStackImage = new TargaImage(Application::instance()->GetMolimageHome() + "/data/images/unminimizeStack.tga");
    closeStackImage = new TargaImage(Application::instance()->GetMolimageHome() + "/data/images/closeStack.tga");

    style.set(RenderStyle::getBallAndStick());

    light = new Light(GL_LIGHT0);
    light->setAmbientColor(0.6f, 0.6f, 0.6f, 1.0f);
    light->setSpecularColor(1.0f, 1.0f, 1.0f, 1.0f);
    Vector4<float> lightPosition(0.0f, 20.0f, 10.0f, 0.0f);
    light->setPosition(lightPosition);

    light2 = new Light(GL_LIGHT1);
    light2->setAmbientColor(0.1f, 0.1f, 0.1f, 1.0f);
    light2->setDiffuseColor(0.0f, 0.0f, 0.0f, 1.0f);
    light2->setSpecularColor(0.0f, 0.0f, 0.0f, 1.0f);
    Vector4<float> lightPosition2(0.0f, 150.0f, 10.0f, 0.0f);
    light2->setPosition(lightPosition2);
}

GLRenderer::~GLRenderer()
{
    delete light;
    delete light2;
    delete popStackImage;
    delete minimizeStackImage;
    delete unminimizeStackImage;
    delete closeStackImage;
    for (std::map<void*, GLUquadric*>::iterator i = quad.begin();
         i!=quad.end(); ++i)
    {
        gluDeleteQuadric(i->second);
    }
}

void GLRenderer::setViewVector(const Vector3<float> &viewVector)
{
    this->viewVector.set(viewVector);
    viewVectorSet = true;
}

void GLRenderer::updateTextScale(float glUnitsPerPixel)
{
    fontSize = ATOMLABEL::defaultSize();
    double textHeight = (double) fontSize * glUnitsPerPixel;
    labelsTextScale = (float) fontSize / 16 * 0.1 * textHeight;
}

void GLRenderer::renderFog()
{
    if (fogEnabled)
    {
        float fogColor[] = { 0.2f, 0.2f, 0.2f, 1.0f };
        glEnable(GL_FOG);
        glFogfv(GL_FOG_COLOR, fogColor);
        glFogf(GL_FOG_START, fogStart);
        glFogf(GL_FOG_END, fogEnd);
        glFogi(GL_FOG_MODE, GL_LINEAR);
    }
    else
    {
        glDisable(GL_FOG);
    }
}

void GLRenderer::initializeForRender()
{

    glEnable(GL_DEPTH_TEST);
    if (useClockwiseFrontFaces)
    {
        glFrontFace(GL_CW);
    }
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
}

void GLRenderer::renderLight()
{
    light->glEnable();
    light->render();
    light2->glEnable();
    light2->render();
    light2->glDisable();
}


void GLRenderer::setQGLWidget(QGLWidget *qcanvas_win)
{
    this->qcanvas_win = qcanvas_win;
    labelsFont_ = createGLFont();
    stackFont_ = createGLFont();
#ifndef Q_OS_LINUX
    // Linux has problem rendering stack font if different than labels.
    stackFont_->setFont(QFont("Helvetica", 10));
#endif
}

void*GLRenderer::getContext()
{
    if (!qcanvas_win)
        return 0;
    return (void*)qcanvas_win->context();
}

void GLRenderer::updateFonts()
{
}

void GLRenderer::swapBuffers()
{
    if (qcanvas_win)
        qcanvas_win->swapBuffers();
}

void GLRenderer::computeBounds(std::list<Molecule*> &molecules)
{
    float minX = (float) INT_MAX;
    float minY = (float) INT_MAX;
    float minZ = (float) INT_MAX;
    float maxX = (float) INT_MIN;
    float maxY = (float) INT_MIN;
    float maxZ = (float) INT_MIN;

    std::list<Molecule*>::iterator molIter = molecules.begin();
    while (molIter != molecules.end())
    {
        Molecule *molecule = *molIter;
        molIter++;
        for (MIIter<Residue> res = molecule->GetResidues(); res; ++res)
        {
            for (int i = 0; i < res->atomCount(); i++)
            {
                MIAtom *atom = res->atom(i);
                minX = std::min(minX, atom->x());
                minY = std::min(minY, atom->y());
                minZ = std::min(minZ, atom->z());
                maxX = std::max(maxX, atom->x());
                maxY = std::max(maxY, atom->y());
                maxZ = std::max(maxZ, atom->z());
            }
        }
    }

    xyDisplayRange = std::max(std::max(maxX-minX, maxY-minY), maxZ-minZ);
    if (xyDisplayRange == 0.0f)
    {
        xyDisplayRange = 1.0f;
    }
}

void GLRenderer::applyProjection(float scale)
{
    float zDisplayRange = xyDisplayRange;

    // Maintain square perspective regardless of window aspect ratio
    float widthRange = xyDisplayRange / scale;
    float aspectRatio = (float)width / (float)height;
    if (width <= height)
    {
        glOrtho(-widthRange, widthRange, -widthRange / aspectRatio, widthRange / aspectRatio, 0.0, 2.0*zDisplayRange);
    }
    else
    {
        glOrtho(-widthRange * aspectRatio, widthRange * aspectRatio, -widthRange, widthRange, 0.0, 2.0*zDisplayRange);
    }
}

void GLRenderer::Draw2(Displaylist *displaylist, bool ShowVus, bool /* drawMapsFirst */, bool showSymmetryAsBackbone)
{

    logOpenGLErrors(__FILE__, __LINE__);
    computeBounds(displaylist->getMolecules());

    //float zDisplayRange = xyDisplayRange;

    float specularMaterial[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, specularMaterial);
    glMateriali(GL_FRONT, GL_SHININESS, 8);

    glShadeModel(GL_SMOOTH);

    glPushMatrix();
    logOpenGLErrors(__FILE__, __LINE__);

    DrawMaps(displaylist->getMaps());
    int moleculeCounter = 0;
    std::list<Molecule*>::iterator node = displaylist->getMolecules().begin();
    while (node != displaylist->getMolecules().end())
    {
        if (dimNonactiveModels)
        {
            if ((*node) == displaylist->CurrentItem())
            {
                currentModel = true;
                light->glEnable();
                light2->glDisable();
            }
            else
            {
                light->glDisable();
                light2->glEnable();
                currentModel = false;
            }
        }
        else
        {
            currentModel = true;
            light->glEnable();
            light2->glDisable();
        }
        drawMolecule(*node);
        drawSymmetryMolecule(*node, showSymmetryAsBackbone);

        moleculeCounter++;
        if (ShowVus)
        {
            drawLines(displaylist->getLines(), 1, true);
        }
        drawSurface(displaylist->getCurrentDots());
        if ((*node)->getSecondaryStructure() != NULL)
        {
            glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT);
            glDisable(GL_CULL_FACE);
            DrawSecondaryStructure((*node)->getSecondaryStructure());
            glPopAttrib();
        }

        node++;
    }
    logOpenGLErrors(__FILE__, __LINE__);

    glPopMatrix();

    logOpenGLErrors(__FILE__, __LINE__);
}

void GLRenderer::DrawSecondaryStructure(SecondaryStructure *secondaryStructure)
{
    RibbonSegment *pNextRSeg = secondaryStructure->m_pRibbonSegmentList;
    while (pNextRSeg)
    {
        DrawRibbonSegment(pNextRSeg);
        pNextRSeg = pNextRSeg->m_pNext;
    }

    std::vector<Helix*>::iterator pHelix = secondaryStructure->m_pHelixList.begin();
    for (; pHelix != secondaryStructure->m_pHelixList.end(); ++pHelix)
    {
        DrawHelix(*pHelix);
    }

    pNextRSeg = secondaryStructure->m_pSheetList;
    while (pNextRSeg)
    {
        DrawRibbonSegment(pNextRSeg);
        pNextRSeg = pNextRSeg->m_pNext;
    }


    pNextRSeg = secondaryStructure->m_pTurnList;
    while (pNextRSeg)
    {
        DrawRibbonSegment(pNextRSeg);
        pNextRSeg = pNextRSeg->m_pNext;
    }

    pNextRSeg = secondaryStructure->m_pRandomList;
    while (pNextRSeg)
    {
        DrawRibbonSegment(pNextRSeg);
        pNextRSeg = pNextRSeg->m_pNext;
    }

}

void GLRenderer::DrawRibbonSegment(RibbonSegment *ribbonSegment)
{
    SetColor(ribbonSegment->m_bRGB);
    RibbonSpan *pSpan = ribbonSegment->m_pRibbonSpanList->m_pNext;
    while (pSpan)
    {
        DrawRibbonSpan(pSpan);
        pSpan = pSpan->m_pNext;
    }
}

void GLRenderer::SetColor(unsigned char *bRGBA)
{
    //Set the new color
    //glColor3b (m_bRGBA[0], m_bRGBA[1], m_bRGBA[2]); //OpenGL uses signed bytes!
    float fColor[3];
    for (int i = 0; i < 3; i++)
    {
        fColor[i] = (float)bRGBA[i]/255.0;
    }
    glColor3f(fColor[0], fColor[1], fColor[2]);
}

void GLRenderer::DrawRibbonSpan(RibbonSpan *ribbonSpan)
{
    DrawPolySurf(ribbonSpan->npr, ribbonSpan->nseg, ribbonSpan->ribbon_pts, ribbonSpan->ribbon_norms);
}

void GLRenderer::DrawHelix(Helix *helix)
{
    SetColor(helix->rgb);
    DrawCappedCylinder(Helix_NSEG, helix->m_dAxisA, helix->m_dAxisB, helix->m_dEndA, helix->m_dEndB, helix->m_dNormA, helix->m_dNormB, helix->m_dCapNormA, helix->m_dCapNormB);
}

void GLRenderer::DrawPolySurf(int npr, int nseg, double dPoints[][3], double dNorms[][3])
{
    //Make sure we are in material mode
    //  if (m_iColorMode != SSI_MATERIAL)
    //    SetColorMode (SSI_MATERIAL);

    //glBegin (GL_QUADS);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < npr; i++) // Begin for each parametric step
    {
        for (int k = 0; k < nseg-1; k++) // Begin for each point in the ribbon
        {
            int ip0 = i*nseg+k;
            int ip1;
            if (i == npr-1)
            {
                ip1 = k;
            }
            else
            {
                ip1 = (i+1)*nseg+k;
            }
            int ip2 = ip1+1;
            int ip3 = ip0+1;

            glNormal3dv(dNorms[ip0]);
            glVertex3dv(dPoints[ip0]);
            glNormal3dv(dNorms[ip1]);
            glVertex3dv(dPoints[ip1]);
            glNormal3dv(dNorms[ip3]);
            glVertex3dv(dPoints[ip3]);
            glNormal3dv(dNorms[ip3]);
            glVertex3dv(dPoints[ip3]);
            glNormal3dv(dNorms[ip1]);
            glVertex3dv(dPoints[ip1]);
            glNormal3dv(dNorms[ip2]);
            glVertex3dv(dPoints[ip2]);
            /*
               glNormal3dv (dNorms[ip0]);
               glVertex3dv (dPoints[ip0]);
               glNormal3dv (dNorms[ip1]);
               glVertex3dv (dPoints[ip1]);
               glNormal3dv (dNorms[ip2]);
               glVertex3dv (dPoints[ip2]);
               glNormal3dv (dNorms[ip3]);
               glVertex3dv (dPoints[ip3]);
             */
        }
    }  // End for each parametric step
    glEnd();
}

void GLRenderer::DrawCappedCylinder(int nSeg, double dAxisA[3], double dAxisB[3], double dEndA[][3], double dEndB[][3],
                                    double dNormA[][3], double dNormB[][3], double dCapNormA[][3], double dCapNormB[][3])
{
    //Make sure we are in material mode
    //  if (m_iColorMode != SSI_MATERIAL)
    //    SetColorMode (SSI_MATERIAL);

    //First draw the cylinder
    //Construct the quadrilaterals
    glBegin(GL_QUADS);
    int i;
    for (i = 0; i < nSeg; i++)
    {
        int j = i+1;
        if (j >= nSeg)
        {
            j = 0;
        }
        glNormal3dv(dNormA[i]);
        glVertex3dv(dEndA[i]);
        glNormal3dv(dNormB[i]);
        glVertex3dv(dEndB[i]);
        glNormal3dv(dNormB[j]);
        glVertex3dv(dEndB[j]);
        glNormal3dv(dNormA[j]);
        glVertex3dv(dEndA[j]);
    }
    glEnd();

    //Now do the caps
    glBegin(GL_TRIANGLES);

    //Top A
    for (i = 0; i < nSeg; i++)
    {
        int j = i+1;
        if (j >= nSeg)
        {
            j = 0;
        }

        glNormal3dv(dCapNormA[i]);
        glVertex3dv(dAxisA);
        glNormal3dv(dCapNormA[j]);
        glVertex3dv(dEndA[j]);
        glNormal3dv(dCapNormA[i]);
        glVertex3dv(dEndA[i]);
    }

    //Bottom B (note reverse order)
    for (i = 0; i < nSeg; i++)
    {
        int j = i+1;
        if (j >= nSeg)
        {
            j = 0;
        }

        glNormal3dv(dCapNormB[i]);
        glVertex3dv(dAxisB);
        glNormal3dv(dCapNormB[i]);
        glVertex3dv(dEndB[i]);
        glNormal3dv(dCapNormB[j]);
        glVertex3dv(dEndB[j]);
    }
    glEnd();
}

void GLRenderer::DrawMaps(std::vector<EMap*> &maps)
{
    if (maps.size() <= 0)
    {
        return;
    }

    for (unsigned int i = 0; i < maps.size(); i++)
    {
        EMap *emap = (EMap*) maps[i];
        if (emap->Visible())
        {
            drawLines(emap->edges, (int) emap->GetMapLinewidth(), true);

        }
    }
}

void GLRenderer::DrawLabels(std::list<Molecule*> &molecules)
{

    // Collect and sort labels for proper transparency
    float view[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, view);
    Matrix4<float> viewmat(view);
    viewmat.transpose();

    typedef std::map<float, ATOMLABEL*> LabelOrderMap;
    LabelOrderMap order;
    std::list<Molecule*>::iterator node;
    for (node = molecules.begin(); node != molecules.end(); ++node)
    {
        if ((*node)->LabelsVisible())
        {
            Molecule::AtomLabelList molLabels = (*node)->getAtomLabels();
            Molecule::AtomLabelList::iterator mlIter;
            for (mlIter = molLabels.begin(); mlIter != molLabels.end(); ++mlIter)
            {
                ATOMLABEL *label = *mlIter;
                if (label->isVisible())
                {
                    const MIAtom &atom = *label->atom();
                    Tuple4<float> vec(atom.x(), atom.y(), atom.z(), 1.0f);
                    viewmat.transform(vec);
                    order.insert(std::make_pair(vec.z, label));
                }
            }
        }
    }
    Molecule::AtomLabelList labels;
    labels.reserve(order.size());
    foreach (const LabelOrderMap::value_type& l, order)
    {
        labels.push_back(l.second);
    }
    drawLabels(labels);
}

void GLRenderer::drawText(const char *text, float x, float y, float z)
{
    float inverseRotation[16];
    getInverseRotation(inverseRotation);
    drawText(text, x, y, z, 0.0f, labelsTextScale, inverseRotation, false, 0, white);
}

std::auto_ptr<mi::opengl::GLFont> GLRenderer::createGLFont()
{
    return std::auto_ptr<mi::opengl::GLFont>(new GLFont(qcanvas_win));
}

void GLRenderer::drawText(const char *text, float x, float y, float z, float offset, float scale, float inverseRotation[16], bool displayNumber, int number, float foreground[4])
{

    if (!labelsFont_.get())
    {
        return;
    }
    glPushMatrix();
    // Position at coordintes.
    glTranslatef(x, y, z);

    // Keep oriented toward camera.
    glMultMatrixf(inverseRotation);

    // Offset position.
    glTranslatef(offset, offset, offset);
    // Scale size.
    glScalef(scale, scale, 1.0f);

    glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_ENABLE_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    if (displayNumber)
    {
        glPushAttrib(GL_CURRENT_BIT);
        static char buffer[10];
        sprintf(buffer, " %d ", number);
        float numberWidth = labelsFont_->fontMetrics().width(QString(buffer));
        glColor3f(1.0f, 1.0f, 0.0f);
        float x1 = 0.0f;
        float y1 = -labelsFont_->fontMetrics().descent();
        float x2 = numberWidth;
        float y2 = labelsFont_->fontMetrics().ascent();
        float z = -0.02f;
        glBegin(GL_POLYGON);
        glVertex3f(x1, y1, z);
        glVertex3f(x2, y1, z);
        glVertex3f(x2, y2, z);
        glVertex3f(x1, y2, z);
        glEnd();
        float xOffset = (x2 - x1) / 10.0f;
        float yOffset = (y2 - y1) / 10.0f;
        x1 += xOffset;
        y1 += yOffset;
        x2 -= xOffset;
        y2 -= yOffset;
        z += 0.01f;
        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_POLYGON);
        glVertex3f(x1, y1, z);
        glVertex3f(x2, y1, z);
        glVertex3f(x2, y2, z);
        glVertex3f(x1, y2, z);
        glEnd();

        glEnable(GL_TEXTURE_2D);
        glColor3f(1.0f, 1.0f, 1.0f);
        glPushMatrix();
        labelsFont_->render(QString(buffer));
        glPopMatrix();

        glTranslatef(numberWidth, 0.0f, 0.0f);
        glPopAttrib();
    }

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_TEXTURE_2D);

    glColor4fv(foreground);
    glPushMatrix();
    labelsFont_->render(QString(text));
    glPopMatrix();

    glPopAttrib();
    glPopMatrix();
}

void GLRenderer::getInverseRotation(float inverseRotation[16])
{
    static GLfloat modelview[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, modelview);
    Matrix4<float> matrix(modelview);
    Quaternion<float> quat;
    quat.set(matrix);
    quat.negate();
    matrix.set(quat);
    matrix.getInColumnMajorOrder(inverseRotation);
}

void GLRenderer::drawLabels(Molecule::AtomLabelList &labels)
{
    if (labels.size() <= 0)
    {
        return;
    }

    float foreground[4];
    float inverseRotation[16];
    getInverseRotation(inverseRotation);

    // Diable lighting to get proper flat shading.
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    Molecule::AtomLabelList::iterator iter;
    for (iter = labels.begin(); iter != labels.end(); ++iter)
    {
        ATOMLABEL &label = **iter;
        const MIAtom *atom = label.atom();
        if (hideHydrogens && MIAtom::MIIsHydrogen(atom))
        {
            continue;
        }

        const std::string &text = label.label();
        float x = atom->x();
        float y = atom->y();
        float z = atom->z();
        float offset = atom->getRadius() * style.getBallPercent();
        foreground[0] = label.red() / 255.0f;
        foreground[1] = label.green() / 255.0f;
        foreground[2] = label.blue() / 255.0f;
        foreground[3] = 1.0f;

        drawText(text.c_str(), x, y, z, offset, labelsTextScale, inverseRotation, false, 0, foreground);
//    }
    }
    glPopAttrib();
}

void GLRenderer::DrawAnnotations(std::list<Molecule*> &molecules)
{
    std::list<Molecule*>::iterator node = molecules.begin();
    while (node != molecules.end())
    {
        drawAnnotations((*node)->getAnnotations(), (*node)->isDrawAnnotationBox());
        node++;
    }
}

void GLRenderer::drawAnnotations(Molecule::AnnotationList &Annotations, bool drawBox)
{
    if (Annotations.size() <= 0)
    {
        return;
    }

    float inverseRotation[16];
    getInverseRotation(inverseRotation);

    // Diable lighting to get proper flat shading.
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    int annotationNumber = 0;
    Molecule::AnnotationList::iterator iter = Annotations.begin();
    while (iter != Annotations.end())
    {
        Annotation *annotation = *iter;
        ++iter;
        ++annotationNumber;
        if (annotation->isHidden())
        {
            continue;
        }
        glPushName(getPickName(*annotation));
        const char *text = annotation->m_text.c_str();
        float x = annotation->m_x;
        float y = annotation->m_y;
        float z = annotation->m_z;
        bool displayNumber = drawBox;
        if (annotation->m_type == Annotation::Plane)
        {
            displayNumber = true;
        }
        float *foreground = getColor(annotation->m_color);
        drawText(text, x, y, z, 0.0f, labelsTextScale, inverseRotation, displayNumber, annotationNumber, foreground);
        glPopName();
    }
    glPopAttrib();
}

void GLRenderer::circleStackAtoms(Stack *stack)
{

    if (stack->empty())
    {
        return;
    }

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LINE_BIT);
    glPushMatrix();

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);

    float inverseRotation[16];
    getInverseRotation(inverseRotation);

    glLineWidth(2.0f);
    Stack::DataContainer::const_iterator iter = stack->getData().begin();
    while (iter != stack->getData().end())
    {
        StackItem item = *iter;
        ++iter;
        MIAtom *atom = item.atom;
        if (MIAtom::isValid(atom))
        {
            if (hideHydrogens && MIAtom::MIIsHydrogen(atom))
            {
                continue;
            }
            if (iter == stack->getData().end())
            {
                glColor3f(1.0f, 1.0f, 1.0f);
            }
            else
            {
                glColor3f(0.7f, 0.7f, 1.0f);
            }
            if (item.molecule->Visible())
            {

                glPushMatrix();
                // Position at coordintes.
                glTranslatef(atom->x(), atom->y(), atom->z());
                // Keep oriented toward camera.
                glMultMatrixf(inverseRotation);

                float radius = 1.2f * atom->getRadius() * style.getBallPercent();
                glBegin(GL_LINE_STRIP);
                float circleX = 0.0f;
                float circleY = radius;
                for (float angle = 0.0f; angle <= (2.0f*M_PI); angle += 0.01f)
                {
                    float x = radius * (float)sin((double)angle);
                    float y = radius * (float)cos((double)angle);
                    glVertex2d(circleX, circleY);
                    circleX = x;
                    circleY = y;
                }
                glEnd();
                glPopMatrix();

            }
        }
    }
    glPopMatrix();
    glPopAttrib();
    logOpenGLErrors(__FILE__, __LINE__);
}

void GLRenderer::drawStackText(int x, int y, const std::string &s)
{
    if (!stackFont_.get())
    {
        return;
    }
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_DEPTH_TEST);
    glPushMatrix();
    glTranslatef(x, y, 0.0f);
    stackFont_->render(QString(s.c_str()));
    glPopMatrix();
    glPopAttrib();
}

void GLRenderer::DrawStack(Stack *stack, int x, int y)
{

    if (!stackFont_.get())
    {
        return;
    }
    if (stack->empty())
    {
        return;
    }

    glPushAttrib(GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT);
    glPushMatrix();
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);

    glColor3f(1.0f, 1.0f, 1.0f);

    Stack::DataContainer::const_iterator iter = stack->getData().begin();
    std::string s;
    const int lineHeight = (int)(stackFont_->fontMetrics().ascent() - stackFont_->fontMetrics().descent() + 5);

    y += 5;
    if (!stack->isMinimized())
    {
        // list the stack in the lower-left corner of the screen
        // stack pushes up from the bottom
        glColor3f(0.7f, 0.7f, 1.0f);
        StackItem item = *iter;
        if (stack->size() > 4)
        {
            for (int n = 0; n < stack->size()-4; ++n)
            {
                ++iter;
            }
            // if more than 4 print a more...
            s = " + ";
            s += ftoa(stack->size()-4);
            s += " more...";
            drawStackText(x, y, s);
            y += lineHeight;
        }
        if (stack->size() > 3)
        {
            item = *iter;
            ++iter;
            if (MIAtom::isValid(item.atom) && Monomer::isValid(item.residue) && MIMoleculeBase::isValid(item.molecule) && item.molecule->Visible())
            {
                s = "4: ";
                s += item.residue->type().c_str();
                s += " ";
                s += item.residue->name().c_str();
                s += " ";
                s += item.atom->name();
                drawStackText(x, y, s);
                y += lineHeight;
            }
        }

        if (stack->size() > 2)
        {
            item = *iter;
            ++iter;
            if (MIAtom::isValid(item.atom) && Monomer::isValid(item.residue) && MIMoleculeBase::isValid(item.molecule) && item.molecule->Visible())
            {
                s = "3: ";
                s += item.residue->type().c_str();
                s += " ";
                s += item.residue->name().c_str();
                s += " ";
                s += item.atom->name();
                drawStackText(x, y, s);
                y += lineHeight;
            }
        }

        if (stack->size() > 1)
        {
            item = *iter;
            ++iter;
            if (MIAtom::isValid(item.atom) && Monomer::isValid(item.residue) && MIMoleculeBase::isValid(item.molecule) && item.molecule->Visible())
            {
                s = "2: ";
                s += item.residue->type().c_str();
                s += " ";
                s += item.residue->name().c_str();
                s += " ";
                s += item.atom->name();
                drawStackText(x, y, s);
                y += lineHeight;
            }
        }

        glColor3f(1.0f, 1.0f, 1.0f);
        item = *iter;
        ++iter;
        if (MIAtom::isValid(item.atom) && Monomer::isValid(item.residue) && MIMoleculeBase::isValid(item.molecule) && item.molecule->Visible())
        {
            s = "1: ";
            s += item.residue->type().c_str();
            s += " ";
            s += item.residue->name().c_str();
            s += " ";
            s += item.atom->name();
            drawStackText(x, y, s);
            y += lineHeight;
        }
    }
    else
    {
        glColor3f(1.0f, 1.0f, 1.0f);
        std::string str("...");
        drawStackText(x, y, str);
    }

    CRect &PopBox = stack->getPopBox();
    CRect &HideBox = stack->getHideBox();
    CRect &ClearBox = stack->getClearBox();

    x += 50;
    const int spacing = 5;
    PopBox.left = x;
    PopBox.right = x + popStackImage->getWidth();
    PopBox.bottom = y + popStackImage->getHeight();
    PopBox.top = y;
    glRasterPos2i(x, y);
    glDrawPixels(popStackImage->getWidth(), popStackImage->getHeight(), popStackImage->getFormat(), GL_UNSIGNED_BYTE, popStackImage->getData());
    x += popStackImage->getWidth();

    x += spacing;
    HideBox.left = x;
    HideBox.right = x + minimizeStackImage->getWidth();
    HideBox.bottom = y + minimizeStackImage->getHeight();
    HideBox.top = y;
    glRasterPos2i(x, y);
    if (stack->isMinimized())
    {
        glDrawPixels(unminimizeStackImage->getWidth(), unminimizeStackImage->getHeight(), unminimizeStackImage->getFormat(), GL_UNSIGNED_BYTE, unminimizeStackImage->getData());
        x += unminimizeStackImage->getWidth();
    }
    else
    {
        glDrawPixels(minimizeStackImage->getWidth(), minimizeStackImage->getHeight(), minimizeStackImage->getFormat(), GL_UNSIGNED_BYTE, minimizeStackImage->getData());
        x += minimizeStackImage->getWidth();
    }

    x += spacing;
    ClearBox.left = x;
    ClearBox.right = x + closeStackImage->getWidth();
    ClearBox.bottom = y + closeStackImage->getHeight();
    ClearBox.top = y;
    glRasterPos2i(x, y);
    glDrawPixels(closeStackImage->getWidth(), closeStackImage->getHeight(), closeStackImage->getFormat(), GL_UNSIGNED_BYTE, closeStackImage->getData());

    glPopMatrix();
    glPopAttrib();
}

void GLRenderer::CircleAtoms(MIAtom **atoms, int natoms)
{

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_LINE_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);

    float inverseRotation[16];
    getInverseRotation(inverseRotation);

    glLineWidth(2.0f);
    logOpenGLErrors(__FILE__, __LINE__);
    glColor3f(1.0f, 0.0f, 0.0f);
    for (int index = 0; index < natoms; index++)
    {
        MIAtom *atom = atoms[index];
        if (hideHydrogens && MIAtom::MIIsHydrogen(atom))
        {
            continue;
        }

        if (atom != NULL)
        {

            glPushMatrix();
            // Position at coordintes.
            glTranslatef(atom->x(), atom->y(), atom->z());
            // Keep oriented toward camera.
            glMultMatrixf(inverseRotation);

            float radius = 1.2f * atom->getRadius() * style.getBallPercent();
            glBegin(GL_LINE_STRIP);
            float circleX = 0.0f;
            float circleY = radius;
            for (float angle = 0.0f; angle <= (2.0f*M_PI); angle += 0.01f)
            {
                float x = radius * (float)sin((double)angle);
                float y = radius * (float)cos((double)angle);
                glVertex2d(circleX, circleY);
                circleX = x;
                circleY = y;
            }
            glEnd();
            glPopMatrix();

        }
    }
    glPopAttrib();
    logOpenGLErrors(__FILE__, __LINE__);

}

void GLRenderer::DrawMessage(const std::string &message, int x, int y)
{
    drawStackText(x, y, message);
}

void GLRenderer::DrawContacts(std::vector<CONTACT> &Contacts)
{
    if (Contacts.size() <= 0)
    {
        return;
    }

    float inverseRotation[16];
    getInverseRotation(inverseRotation);
    char text[20];

    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_LINE_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glDisable(GL_LIGHTING);

    glLineWidth(2.0f);
    for (unsigned int i = 0; i < Contacts.size(); i++)
    {
        MIAtom &a1 = *Contacts[i].line.getAtom1();
        MIAtom &a2 = *Contacts[i].line.getAtom2();
        if (Contacts[i].color > 0)
        {
            sprintf(text, "%0.2f", Contacts[i].d);
            glColor4fv(getColor(Contacts[i].color));
            glBegin(GL_LINES);
            glVertex3f(a1.x(), a1.y(), a1.z());
            glVertex3f(a2.x(), a2.y(), a2.z());
            glEnd();
            float x = (a1.x() + a2.x()) * 0.5f;
            float y = (a1.y() + a2.y()) * 0.5f;
            float z = (a1.z() + a2.z()) * 0.5f;
            drawText(text, x, y, z, 0.0f, labelsTextScale, inverseRotation, false, 0, white);
        }
    }
    glPopAttrib();
}

void GLRenderer::DrawVus(std::vector<PLINE> &Vus, int w)
{
    drawLines(Vus, w, false);
}

void GLRenderer::drawLines(std::vector<PLINE> &Vus, int w, bool withDepthTest)
{
    if (Vus.size() <= 0)
    {
        return;
    }

    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_LINE_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    if (!withDepthTest)
    {
        glDisable(GL_DEPTH_TEST);
    }
    glDisable(GL_LIGHTING);

    glLineWidth(w);
    for (unsigned int i = 0; i < Vus.size(); i++)
    {
        APOINT &a1 = Vus[i].p1;
        APOINT &a2 = Vus[i].p2;
        glColor4fv(getColor(a1.color));
        glBegin(GL_LINES);
        glVertex3f(a1.x, a1.y, a1.z);
        glVertex3f(a2.x, a2.y, a2.z);
        glEnd();
    }
    glPopAttrib();
}

float*GLRenderer::getColor(int colorIndex, bool adjustForCurrentModel)
{
    static float color[4];
    float scaling = 1.0f;
    if (adjustForCurrentModel && !currentModel)
    {
        scaling = 1.0 - amountToDimNonactiveModels;
    }
    int c = PaletteIndex(colorIndex);
    color[0] = Colors::RPallette[c] / 255.0f * scaling;
    color[1] = Colors::GPallette[c] / 255.0f * scaling;
    color[2] = Colors::BPallette[c] / 255.0f * scaling;
    color[3] = 1.0f;
    return color;
}

float*GLRenderer::getColor(PaletteColor c, bool adjustForCurrentModel)
{
    static float color[4];
    float scaling = 1.0f;
    if (adjustForCurrentModel && !currentModel)
    {
        scaling = 1.0 - amountToDimNonactiveModels;
    }
    color[0] = c.red / 255.0f * scaling;
    color[1] = c.green / 255.0f * scaling;
    color[2] = c.blue / 255.0f * scaling;
    color[3] = 1.0f;
    return color;
}

void GLRenderer::drawAngles(std::vector<ANGLE> &angles)
{
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(4.0f);
    for (unsigned int i = 0; i < angles.size(); i++)
    {
        ANGLE &angle = angles[i];
        glPushName(getPickName(angle));
        MIAtom *a1 = angles[i].getAtom1();
        MIAtom *a2 = angles[i].atom3;
        if (hideHydrogens && (MIAtom::MIIsHydrogen(a1) || MIAtom::MIIsHydrogen(a2)))
        {
            continue;
        }
        glBegin(GL_LINES);
        glVertex3f(a1->x(), a1->y(), a1->z());
        glVertex3f(a2->x(), a2->y(), a2->z());
        glEnd();
        glPopName();
    }
    glPopAttrib();
}

void GLRenderer::drawBonds(std::vector<Bond> &bonds, int *color)
{
    Vector3<float> eye;
    if (camera != NULL)
    {
        eye.set(camera->getEye());
    }
    int lineDepthCueSteps = 4;
    float lineDepthCueIncrement = (frustum->getFarClipping() - frustum->getNearClipping()) / lineDepthCueSteps;
    float lineDepthCueStart = frustum->getNearClipping();
    std::vector<Bond>::iterator bondsIter = bonds.begin();
    std::vector<Bond>::iterator bondsEnd = bonds.end();
    Vector3<float> pos1;
    Vector3<float> pos2;
    Vector3<float> diff;
    while (bondsIter != bondsEnd)
    {
        Bond &bond = *bondsIter;
        ++bondsIter;
        MIAtom *a1 = bond.getAtom1();
        MIAtom *a2 = bond.getAtom2();
        if (hideHydrogens && (MIAtom::MIIsHydrogen(a1) || MIAtom::MIIsHydrogen(a2)))
        {
            continue;
        }

        pos1.set(a1->x(), a1->y(), a1->z());
        pos2.set(a2->x(), a2->y(), a2->z());
        diff = pos1 - pos2;
        float separation = diff.length();
        if (a1->color() > 0 && a2->color() > 0)
        {
            if (frustum->sphereInFrustum(pos1, separation) != Frustum::OUTSIDE
                && frustum->sphereInFrustum(pos2, separation) != Frustum::OUTSIDE)
            {

                int c1 = ::getColor(a1);
                int c2 = ::getColor(a2);
                if (color != NULL)
                {
                    c1 = *color;
                    c2 = *color;
                }
                unsigned char bondOrder = bond.getOrder();
                if (!showBondOrders)
                {
                    bondOrder = NORMALBOND;
                }
                glPushName(getPickName(bond));
                if (bond.type == B_POINT || bond.type == B_SYMM_POINT)
                {
                    if (!style.isAtomBall())
                    {
                        float radius = a1->getRadius() * style.getBallPercent();
                        if (radius > 0.0f && a1->color() > 0)
                        {
                            int c1 = ::getColor(a1);
                            glColor4fv(getColor(c1));
                            float p1[3];
                            p1[0] = a1->x();
                            p1[1] = a1->y();
                            p1[2] = a1->z();
                            DrawSphere(p1, radius, 16);
                        }
                    }
                }
                else if (style.isBondCylinder())
                {
                    float radius = std::min(a1->getRadius(), a2->getRadius());
                    radius *= style.getStickPercent();

                    if (viewVectorSet && (bondOrder == DOUBLEBOND || bondOrder == PARTIALDOUBLEBOND || bondOrder == TRIPLEBOND))
                    {
                        radius /= 2.0f;
                    }
                    if (!viewVectorSet || bondOrder == NORMALBOND || bondOrder == SINGLEBOND || bondOrder == TRIPLEBOND
                        || bondOrder == IONICBOND || bondOrder == METALLIGANDBOND)
                    {
                        drawBondCylinder(pos1, c1, pos2, c2, radius, true);
                    }
                    if (viewVectorSet && (bondOrder == DOUBLEBOND || bondOrder == PARTIALDOUBLEBOND || bondOrder == TRIPLEBOND))
                    {
                        Vector3<float> offset;
                        Vector3<float> bondVector(pos2 - pos1);
                        offset.cross(viewVector, bondVector);
                        float scaleFactor = 2.0f*radius;
                        if (bondOrder == TRIPLEBOND)
                        {
                            scaleFactor *= 1.5;
                        }
                        offset.scale(scaleFactor / offset.length());
                        Vector3<float> p1 = pos1 + offset;
                        Vector3<float> p2 = pos2 + offset;
                        drawBondCylinder(p1, c1, p2, c2, radius, true);
                        p1 = pos1 - offset;
                        p2 = pos2 - offset;
                        if (bondOrder == PARTIALDOUBLEBOND)
                        {
                            bondVector.scale(0.3f);
                            p1.add(bondVector);
                            p2.subtract(bondVector);
                        }
                        drawBondCylinder(p1, c1, p2, c2, radius, true);
                    }
                }
                else if (style.isBondLine())
                {
                    float lineWidth = style.getBondLineWidth();
                    if (lineWidthDepthCued)
                    {
                        float bondDist = std::min(Vector3<float>(pos1 - eye).length(), Vector3<float>(pos2 - eye).length());
                        for (int i = 1; i <= lineDepthCueSteps; ++i)
                        {
                            if (bondDist < (lineDepthCueStart + lineDepthCueIncrement * i))
                            {
                                lineWidth += (float) (lineDepthCueSteps - i);
                                break;
                            }
                        }
                    }
                    if (bondOrder == NORMALBOND || bondOrder == SINGLEBOND || bondOrder == TRIPLEBOND
                        || bondOrder == IONICBOND || bondOrder == METALLIGANDBOND)
                    {
                        drawBondLine(pos1, c1, pos2, c2, lineWidth);
                    }
                    if (bondOrder == DOUBLEBOND || bondOrder == PARTIALDOUBLEBOND || bondOrder == TRIPLEBOND)
                    {
                        Vector3<float> offset;
                        Vector3<float> bondVector(pos2 - pos1);
                        offset.cross(viewVector, bondVector);
                        offset.scale(0.08 / offset.length());
                        Vector3<float> p1 = pos1 + offset;
                        Vector3<float> p2 = pos2 + offset;
                        drawBondLine(p1, c1, p2, c2, lineWidth);
                        p1 = pos1 - offset;
                        p2 = pos2 - offset;
                        if (bondOrder == PARTIALDOUBLEBOND)
                        {
                            bondVector.scale(0.3f);
                            p1.add(bondVector);
                            p2.subtract(bondVector);
                        }
                        drawBondLine(p1, c1, p2, c2, lineWidth);
                    }
                }
                glPopName();
            }
        }
    }
}

void GLRenderer::drawBondCylinder(const Vector3<float> &pos1, int color1, const Vector3<float> &pos2, int color2, float radius, bool capped)
{
    if (radius <= 0.0f)
    {
        return;
    }
    Vector3<float> midPoint = pos1 + pos2;
    midPoint.scale(0.5);

    glPushName(1);
    glColor4fv(getColor(color1));
    DrawCylinder(pos1, midPoint, radius);
    float p[3];
    if (capped)
    {
        p[0] = pos1.x;
        p[1] = pos1.y;
        p[2] = pos1.z;
        DrawSphere(p, radius, 16);
    }
    glPopName();

    glPushName(2);
    glColor4fv(getColor(color2));
    DrawCylinder(midPoint, pos2, radius);
    if (capped)
    {
        p[0] = pos2.x;
        p[1] = pos2.y;
        p[2] = pos2.z;
        DrawSphere(p, radius, 16);
    }
    glPopName();
}

void GLRenderer::drawBondLine(const Vector3<float> &pos1, int color1, const Vector3<float> &pos2, int color2, float lineWidth)
{
    if (lineWidth <= 0.0f)
    {
        return;
    }
    glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_LIGHTING_BIT);
    if (antialiasLines)
    {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    }
    else
    {
        //    glDisable(GL_BLEND);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDisable(GL_LINE_SMOOTH);
    }
    glDisable(GL_LIGHTING);
    glLineWidth(lineWidth);

    if (joinBondsOfSameColor && color1 == color2)
    {
        glBegin(GL_LINES);
        glColor4fv(getColor(color1, true));
        glVertex3f(pos1.x, pos1.y, pos1.z);
        glVertex3f(pos2.x, pos2.y, pos2.z);
        glEnd();
    }
    else
    {
        Vector3<float> midPoint = pos1 + pos2;
        midPoint.scale(0.5f);
        glPushName(1);
        glBegin(GL_LINES);
        glColor4fv(getColor(color1, true));
        glVertex3f(pos1.x, pos1.y, pos1.z);
        glVertex3f(midPoint.x, midPoint.y, midPoint.z);
        glEnd();
        glPopName();
        glPushName(2);
        glBegin(GL_LINES);
        glColor4fv(getColor(color2, true));
        glVertex3f(midPoint.x, midPoint.y, midPoint.z);
        glVertex3f(pos2.x, pos2.y, pos2.z);
        glEnd();
        glPopName();
    }
    glPopAttrib();
}

void GLRenderer::drawResidueAtoms(MIIter<Residue> &res)
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Vector3<float> pos;
    for (; res; ++res)
    {
        for (unsigned int i = 0; i < (unsigned int) res->atomCount(); i++)
        {
            MIAtom *a1 = res->atom(i);
            if (hideHydrogens && MIAtom::MIIsHydrogen(a1))
            {
                continue;
            }
            float radius = a1->getRadius() * style.getBallPercent();
            if (radius > 0.0f && a1->color() > 0)
            {
                pos.set(a1->x(), a1->y(), a1->z());
                if (frustum->sphereInFrustum(pos, radius) != Frustum::OUTSIDE)
                {
                    glPushName(getPickName(a1));
                    int c1 = ::getColor(a1);
                    glColor4fv(getColor(c1));
                    float p1[3];
                    p1[0] = a1->x();
                    p1[1] = a1->y();
                    p1[2] = a1->z();
                    DrawSphere(p1, radius, 16);
                    glPopName();
                }
            }
        }
    }
}

void GLRenderer::drawMolecule(Molecule *molecule)
{
    if (molecule->Visible())
    {
        drawBonds(molecule->getBonds());
        int hBondColor = Colors::WHITE;
        if (molecule->hbonds.size() != molecule->hbondContacts.size())
        {
            std::vector<CONTACT>().swap(molecule->hbondContacts); // was molecule->hbondContacts.clear();
            std::vector<Bond>::iterator hbond = molecule->hbonds.begin();
            while (hbond != molecule->hbonds.end())
            {
                Bond &bond = *hbond;
                ++hbond;
                CONTACT contact;
                contact.line = bond;
                contact.color = hBondColor;
                contact.d = bond.getAtom1()->distance(bond.getAtom2());
                molecule->hbondContacts.push_back(contact);
            }
        }
        DrawContacts(molecule->hbondContacts);
        if (style.isAtomBall())
        {
            MIIter<Residue> res = molecule->GetResidues();
            drawResidueAtoms(res);
        }
        if (molecule->DotsVisible())
        {
            drawSurface(molecule->getDots());
        }
    }
}

void GLRenderer::drawSymmetryMolecule(Molecule *molecule, bool showSymmetryAsBackbone)
{
    if (molecule->Visible())
    {
        drawBonds(molecule->getSymmetryBonds());
        if (!showSymmetryAsBackbone && style.isAtomBall())
        {
            MIIter<Residue> res = molecule->GetSymmResidues();
            drawResidueAtoms(res);
        }
    }
}

GLuint GLRenderer::getPickName(MIAtom *atom)
{
    GLuint name = atomPickNameMap.size()+1;
    if (pickingEnabled)
    {
        atomPickNameMap.insert(std::map<GLuint, MIAtom*>::value_type(name, atom));
    }
    return name;
}

GLuint GLRenderer::getPickName(Annotation &annotation)
{
    GLuint name = annotationPickNameMap.size()+1;
    if (pickingEnabled)
    {
        annotationPickNameMap.insert(std::map<GLuint, Annotation*>::value_type(name, &annotation));
    }
    return name;
}

GLuint GLRenderer::getPickName(Bond &bond)
{
    GLuint name = bondPickNameMap.size()+1;
    if (pickingEnabled)
    {
        bondPickNameMap.insert(std::map<GLuint, Bond*>::value_type(name, &bond));
    }
    return name;
}

GLuint GLRenderer::getPickName(ANGLE &angle)
{
    GLuint name = anglePickNameMap.size()+1;
    if (pickingEnabled)
    {
        anglePickNameMap.insert(std::map<GLuint, ANGLE*>::value_type(name, &angle));
    }
    return name;
}

void GLRenderer::DrawCylinder(const Vector3<float> &start, const Vector3<float> &end, float radius)
{
#define   c45  0.70711f

    float vp[16][6];
    float mat[3][3];
    float a[3];
    float b[3];
    float c[3];
    float d[3];
    int i, j;
    static float cdata[8][3] =  {{0.0f, 1.0f, 0.0f},
                                 {c45, c45, 0.0f},
                                 {1.0f, 0.0f, 0.0f},
                                 {c45, -c45, 0.0f},
                                 {0.0f, -1.0f, 0.0f},
                                 {-c45, -c45, 0.0f},
                                 {-1.0f, 0.0f, 0.0f},
                                 {-c45, c45, 0.0f}};

    //Construct a transformation matrix to map the cylinder to the endpoints given
    for (i = 0; i < 3; i++)
    {
        a[i] = end[i] - start[i];
        b[i] = 0.0;
    }
    normalize(a); //Axis direction

    if (fabs(a[0]) < 0.9)
    {
        b[0] = 1.0;
    }
    else
    {
        b[1] = 1.0; // random unit vector
    }

    cross(a, b, c);  // perpendicular to a
    normalize(c);
    cross(a, c, b);  // perpendicular to a, c
    normalize(b);
    for (i = 0; i < 3; i++)
    {
        mat[0][i] = c[i];
        mat[1][i] = b[i];
        mat[2][i] = a[i];
    }

    //Construct the transformed coordinates for the cyclinder
    for (i = 0; i < 8; i++)
    {
        vecmat(cdata[i], mat, d);
        for (j = 0; j < 3; j++)
        {
            vp[i][j] = d[j]*radius + start[j];
            vp[i][j+3] = d[j];
            vp[i+8][j] = d[j]*radius + end[j];
            vp[i+8][j+3] = d[j];
        }
    }

    //Construct the quadrilaterals
    glBegin(GL_QUADS);
    for (i = 0; i < 8; i++)
    {
        j = i+1;
        if (j >= 8)
        {
            j = 0;
        }
        if (useClockwiseFrontFaces)
        {
            glNormal3fv(&vp[j][3]);
            glVertex3fv(&vp[j][0]);
            glNormal3fv(&vp[j+8][3]);
            glVertex3fv(&vp[j+8][0]);
            glNormal3fv(&vp[i+8][3]);
            glVertex3fv(&vp[i+8][0]);
            glNormal3fv(&vp[i][3]);
            glVertex3fv(&vp[i][0]);
        }
        else
        {
            glNormal3fv(&vp[i][3]);
            glVertex3fv(&vp[i][0]);
            glNormal3fv(&vp[i+8][3]);
            glVertex3fv(&vp[i+8][0]);
            glNormal3fv(&vp[j+8][3]);
            glVertex3fv(&vp[j+8][0]);
            glNormal3fv(&vp[j][3]);
            glVertex3fv(&vp[j][0]);
        }
    }
    glEnd();
}

void GLRenderer::normalize(float a[3])
{
    // Find length and divide each term
    float len = sqrt((a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
    if (len > .000001)
    {
        for (int i = 0; i < 3; i++)
        {
            a[i] = a[i]/len;
        }
    }

}

void GLRenderer::cross(float a[3], float b[3], float c[3])
{

    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];

}

void GLRenderer::vecmat(float a[3], float mat[3][3], float b[3])
{
    for (int i = 0; i < 3; i++)
    {
        b[i] = a[0]*mat[0][i] + a[1]*mat[1][i] + a[2]*mat[2][i];
    }
}

void GLRenderer::setIdentity(float rotation[])
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

void GLRenderer::DrawSphere(float c[3], float radius, int tess)
{

    if (!quad[getContext()])
        quad[getContext()] = gluNewQuadric();

    glPushMatrix();
    glTranslated(c[0], c[1], c[2]);
    glScaled(radius, radius, radius);
    //  static float center[] = { 0.0, 0.0, 0.0 };
    //  DrawSphere2(0, center, 1.0, tess);
    gluQuadricNormals(quad[getContext()], GLU_SMOOTH);
    gluSphere(quad[getContext()], 1.0, tess, tess);

    glPopMatrix();
}

void GLRenderer::DrawSphere2(float c[3], float radius, int tess)
{
    //Translate and draw
    if (abs(tess) <= 4)
    {
        //if (!m_quad)
        //m_quad = gluNewQuadric ();
        //gluSphere (m_quad, radius, tess, tess);
        //    RenderCube (c, radius*0.8, tess); //Full cubes look too big
    }
    else if (tess < 0)
    {
        //    RenderTetra (c, radius, tess);
    }
    else
    {
        glPushMatrix();
        glTranslated(c[0], c[1], c[2]);
        float dradius[3];
        for (int i = 0; i < 3; i++)
        {
            dradius[i] = radius;
        }
        RenderSphereoid(c, dradius, tess);
        glPopMatrix();
    }
}

#define sphereoidX 0.525731112119133606f
#define sphereoidZ 0.850650808352039932f

void GLRenderer::RenderSphereoid(float*, float radius[3], int tess)
{
    static float fNormal[12][3] =
    {
        {-sphereoidX, 0.0f, sphereoidZ}, {sphereoidX, 0.0f, sphereoidZ}, {-sphereoidX, 0.0f, -sphereoidZ}, {sphereoidX, 0.0f, -sphereoidZ},
        {0.0f, sphereoidZ, sphereoidX}, {0.0f, sphereoidZ, -sphereoidX}, {0.0f, -sphereoidZ, sphereoidX}, {0.0f, -sphereoidZ, -sphereoidX},
        {sphereoidZ, sphereoidX, 0.0f}, {-sphereoidZ, sphereoidX, 0.0f}, {sphereoidZ, -sphereoidX, 0.0f}, {-sphereoidZ, -sphereoidX, 0.0f}
    };

    // Clockwise front faces
    //  static int tindeces[20][3] = {
    //    {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
    //    {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
    //    {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
    //    {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
    //  };
    // Counter-clockwise front faces
    static int tindeces[20][3] =
    {
        {1, 4, 0}, {4, 9, 0}, {4, 5, 9}, {8, 5, 4}, {1, 8, 4},
        {1, 10, 8}, {10, 3, 8}, {8, 3, 5}, {3, 2, 5}, {3, 7, 2},
        {3, 10, 7}, {10, 6, 7}, {6, 11, 7}, {6, 0, 11}, {6, 1, 0},
        {10, 1, 6}, {11, 0, 9}, {2, 11, 9}, {5, 2, 9}, {11, 2, 7}
    };
    int i;

    //if tess is negative then draw as a wire frame
    bool bWireFrame = false;
    if (tess < 0)
    {
        tess = -tess;
        bWireFrame = true;
    }

    //Set the normals and adjust for radius
    float fData[12][3];
    float fRadius[3];
    for (i = 0; i < 3; i++)
    {
        fRadius[i] = (float)radius[i];
    }
    for (i = 0; i < 12; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            fData[i][j] = fNormal[i][j]*fRadius[j];
        }
    }

    //Since this is a sphere, all the midpoint vectors will have the same length.
    float v[3];
    for (i = 0; i < 3; i++)
    {
        v[i] = fNormal[0][i] + fNormal[4][i];
    }
    float fnorm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    int level = 2;
    if (tess < 10)
    {
        level = 1;
    }
    if (bWireFrame)
    {
        for (i = 0; i < 20; i++)
        {
            TriangleRender(bWireFrame, fNormal[tindeces[i][0]], fNormal[tindeces[i][1]], fNormal[tindeces[i][2]],
                           fData[tindeces[i][0]], fData[tindeces[i][1]], fData[tindeces[i][2]]);
        }
    }
    else
    {
        for (i = 0; i < 20; i++)
        {
            Subdivide(level, bWireFrame, fnorm, fRadius,
                      fNormal[tindeces[i][0]], fNormal[tindeces[i][1]], fNormal[tindeces[i][2]],
                      fData[tindeces[i][0]], fData[tindeces[i][1]], fData[tindeces[i][2]]);
        }
    }
}

void GLRenderer::Subdivide(int level, bool bWireFrame, float fnorm, float radius[3],
                           float n1[3], float n2[3], float n3[3],
                           float v1[3], float v2[3], float v3[3])
{
    float v12[3], v23[3], v31[3], n12[3], n23[3], n31[3];
    for (int i = 0; i < 3; i++)
    {
        //This is effectively a normalize since these are all the same lenght.
        n12[i] = (n1[i]+n2[i])/fnorm;
        n23[i] = (n2[i]+n3[i])/fnorm;
        n31[i] = (n3[i]+n1[i])/fnorm;
    }
    /*
       float fnorm12 = sqrt ((double)(n12[0]*n12[0] + n12[1]*n12[1] + n12[2]*n12[2]));
       float fnorm23 = sqrt ((double)(n23[0]*n23[0] + n23[1]*n23[1] + n23[2]*n23[2]));
       float fnorm31 = sqrt ((double)(n31[0]*n31[0] + n31[1]*n31[1] + n31[2]*n31[2]));

       for (i=0; i<3; i++)
       {
        n12[i] = n12[i]/fnorm12;
        n23[i] = n23[i]/fnorm23;
        n31[i] = n31[i]/fnorm31;
       }
     */

    for (int i2 = 0; i2 < 3; i2++)
    {
        v12[i2] = n12[i2]*radius[i2];
        v23[i2] = n23[i2]*radius[i2];
        v31[i2] = n31[i2]*radius[i2];
    }

    if (level == 1)
    {
        TriangleRender(bWireFrame, n1, n12, n31, v1, v12, v31);
        TriangleRender(bWireFrame, n2, n23, n12, v2, v23, v12);
        TriangleRender(bWireFrame, n3, n31, n23, v3, v31, v23);
        TriangleRender(bWireFrame, n12, n23, n31, v12, v23, v31);
    }
    else
    {
        float v[3];
        for (int i = 0; i < 3; i++)
        {
            v[i] = n1[i] + n12[i];
        }
        fnorm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        level--;
        Subdivide(level, bWireFrame, fnorm, radius, n1, n12, n31, v1, v12, v31);
        Subdivide(level, bWireFrame, fnorm, radius, n2, n23, n12, v2, v23, v12);
        Subdivide(level, bWireFrame, fnorm, radius, n3, n31, n23, v3, v31, v23);
        Subdivide(level, bWireFrame, fnorm, radius, n12, n23, n31, v12, v23, v31);
    }
}

void GLRenderer::TriangleRender(bool bWireFrame, float n1[3], float n2[3], float n3[3],
                                float v1[3], float v2[3], float v3[3])
{
    if (bWireFrame)
    {
        glBegin(GL_LINE_LOOP);
        glVertex3fv(v1);
        glVertex3fv(v2);
        glVertex3fv(v3);
        glEnd();
    }
    else
    {
        glBegin(GL_TRIANGLES);
        glNormal3fv(n1);
        glVertex3fv(v1);
        glNormal3fv(n2);
        glVertex3fv(v2);
        glNormal3fv(n3);
        glVertex3fv(v3);
        glEnd();
    }

    /* Display normals
       float vn1[3], vn2[3], vn3[3];
       for (int i=0; i<3; i++)
       {
       vn1[i] = v1[i] + n1[i];
       vn2[i] = v2[i] + n2[i];
       vn3[i] = v3[i] + n3[i];
       }
       unsigned char bWhite[4] = {255,255,255, 0};
       SetColorMode (SSI_LINES);
       SetColor (bWhite);
       glBegin (GL_LINES);
       glVertex3fv (v1);
       glVertex3fv (vn1);
       glVertex3fv (v2);
       glVertex3fv (vn2);
       glVertex3fv (v3);
       glVertex3fv (vn3);
       glEnd();
     */
}

bool GLRenderer::isFogEnabled()
{
    return fogEnabled;
}

void GLRenderer::setFogEnabled(bool on)
{
    fogEnabled = on;
}

float GLRenderer::getFogStart()
{
    return fogStart;
}

void GLRenderer::setFogStart(float z)
{
    fogStart = z;
}

float GLRenderer::getFogEnd()
{
    return fogEnd;
}

void GLRenderer::setFogEnd(float z)
{
    fogEnd = z;
}

RenderStyle&GLRenderer::getRenderStyle()
{
    return style;
}

void GLRenderer::setRenderStyle(const RenderStyle &style)
{
    this->style.set(style);
}

void GLRenderer::setHideHydrogens(bool on)
{
    hideHydrogens = on;
}

bool GLRenderer::isHideHydrogens()
{
    return hideHydrogens;
}

void GLRenderer::setFrustum(Frustum *frustum)
{
    this->frustum = frustum;
}

void GLRenderer::setCamera(Camera *camera)
{
    this->camera = camera;
}

MIAtom*GLRenderer::getAtom(int id)
{
    MIAtom *atom = NULL;
    if (id > 0)
    {
        std::map<GLuint, MIAtom*>::iterator pos = atomPickNameMap.find(id);
        if (pos != atomPickNameMap.end())
        {
            atom = pos->second;
        }
    }
    return atom;
}

Bond*GLRenderer::getBond(int id)
{
    Bond *bond = NULL;
    if (id > 0)
    {
        std::map<GLuint, Bond*>::iterator pos = bondPickNameMap.find(id);
        if (pos != bondPickNameMap.end())
        {
            bond = pos->second;
        }
    }
    return bond;
}

ANGLE*GLRenderer::getAngle(int id)
{
    ANGLE *angle = NULL;
    if (id > 0)
    {
        std::map<GLuint, ANGLE*>::iterator pos = anglePickNameMap.find(id);
        if (pos != anglePickNameMap.end())
        {
            angle = pos->second;
        }
    }
    return angle;
}

Annotation*GLRenderer::getAnnotation(int id)
{
    Annotation *annotation = NULL;
    if (id > 0)
    {
        std::map<GLuint, Annotation*>::iterator pos = annotationPickNameMap.find(id);
        if (pos != annotationPickNameMap.end())
        {
            annotation = pos->second;
        }
    }
    return annotation;
}

void GLRenderer::clearPickNames()
{
    atomPickNameMap.clear();
    bondPickNameMap.clear();
    anglePickNameMap.clear();
    annotationPickNameMap.clear();
}

bool GLRenderer::isAntialiasLines()
{
    return antialiasLines;
}

void GLRenderer::setAntialiasLines(bool antialiasLines)
{
    this->antialiasLines = antialiasLines;
}

bool GLRenderer::isLineWidthDepthCued()
{
    return lineWidthDepthCued;
}

void GLRenderer::setLineWidthDepthCued(bool lineWidthDepthCued)
{
    this->lineWidthDepthCued = lineWidthDepthCued;
}

bool GLRenderer::isDimNonactiveModels()
{
    return dimNonactiveModels;
}

void GLRenderer::setDimNonactiveModels(bool dimNonactiveModels)
{
    this->dimNonactiveModels = dimNonactiveModels;
}

void GLRenderer::setAmountToDimNonactiveModels(float percent)
{
    amountToDimNonactiveModels = percent;
    percent = 1.0f - percent;
    light2->setAmbientColor(percent, percent, percent, 1.0f);
}

int GLRenderer::getFontSize()
{
    return fontSize;
}

void GLRenderer::setFontSize(int fontSize)
{
    this->fontSize = fontSize;
}

void GLRenderer::drawSurface(std::vector<SURFDOT> &dots)
{
    if (dots.size() <= 0)
    {
        return;
    }
    glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glPointSize(2.0f);
    glBegin(GL_POINTS);
    for (unsigned int i = 0; i < dots.size(); i++)
    {
        SURFDOT &dot = dots[i];
        if (dot.color > 0)
        {
            glColor4fv(getColor(dot.color));
            glVertex3f(dot.x, dot.y, dot.z);
        }
    }
    glEnd();
    glPopAttrib();
}

void GLRenderer::setJoinBondsOfSameColor(bool on)
{
    joinBondsOfSameColor = on;
}

void GLRenderer::setShowBondOrders(bool on)
{
    showBondOrders = on;
}

void GLRenderer::setPickingEnabled(bool value)
{
    pickingEnabled = value;
}

