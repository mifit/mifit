#ifndef mifit_ui_GLRenderer_h
#define mifit_ui_GLRenderer_h

#include "core/corelib.h"

#include <math/Vector3.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

class Displaylist;
class EMap;
class Helix;
class RibbonSegment;
class RibbonSpan;
class SecondaryStructure;
class TargaImage;
class QGLWidget;
class QVector4D;
namespace mi
{
    namespace opengl
    {
        class Light;
        class Frustum;
        class Camera;
        class Text;
    }
}

class GLRenderer
{
    std::map<void*, GLUquadric*> quad;

    std::map<GLuint, chemlib::MIAtom*> atomPickNameMap;
    std::map<GLuint, Annotation*> annotationPickNameMap;
    std::map<GLuint, chemlib::Bond*> bondPickNameMap;
    std::map<GLuint, chemlib::ANGLE*> anglePickNameMap;
    bool pickingEnabled;

    QGLWidget *qcanvas_win;

    int width;
    int height;
    float xyDisplayRange;
    std::auto_ptr<mi::opengl::Text> labelsFont_;
    std::auto_ptr<mi::opengl::Text> stackFont_;
    TargaImage *popStackImage;
    TargaImage *minimizeStackImage;
    TargaImage *unminimizeStackImage;
    TargaImage *closeStackImage;

    RenderStyle style;
    bool hideHydrogens;

    mi::opengl::Light *light;
    mi::opengl::Light *light2;

    mi::math::Vector3<float> viewVector;

    mi::opengl::Frustum *frustum;

    mi::opengl::Camera *camera;

    bool currentModel;

    bool fogEnabled;
    float fogStart;
    float fogEnd;
    bool antialiasLines;
    bool lineWidthDepthCued;
    bool dimNonactiveModels;
    float amountToDimNonactiveModels;

    bool joinBondsOfSameColor;

    bool showBondOrders;

    int fontSize;
    float labelsTextScale;

    bool viewVectorSet;

    void drawBondLine(const mi::math::Vector3<float> &pos1, int color1, const mi::math::Vector3<float> &pos2, int color2, float lineWidth);
    void drawBondCylinder(const mi::math::Vector3<float> &pos1, int color1, const mi::math::Vector3<float> &pos2, int color2, float radius, bool capped);

    /**
     * Draw a cylinder (used for the stick part of ball and stick).
     * @param start one end of the cylinder
     * @param end the other end
     * @param radius radius of the cylinder
     */
    static void DrawCylinder(const mi::math::Vector3<float> &start, const mi::math::Vector3<float> &end, float radius);

    static void normalize(float a[3]);

    static void cross(float a[3], float b[3], float c[3]);

    /**
     * Multiply a vector (1x3) times a matrix (3x3)
     * @param a vector
     * @param mat matrix
     * @param b transformed vector
     */
    static void vecmat(float a[3], float mat[3][3], float b[3]);

    static void setIdentity(float rotation[]);

    void DrawSphere(float c[3], float radius, int tess);
    static void DrawSphere2(float c[3], float radius, int tess);

    /**
     * Draw a sphereoid (sphere or ellipse)
     * @param c Center
     * @param radius axes of the ellipse
     * @param tess number of tesselations (<10 low, >10 high, negative wire frame)
     */
    static void RenderSphereoid(float c[3], float radius[3], int tess);

    static void Subdivide(int level, bool bWireFrame, float fnorm, float radius[3],
                          float n1[3], float n2[3], float n3[3], float v1[3], float v2[3], float v3[3]);

    static void TriangleRender(bool bWireFrame, float n1[3], float n2[3], float n3[3],
                               float v1[3], float v2[3], float v3[3]);

    void getInverseRotation(float inverseRotation[16]);

    void drawText(const char *text, float x, float y, float z, float offset, float scale, float inverseRotation[16], bool displayNumber, int number, float foreground[4]);

    GLuint getPickName(chemlib::MIAtom *atom);
    GLuint getPickName(Annotation &annotation);
    GLuint getPickName(chemlib::Bond &bond);
    GLuint getPickName(chemlib::ANGLE &angle);

    void drawStackText(int x, int y, const QString &s);

    void computeBounds(std::list<Molecule*> &molecules);

    void applyProjection(float scale);

    void drawResidueAtoms(chemlib::ResidueListIterator res, chemlib::ResidueListIterator resEnd);

    void DrawSecondaryStructure(SecondaryStructure *secondaryStructure);
    void DrawRibbonSegment(RibbonSegment *ribbonSegment);
    void DrawRibbonSpan(RibbonSpan *ribbonSpan);
    void DrawHelix(Helix *helix);
    void DrawPolySurf(int npr, int nseg, double dPoints[][3], double dNorms[][3]);
    void DrawCappedCylinder(int nSeg, double dAxisA[3], double dAxisB[3], double dEndA[][3], double dEndB[][3],
                            double dNormA[][3], double dNormB[][3], double dCapNormA[][3], double dCapNormB[][3]);

    void SetColor(unsigned char *bRGBA);

public:

    GLRenderer();
    virtual ~GLRenderer();

    void setQGLWidget(QGLWidget *qcanvas_win);

    void setViewVector(const mi::math::Vector3<float> &viewVector);

    void setFrustum(mi::opengl::Frustum *frustum);

    void setCamera(mi::opengl::Camera *camera);

    void swapBuffers();

    void Draw2(Displaylist *displaylist, bool ShowVus, bool drawMapsFirst, bool showSymmetryAsBackbone = false);
    void DrawLabels(std::list<Molecule*> &molecules);
    void DrawAnnotations(std::list<Molecule*> &molecules);
    void DrawStack(Stack *stack, int x, int y);
    void circleStackAtoms(Stack *stack);
    void CircleAtoms(chemlib::MIAtom **atoms, int natoms);
    void DrawMessage(const std::string&, int x, int y);
    void DrawContacts(std::vector<CONTACT> &Contacts);
    void DrawVus(std::vector<PLINE> &Vus, int w);
    void DrawMaps(std::vector<EMap*> &maps);

    void drawSurface(std::vector<SURFDOT> &dots);

    bool isFogEnabled();
    void setFogEnabled(bool on);

    float getFogStart();
    void setFogStart(float z);

    float getFogEnd();
    void setFogEnd(float z);

    void renderFog();

    bool isAntialiasLines();

    void setAntialiasLines(bool antialiasLines);

    bool isLineWidthDepthCued();

    void setLineWidthDepthCued(bool lineWithDepthCued);

    bool isDimNonactiveModels();

    void setDimNonactiveModels(bool dimNonactiveModels);
    void setAmountToDimNonactiveModels(float percent);

    int getFontSize();

    void setFontSize(int fontSize);

    RenderStyle&getRenderStyle();
    void setRenderStyle(const RenderStyle &style);

    void setHideHydrogens(bool on);
    bool isHideHydrogens();

    void updateTextScale(float glUnitsPerPixel);
    void initializeForRender();
    void renderLight();

    void drawBonds(std::vector<chemlib::Bond> &bonds, int *color = NULL);
    void drawAngles(std::vector<chemlib::ANGLE> &angles);

    void drawLines(std::vector<PLINE> &Vus, int w, bool withDepthTest);


    chemlib::MIAtom *getAtom(int id);
    chemlib::Bond *getBond(int id);
    chemlib::ANGLE *getAngle(int id);
    Annotation *getAnnotation(int id);

    void clearPickNames();

    float *getColor(int colorIndex, bool adjustForCurrentModel = false);
    QVector4D getColorVector(int colorIndex, bool adjustForCurrentModel);
    float *getColor(PaletteColor c, bool adjustForCurrentModel = false);

    void drawLabels(Molecule::AtomLabelList &labels);
    void drawAnnotations(Molecule::AnnotationList &Annotations, bool drawBox);
    void drawMolecule(Molecule *molecule);
    void drawSymmetryMolecule(Molecule *molecule, bool showSymmetryAsBackbone);

    void setJoinBondsOfSameColor(bool on);

    void setShowBondOrders(bool on);

    void drawText(const char *text, float x, float y, float z);

    void *getContext();

    void setPickingEnabled(bool value);
};

#endif // ifndef mifit_ui_GLRenderer_h
