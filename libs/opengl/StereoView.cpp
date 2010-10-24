#include <opengl/StereoView.h>

#include <math/Vector3.h>
#include <opengl/Camera.h>
#include <opengl/Frustum.h>
#include <opengl/Viewport.h>
#include <opengl/Scene.h>
#include <opengl/OpenGL.h>
#include <opengl/QuatUtil.h>

using namespace mi::math;

#ifdef _DEBUG
static void printOpenGLErrors(const char *file, int line)
{
    int error = glGetError();
    if (error != 0)
    {
        printf("OpenGL error %s (%d) at line %d in %s\n", gluErrorString(error), error, line, file);
        fflush(stdout);
    }
}

#else
static void printOpenGLErrors(const char*, int)
{
}

#endif

namespace mi
{
namespace opengl
{

StereoView::StereoView(Frustum *frustum, Camera *camera)
{
    viewport = new Viewport();
    rightViewport = new Viewport();
    leftViewport = new Viewport();
    this->frustum = frustum;
    this->camera = camera;
    stereoCamera = new Camera();
    overlayFrustum = new Frustum();
    overlayFrustum->setPerspective(false);
    overlayFrustum->setNearClipping(-1000.0f);
    overlayFrustum->setFarClipping(1000.0f);
    overlayCamera = new Camera();
    eyeSeparation = 1.0f;
    imageSeparation = 20.0f;

    stereo = false;
    hardwareStereo = false;
    glStereo = false;
}

StereoView::~StereoView()
{
    delete viewport;
    delete stereoCamera;
    delete rightViewport;
    delete leftViewport;
    delete overlayFrustum;
    delete overlayCamera;
}

bool StereoView::isGLStereo()
{
    return glStereo;
}

bool StereoView::isStereo()
{
    return stereo;
}

void StereoView::setStereo(bool stereo)
{
    this->stereo = stereo;
}

bool StereoView::isHardwareStereo()
{
    return hardwareStereo && glStereo;
}

void StereoView::setHardwareStereo(bool hardwareStereo)
{
    this->hardwareStereo = hardwareStereo;
}

Viewport*StereoView::getViewport()
{
    return viewport;
}

Viewport*StereoView::getLeftViewport()
{
    return leftViewport;
}

Viewport*StereoView::getRightViewport()
{
    return rightViewport;
}

void StereoView::render(Scene &scene)
{
    printOpenGLErrors(__FILE__, __LINE__);
    GLboolean b;
    glGetBooleanv(GL_STEREO, &b);
    glStereo = (b == GL_TRUE);
    printOpenGLErrors(__FILE__, __LINE__);

    frustum->updateFrustum(*viewport);
    overlayFrustum->setFocalLength(frustum->getFocalLength());
    overlayFrustum->setHeight((float)viewport->getHeight());

    resetView(GL_DEPTH_BUFFER_BIT);

    Viewport *currentViewport = viewport;
    if (stereo)
    {
        stereoCamera->setRotation(camera->getRotation());
        float cameraOffset;
        float frustumAsymmetry;
        calcStereoParameters(&cameraOffset, &frustumAsymmetry);

        currentViewport = prepareRightEye(cameraOffset, frustumAsymmetry);
        scene.setViewport(currentViewport);

        scene.preCameraRender();

        stereoCamera->render();
        frustum->updatePlanes();
        scene.render();

        currentViewport = prepareLeftEye(cameraOffset, frustumAsymmetry);
        scene.setViewport(currentViewport);

        stereoCamera->render();
        frustum->updatePlanes();
        scene.render();

        resetView(GL_DEPTH_BUFFER_BIT);

        if (isHardwareStereo())
        {
            glDrawBuffer(GL_BACK_RIGHT);
            scene.setViewport(viewport);
            viewport->render();
            overlayFrustum->render(*viewport);
        }
        else
        {
            scene.setViewport(rightViewport);
            rightViewport->render();
            overlayFrustum->render(*rightViewport);
        }
        scene.renderOverlay();

        if (isHardwareStereo())
        {
            glDrawBuffer(GL_BACK_LEFT);
            scene.setViewport(viewport);
            viewport->render();
            overlayFrustum->render(*viewport);
        }
        else
        {
            scene.setViewport(leftViewport);
            leftViewport->render();
            overlayFrustum->render(*leftViewport);
        }
        scene.renderOverlay();

    }
    else
    {
        scene.setViewport(viewport);
        viewport->render();
        frustum->render(*viewport);

        scene.preCameraRender();
        camera->render();
        frustum->updatePlanes();
        scene.render();

        resetView(GL_DEPTH_BUFFER_BIT);

        overlayFrustum->render(*viewport);
        scene.renderOverlay();

    }
}

void StereoView::resetView(int flags)
{
    printOpenGLErrors(__FILE__, __LINE__);
    if (glStereo)
    {
        glDrawBuffer(GL_BACK_RIGHT);
        glClear(flags);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glDrawBuffer(GL_BACK_LEFT);
        glClear(flags);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    else
    {
        glClear(flags);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    printOpenGLErrors(__FILE__, __LINE__);
}

Viewport*StereoView::prepareRightEye(float cameraOffset, float frustumAsymmetry)
{
    printOpenGLErrors(__FILE__, __LINE__);

    if (glStereo)
    {
        glDrawBuffer(GL_BACK_RIGHT);
        glClear(GL_DEPTH_BUFFER_BIT);
    }

    Vector3<float> newEye(1.0f, 0.0f, 0.0f);
    QuatUtil::rotateVector(camera->getRotation(), newEye);
    newEye.scale(cameraOffset);
    newEye.add(camera->getEye());
    stereoCamera->setEye(newEye);

    Viewport *currentViewport;
    if (isHardwareStereo())
    {
        glDrawBuffer(GL_BACK_RIGHT);
        printOpenGLErrors(__FILE__, __LINE__);
        currentViewport = viewport;
        viewport->render();
        frustum->render(*viewport, frustumAsymmetry);
    }
    else
    {
        currentViewport = rightViewport;
        rightViewport->render();
        frustum->render(*rightViewport, frustumAsymmetry);
    }
    printOpenGLErrors(__FILE__, __LINE__);
    return currentViewport;
}

Viewport*StereoView::prepareLeftEye(float cameraOffset, float frustumAsymmetry)
{
    printOpenGLErrors(__FILE__, __LINE__);

    if (glStereo)
    {
        glDrawBuffer(GL_BACK_LEFT);
        glClear(GL_DEPTH_BUFFER_BIT);
    }

    Vector3<float> newEye(1.0f, 0.0f, 0.0f);
    QuatUtil::rotateVector(camera->getRotation(), newEye);
    newEye.scale(-cameraOffset);
    newEye.add(camera->getEye());
    stereoCamera->setEye(newEye);

    Viewport *currentViewport;
    if (isHardwareStereo())
    {
        glDrawBuffer(GL_BACK_LEFT);
        printOpenGLErrors(__FILE__, __LINE__);
        currentViewport = viewport;
        viewport->render();
        frustum->render(*viewport, -frustumAsymmetry);
    }
    else
    {
        currentViewport = leftViewport;
        leftViewport->render();
        frustum->render(*leftViewport, -frustumAsymmetry);
    }
    return currentViewport;
}

void StereoView::render(Renderable &renderable)
{
    frustum->updateFrustum(*viewport);

    resetView(GL_DEPTH_BUFFER_BIT);

    if (stereo)
    {
        stereoCamera->setRotation(camera->getRotation());
        float cameraOffset;
        float frustumAsymmetry;
        calcStereoParameters(&cameraOffset, &frustumAsymmetry);

        prepareRightEye(cameraOffset, frustumAsymmetry);

        stereoCamera->render();
        frustum->updatePlanes();
        renderable.render();

        prepareLeftEye(cameraOffset, frustumAsymmetry);

        stereoCamera->render();
        frustum->updatePlanes();
        renderable.render();

    }
    else
    {
        viewport->render();
        frustum->render(*viewport);

        camera->render();
        frustum->updatePlanes();
        renderable.render();

    }

}

Frustum*StereoView::getFrustum()
{
    return frustum;
}

void StereoView::setFrustum(Frustum *frustum)
{
    this->frustum = frustum;
}

Camera*StereoView::getCamera()
{
    return camera;
}

void StereoView::setCamera(Camera *camera)
{
    this->camera = camera;
}

void StereoView::setSize(int width, int height)
{
    viewport->set(0, 0, width, height);
    rightViewport->set(0, 0, width / 2, height);
    leftViewport->set(width / 2, 0, width / 2, height);
}

bool StereoView::isCrossEyed()
{
    return crossEyed;
}

void StereoView::setCrossEyed(bool crossEyed)
{
    this->crossEyed = crossEyed;
}

float StereoView::getEyeSeparation()
{
    return eyeSeparation;
}

void StereoView::setEyeSeparation(float eyeSeparation)
{
    this->eyeSeparation = eyeSeparation;
}

float StereoView::getImageSeparation()
{
    return imageSeparation;
}

void StereoView::setImageSeparation(float imageSeparation)
{
    this->imageSeparation = imageSeparation;
}

void StereoView::calcStereoParameters(float *cameraOffset, float *frustumAsymmetry)
{

    *cameraOffset = -frustum->getPerspectiveHeight(frustum->getFocalLength());
    *cameraOffset *= 0.025f * eyeSeparation;
    *frustumAsymmetry = *cameraOffset * frustum->getNearClipping() / frustum->getFocalLength();
    *frustumAsymmetry -= imageSeparation * 0.10f * *frustumAsymmetry;

    if (crossEyed)
    {
        *frustumAsymmetry *= -1.0;
        *cameraOffset *= -1.0f;
    }

}

}

}
