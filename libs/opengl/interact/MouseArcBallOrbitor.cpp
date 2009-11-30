#include <opengl/interact/MouseArcBallOrbitor.h>

#include <opengl/Viewpoint.h>
#include <math/Vector3.h>

using namespace mi::math;
using namespace mi::opengl;

namespace mi
{
namespace opengl
{
namespace interact
{

MouseArcBallOrbitor::MouseArcBallOrbitor(Viewpoint *viewpoint, float distanceToTarget)
    : MouseArcBall(viewpoint),
      startQuat(0.0f, 0.0f, 0.0f, 1.0f)
{

    this->distanceToTarget = distanceToTarget;
}

MouseArcBallOrbitor::MouseArcBallOrbitor(Viewpoint *viewpoint, Point3<float> &center, float radius, float distanceToTarget)
    : MouseArcBall(viewpoint, center, radius),
      startQuat(0.0f, 0.0f, 0.0f, 1.0f)
{

    this->distanceToTarget = distanceToTarget;
}

void MouseArcBallOrbitor::initializeRotation()
{
    startQuat.set(viewpoint->getRotation());
    historicalQuat.set(0.0f, 0.0f, 0.0f, 1.0f);
    currentQuat.set(historicalQuat);
}

void MouseArcBallOrbitor::applyRotation()
{
    Quaternion<float> q(startQuat);
    q.multiply(currentQuat);
    Vector3<float> view = viewpoint->getTarget(-distanceToTarget);
    viewpoint->setRotation(q);
    Vector3<float> view2 = viewpoint->getTarget(-distanceToTarget);
    view2.subtract(view);
    viewpoint->translate(view2);
}

void MouseArcBallOrbitor::setDistanceToTarget(float distanceToTarget)
{
    this->distanceToTarget = distanceToTarget;
}

}
}
}
