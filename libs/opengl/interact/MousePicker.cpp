#include <opengl/interact/MousePicker.h>

#include <opengl/Frustum.h>
#include <opengl/Renderable.h>

#include <climits>

namespace mi
{
namespace opengl
{
namespace interact
{

MousePicker::MousePicker()
{
    selectBufferLength = 1024;
    selectBuffer = new GLuint[selectBufferLength];
}

MousePicker::~MousePicker()
{
    delete[] selectBuffer;
}

std::vector<GLuint> MousePicker::pick(int x, int y, Frustum *frustum, Renderable *scene)
{

    glSelectBuffer(selectBufferLength, selectBuffer);
    glRenderMode(GL_SELECT);

    glInitNames();

    frustum->beginPicking(x, y);

    scene->render();

    frustum->endPicking();

    int numberOfHits = glRenderMode(GL_RENDER);

    std::vector<GLuint> selectedIds;
    if (numberOfHits > 0)
    {
        GLuint smallestMinZ = UINT_MAX;
        int i = 0;
        for (int index = 0; index < numberOfHits; ++index)
        {
            GLuint numberOfNames = selectBuffer[i];
            GLuint minZ = selectBuffer[i + 1];
            // GLuint maxZ = selectBuffer[i+2];
            if (minZ < smallestMinZ)
            {
                smallestMinZ = minZ;
                selectedIds.clear();
                for (GLuint j = 0; j < numberOfNames; ++j)
                {
                    GLuint name = selectBuffer[i + 3 + j];
                    selectedIds.push_back(name);
                }
            }
            i += numberOfNames + 3;
        }
    }

    return selectedIds;
}

}

}
}
