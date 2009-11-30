#ifndef mi_opengl_Object3d_h
#define mi_opengl_Object3d_h

#include <vector>
#include <string>
#include <iostream>

namespace mi
{
    namespace opengl
    {

        class Object3d
        {

            std::vector<std::vector<float> > vertexsets;
            std::vector<std::vector<float> > vertexsetsnorms;
            std::vector<std::vector<float> > vertexsetstexs;
            std::vector<std::vector<int> > faces;

            std::vector<std::vector<std::string> > facesGroups;
            std::vector<std::vector<int> > facestexs;
            std::vector<std::vector<int> > facesnorms;

            int objectlist;

            int numpolys;

            void cleanup();

            void loadobject(std::istream &input);

            void centerObject();

        public:

            float toppoint;
            float bottompoint;
            float leftpoint;
            float rightpoint;
            float farpoint;
            float nearpoint;

            Object3d(std::istream &input, bool centerit);

            float getXWidth();

            float getYHeight();

            float getZDepth();

            int numpolygons();

            void createList();

            void render();

            std::vector<std::string>&getGroups(int index);

        };

    }
}

#endif // ifndef mi_opengl_Object3d_h

