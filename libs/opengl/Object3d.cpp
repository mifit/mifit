#include <opengl/Object3d.h>

#include <opengl/OpenGL.h>
#include <cstdio>

#include <stdlib.h>

using namespace std;

namespace mi
{
namespace opengl
{

Object3d::Object3d(istream &input, bool center)
    : numpolys(0),
      toppoint(0.0f),
      bottompoint(0.0f),
      leftpoint(0.0f),
      rightpoint(0.0f),
      farpoint(0.0f),
      nearpoint(0.0f)
{

    loadobject(input);
    if (center)
    {
        centerObject();
    }
    createList();
    numpolys = faces.size();
    cleanup();
}

void Object3d::cleanup()
{
    vertexsets.clear();
    vertexsetsnorms.clear();
    vertexsetstexs.clear();
    faces.clear();
    facestexs.clear();
    facesnorms.clear();
}

void stringSplit(string str, string delim, vector<string> &results)
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

string stringReplaceAll(string s, string target, string replacement)
{
    string result(s);
    unsigned int pos;
    while ((pos = result.find(target)) != string::npos)
    {
        result.replace(pos, target.size(), replacement);
    }
    return result;
}

string stringTrim(string &s, const string &drop = " ")
{
    std::string r = s.erase(s.find_last_not_of(drop)+1);
    return r.erase(0, r.find_first_not_of(drop));
}

void Object3d::loadobject(istream &br)
{
    int linecounter = 0;

    string newline;
    bool firstpass = true;
    vector<string> currentGroups;

    getline(br, newline);
    while (br)
    {
        linecounter++;
        newline = stringTrim(newline);
        if (newline.length() > 0)
        {
            if (newline[0] == 'v' && newline[1] == ' ')
            {
                vector<float> coords;
                vector<string> coordstext;
                stringSplit(newline, " ", coordstext);
                for (unsigned int i = 1; i < coordstext.size(); i++)
                {
                    coords.push_back((float)atof(coordstext[i].c_str()));
                }
                if (firstpass)
                {
                    rightpoint = coords[0];
                    leftpoint = coords[0];
                    toppoint = coords[1];
                    bottompoint = coords[1];
                    nearpoint = coords[2];
                    farpoint = coords[2];
                    firstpass = false;
                }
                if (coords[0] > rightpoint)
                {
                    rightpoint = coords[0];
                }
                if (coords[0] < leftpoint)
                {
                    leftpoint = coords[0];
                }
                if (coords[1] > toppoint)
                {
                    toppoint = coords[1];
                }
                if (coords[1] < bottompoint)
                {
                    bottompoint = coords[1];
                }
                if (coords[2] > nearpoint)
                {
                    nearpoint = coords[2];
                }
                if (coords[2] < farpoint)
                {
                    farpoint = coords[2];
                }
                vertexsets.push_back(coords);
            }
            if (newline[0] == 'v' && newline[1] == 't')
            {
                vector<float> coords;
                vector<string> coordstext;
                stringSplit(newline, " ", coordstext);
                for (unsigned int i = 1; i < coordstext.size(); i++)
                {
                    coords.push_back((float)atof(coordstext[i].c_str()));
                }
                vertexsetstexs.push_back(coords);
            }
            if (newline[0] == 'v' && newline[1] == 'n')
            {
                vector<float> coords;
                vector<string> coordstext;
                stringSplit(newline, " ", coordstext);
                for (unsigned int i = 1; i < coordstext.size(); i++)
                {
                    coords.push_back((float)atof(coordstext[i].c_str()));
                }
                vertexsetsnorms.push_back(coords);
            }
            if (newline[0] == 'f' && newline[1] == ' ')
            {
                vector<string> coordstext;
                stringSplit(newline, " ", coordstext);
                vector<int> v;
                vector<int> vt;
                vector<int> vn;

                for (unsigned int i = 1; i < coordstext.size(); i++)
                {
                    string fixstring = stringReplaceAll(coordstext[i], "//", "/0/");
                    vector<string> tempstring;
                    stringSplit(fixstring, "/", tempstring);
                    v.push_back(atoi(tempstring[0].c_str()));
                    if (tempstring.size() > 1)
                    {
                        vt.push_back(atoi(tempstring[1].c_str()));
                    }
                    else
                    {
                        vt.push_back(0);
                    }
                    if (tempstring.size() > 2)
                    {
                        vn.push_back(atoi(tempstring[2].c_str()));
                    }
                    else
                    {
                        vn.push_back(0);
                    }
                }
                faces.push_back(v);
                facesGroups.push_back(currentGroups);
                facestexs.push_back(vt);
                facesnorms.push_back(vn);
            }
            if (newline[0] == 'g' && newline[1] == ' ')
            {
                vector<string> groupsText;
                stringSplit(newline, " ", groupsText);
                currentGroups.clear();
                vector<string>::iterator i = groupsText.begin();
                ++i;
                for (; i != groupsText.end(); ++i)
                {
                    currentGroups.push_back(*i);
                }
            }
        }
        getline(br, newline);
    }
}

void Object3d::centerObject()
{
    float xshift = (rightpoint - leftpoint) / 2.0f;
    float yshift = (toppoint - bottompoint) / 2.0f;
    float zshift = (nearpoint - farpoint) / 2.0f;

    for (unsigned int i = 0; i < vertexsets.size(); i++)
    {
        vector<float> coords;

        coords.push_back(vertexsets[i][0] - leftpoint - xshift);
        coords.push_back(vertexsets[i][1] - bottompoint - yshift);
        coords.push_back(vertexsets[i][2] - farpoint - zshift);

        vertexsets[i] = coords;
    }

}

float Object3d::getXWidth()
{
    float returnval = 0.0f;
    returnval = rightpoint - leftpoint;
    return returnval;
}

float Object3d::getYHeight()
{
    float returnval = 0.0f;
    returnval = toppoint - bottompoint;
    return returnval;
}

float Object3d::getZDepth()
{
    float returnval = 0.0f;
    returnval = nearpoint - farpoint;
    return returnval;
}

int Object3d::numpolygons()
{
    return numpolys;
}

void Object3d::createList()
{

    this->objectlist = glGenLists(1);

    glNewList(objectlist, GL_COMPILE);
    for (unsigned int i = 0; i < faces.size(); i++)
    {
        vector<int> &tempfaces = faces[i];
        vector<int> &tempfacesnorms = facesnorms[i];
        vector<int> &tempfacestexs = facestexs[i];

        int polytype;
        if (tempfaces.size() == 3)
        {
            polytype = GL_TRIANGLES;
        }
        else if (tempfaces.size() == 4)
        {
            polytype = GL_QUADS;
        }
        else
        {
            polytype = GL_POLYGON;
        }
        glPushName(i);
        glBegin(polytype);

        for (unsigned int w = 0; w < tempfaces.size(); w++)
        {
            if (tempfacesnorms[w] != 0)
            {
                float normtempx = vertexsetsnorms[tempfacesnorms[w] - 1][0];
                float normtempy = vertexsetsnorms[tempfacesnorms[w] - 1][1];
                float normtempz = vertexsetsnorms[tempfacesnorms[w] - 1][2];
                glNormal3f(normtempx, normtempy, normtempz);
            }

            if (tempfacestexs[w] != 0)
            {
                float textempx = vertexsetstexs[tempfacestexs[w] - 1][0];
                float textempy = vertexsetstexs[tempfacestexs[w] - 1][1];
                float textempz = vertexsetstexs[tempfacestexs[w] - 1][2];
                glTexCoord3f(textempx, 1.0f - textempy, textempz);
            }

            float tempx = vertexsets[tempfaces[w] - 1][0];
            float tempy = vertexsets[tempfaces[w] - 1][1];
            float tempz = vertexsets[tempfaces[w] - 1][2];
            glVertex3f(tempx, tempy, tempz);
        }

        glEnd();
        glPopName();

    }
    glEndList();
}

void Object3d::render()
{
    glCallList(objectlist);
}

vector<string>&Object3d::getGroups(int index)
{
    static vector<string> emptyVector;
    if (index < 0 || index > (int)facesGroups.size())
    {
        return emptyVector;
    }
    return facesGroups[index];
}

}
}

