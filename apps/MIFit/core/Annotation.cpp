#include "Annotation.h"

using namespace chemlib;

const float Annotation::NO_COORD = -99999.9f;

char ann_typenames[][25] =
{
    "Unknown",
    "Note_to_self",
    "Commands",
    "Structure_annotation",
    "Discussion_thread",
    "URL",
    "Plane"
};
Annotation::Annotation(const char *text, float x, float y, float z, int id)
{
    SetPosition(x, y, z);
    SetText(text);
    PaletteColor white(255, 255, 255);
    m_color = white;
    m_type = Unknown;
    m_id = id;
    hidden = false;
}

Annotation::Annotation()
{
    m_x = m_y = m_z = 0;
    PaletteColor white(255, 255, 255);
    m_color = white;
    m_type = Unknown;
    m_id = 0;
    hidden = false;
}

Annotation::~Annotation()
{

}

void Annotation::WriteXML(XMLArchive *ar)
{
    std::string s = format("Id=\"%d\" x=\"%f\" y=\"%f\" z=\"%f\"", (int)m_id, m_x, m_y, m_z);
    ar->BeginTag("Annotation", s.c_str());
    ar->WriteField("Type", ann_typenames[m_type]);
    ar->WriteField("Text", m_text.c_str());
    ar->WriteColor(m_color);
    ar->EndTag();
}

bool Annotation::isHidden()
{
    return hidden;
}

void Annotation::setHidden(bool hidden)
{
    this->hidden = hidden;
    annotationChanged(this);
}

void Annotation::setText(const std::string &text)
{
    m_text = text;
    annotationChanged(this);
}

void Annotation::setColor(PaletteColor color)
{
    m_color = color;
    annotationChanged(this);
}

