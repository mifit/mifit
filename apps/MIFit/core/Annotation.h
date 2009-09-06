#ifndef ANNOTATION_H_
#define ANNOTATION_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <QObject>
#include "xmlarchive.h"

extern char ann_typenames[][25];

/**
 * Annotations function like yellow stickies on the model.
 * They get saved with the session and can be used to add arbitrary information
 * in a positional way.
 */
class Annotation : public QObject {
  Q_OBJECT

  static const float NO_COORD;

  bool hidden;


public:
  //@{
  // Enumeration of different annotation types
  //@}
  enum
  {
    Unknown,
    Note_to_self,
    Commands,
    Structure_annotation,
    Discussion_thread,
    Geom_error,
    URL,
    Plane
  };
  //@{
  // The type - one of the enumeration values
  //@}
  unsigned int m_type;   // one of enum
  //@{
  // An id for this annotation
  //@}
  unsigned long m_id;
  //@{
  // Coordinates in Cartesian space
  //@}
  float m_x, m_y, m_z;
  //@{
  // sx, sy, sz Coordinates in screen space.
  // sw, sh
  //  yadd is the offset used to prevent overlaps.
  //@}
  int sx, sy, sz, sw, sh, yadd;
  //@{
  // The text of the annotation
  //@}
  std::string m_text;
  void setText(const std::string& text);
  //@{
  // The color of the annotation.
  //@}
  PaletteColor m_color;
  void setColor(PaletteColor color);
  //@{
  // The guy who wrote this annotation
  //@}
  std::string m_author;
  //@{
  // A subject fo this annotation for classification purposes
  //@}
  std::string m_subject;

  bool isHidden();
  void setHidden(bool hidden);

  //@{
  // The constructor.  Only the text is required.
  //@}
  Annotation(const char* text, float x = NO_COORD, float y = NO_COORD, float z = NO_COORD, int id = 0);
  //@{
  // Default constructor.  If you use this you must set all the values yourself.
  //@}
  Annotation();
  //@{
  // Destructor.
  //@}
  virtual ~Annotation();
  //@{
  // Return the x coordinate.
  //@}
  float GetX() {
    return m_x;
  }

  //@{
  // Return the y coordinate.
  //@}
  float GetY() {
    return m_y;
  }

  //@{
  // Return the z coordinate.
  //@}
  float GetZ() {
    return m_z;
  }

  //@{
  // Return the text string.
  //@}
  const char* GetText() {
    return (const char*)m_text.c_str();
  }

  //@{
  // Change the position of the annotation
  //@}
  void SetPosition(float x, float y, float z) {
    m_x = x; m_y = y; m_z = z;
  }

  //@{
  // Change the text string
  //@}
  void SetText(const char* t) {
    m_text = t;
  }

  //@{
  // Write the annotation to an XML file
  // @param ar the XMLArchive to write to.  See also XMLArchive.
  //@}
  void WriteXML(XMLArchive* ar);
  //@{
  // Clears all the fields and resets to the default values.
  //@}
  void Clear() {
    m_x = m_y = m_z = 0;
    m_color = PaletteColor(255, 255, 255);
    m_text.clear();
    m_type = Annotation::Unknown;
    m_id = 0;
    m_author.clear();
    m_subject.clear();
  }

signals:
  void annotationChanged(Annotation*);

};

#endif
