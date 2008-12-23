//TODO: move to ui directory?

#include "GLDrawingCanvas.h"
#include "../ui/id.h"
#include "../ui/MIHistory.h"
#include "../wxdr/MIDialog.h"

#include "../ui/MIMenuBar.h"
#include "../ui/MIMenu.h"
#include "../ui/MIEventHandlerMacros.h"

#include <QPrinter>
#include <QPrintDialog>
#include <QClipboard>
#include <QApplication>
#include <QPainter>
#include <QDir>

const unsigned int NSPOKES=9;
const unsigned int SPOKESPACING=18;
//const unsigned int SPOKELENGTH=8;

enum
{
  ID_ASP_PRINT = 800,
  ID_ASP_EDITCOPY = 801,
  ID_ASP_ADD_ANNOTATION = 803,
  ID_ASP_EXPORT_IMAGE = 804,
  ID_ASP_CURRENT_ITEM_FONT = 805,
  ID_ASP_ALL_FONTS = 806,
  ID_ASP_ANNOTATION_FONTS = 807,
  ID_ASP_ATOM_FONTS = 808,
  ID_ASP_RESIDUE_FONTS = 809,
  ID_ASP_HYDROGEN_BOND_FONTS = 810,
  ID_ASP_SET_CURRENT_FONT = 811
};



using namespace moldraw;

GLDrawingCanvas::GLDrawingCanvas(QWidget *parent) :
  QGLWidget(parent), MIEventHandler(this),
  orig_x(0.0f), orig_y(0.0f), dim_x(1.0f), dim_y(1.0f), _selected_item(UINT_MAX) {

  menuBar = new MIMenuBar(new QMenuBar(parent));
  MIGetHistory()->AddMenuBar(menuBar, "ActiveSitePlot");

  MIMenu* menuFile = new MIMenu(*this,this);
  MIMenu* menuEdit = new MIMenu(*this,this);
  MIMenu* menuFont = new MIMenu(*this,this);
  //MIMenu* menuHelp = new MIMenu(*this,this);

  menuFile->Append(ID_ASP_PRINT, "Print...", "Print the canvas");
  menuFile->Append(ID_ASP_EXPORT_IMAGE, "Export Image", "Export a PNG file");
  //menuFile->AppendSeparator();
  menuBar->Append(menuFile, "&File");

  menuEdit->Append(ID_ASP_EDITCOPY, "&Copy");
  menuEdit->Append(ID_ASP_ADD_ANNOTATION, "&Add Annotation");
  menuEdit->Append(ID_ASP_SET_CURRENT_FONT, "&Select Default Font");
  menuBar->Append(menuEdit, "&Edit");

  menuFont->Append(ID_ASP_CURRENT_ITEM_FONT, "Change Current Item's Font");
  menuFont->Append(ID_ASP_ALL_FONTS, "Change All Fonts");
  menuFont->Append(ID_ASP_ANNOTATION_FONTS, "Change All Annotation Fonts");
  menuFont->Append(ID_ASP_ATOM_FONTS, "Change Atom Label Font");
  menuFont->Append(ID_ASP_RESIDUE_FONTS, "Change Residue Label Font");
  menuFont->Append(ID_ASP_HYDROGEN_BOND_FONTS, "Change Hydrogen Bond Font");
  menuBar->Append(menuFont, "&Font");


  BEGIN_EVENT_TABLE(this, none);
  EVT_MENU(ID_ASP_EDITCOPY, ASPView::OnEditCopy);
  EVT_MENU(ID_ASP_PRINT, ASPView::OnPrint);
  EVT_MENU(ID_ASP_EXPORT_IMAGE, ASPView::OnExportImage);

  EVT_MENU(ID_ASP_ADD_ANNOTATION, ASPView::OnAddAnnotation);
  EVT_MENU(ID_ASP_SET_CURRENT_FONT, ASPView::OnSelectFont);


  EVT_MENU(ID_ASP_CURRENT_ITEM_FONT,   ASPView::OnChangeFonts);
  EVT_MENU(ID_ASP_ALL_FONTS,           ASPView::OnChangeFonts);
  EVT_MENU(ID_ASP_ANNOTATION_FONTS,    ASPView::OnChangeFonts);
  EVT_MENU(ID_ASP_ATOM_FONTS,          ASPView::OnChangeFonts);
  EVT_MENU(ID_ASP_RESIDUE_FONTS,       ASPView::OnChangeFonts);
  EVT_MENU(ID_ASP_HYDROGEN_BOND_FONTS, ASPView::OnChangeFonts);
  END_EVENT_TABLE();
}

GLDrawingCanvas::~GLDrawingCanvas() {
  MIGetHistory()->RemoveMenuBar(menuBar);
}


void GLDrawingCanvas::initializeGL() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
}


void GLDrawingCanvas::resizeGL(int width, int height) {
  glViewport(0,0,(GLint)width, (GLint)height);

  float haspect=(float)width/(float)height;
  haspect=1.0f;


  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (haspect < 1.0) {
    glOrtho(orig_x, (orig_x + dim_x), orig_y/haspect, (orig_y + dim_y)/haspect, -1.0f, 1.0f);
  } else {
    glOrtho(orig_x * haspect, (orig_x + dim_x) * haspect, orig_y, (orig_y + dim_y), -1.0f, 1.0f);
  }
}

void GLDrawingCanvas::toWorldCoords(const QPoint &pos, float &x, float &y)
{
  float fx=(float)pos.x()/(float)width();
  float fy=1.0f-(float)pos.y()/(float)height();
  float haspect=(float)width()/(float)height();
  haspect=1.0f;

  if (haspect < 1.0 ) {
    x=orig_x + fx * dim_x;
    y=orig_y + fy * dim_y/haspect;
  } else {
    x=orig_x + fx * dim_x * haspect;
    y=orig_y + fy * dim_y;
  }
}



static void drawCircle (float x, float y, float r) {
  float cx, cy;
  glBegin(GL_POLYGON);
  for (int i = 0; i <= 360; i+=18)
  {
    cx = r * cos(i*DEG2RAD);
    cy = r * sin(i*DEG2RAD);
    glVertex3f(cx + x, cy + y, 0.0f);
  }
  glEnd();
}

static void DrawSpokesGL(float px, float py, float radius, float linewidth, int direction)
{
  float spokelength = radius * 0.4F;
  //Direction (in radians) of first spoke...other spokes proceed counterclockwise
  int theta;
  int init_theta = (direction - SPOKESPACING * (NSPOKES - 1) / 2);
  int end_theta = init_theta + SPOKESPACING * (NSPOKES - 1);

  float old_startx=0.0f;
  float old_starty=0.0f;

  glLineWidth(linewidth);
  for (theta = init_theta; theta <= end_theta; theta += SPOKESPACING) {
    float vx = cos(DEG2RAD * theta);
    float vy = sin(DEG2RAD * theta);

    float startx=radius * vx;
    float starty=radius * vy;
    float endx=(radius + spokelength) * vx;
    float endy=(radius + spokelength) * vy;

    glBegin(GL_LINE_STRIP);

    // ray
    glVertex2f(px + endx,   py + endy);
    glVertex2f(px + startx, py + starty);

    // arc segment
    if (theta != init_theta) {
      glVertex2f(px + old_startx, py + old_starty);
    }
    old_startx=startx;
    old_starty=starty;

    glEnd();
  }
}

//taken from Qt 4.2.3 to avoid assertion failure with faulty glx driver(s) from xorg 6.9/7.0 on linux.
void GLDrawingCanvas::myRenderText(int x, int y, const QString & str, const QFont & font, int listBase)
{
#ifdef __linux__
  static bool use_alternate_rendering=false;
  static bool firsttime=true;
  if (firsttime) {
    firsttime=false;
    use_alternate_rendering=
      MIConfig::Instance()->GetProfileInt("View Parameters", "AlternateTextRendering", 0) != 0;
  }

  if (use_alternate_rendering) {
    makeCurrent();
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, width(), height(), 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glAlphaFunc(GL_GREATER, 0.0);
    glEnable(GL_ALPHA_TEST);
    glRasterPos2i(0, 0);
    glBitmap(0, 0, 0, 0, x, -y, NULL);
    glListBase(fontDisplayListBase(font, listBase));
    QByteArray cstr(str.toLatin1());
    glCallLists(cstr.size(), GL_UNSIGNED_BYTE, cstr.constData());

    // can't support this, but not needed
    //if (font.underline() || font.strikeOut() || font.overline())
    // qt_drawFontLining(x, y, str, font);

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glPopAttrib();
    return;
  }
#endif
  renderText(x, y, str, font, listBase);
}

//taken from Qt 4.2.3 to avoid assertion failure with faulty glx driver(s) from xorg 6.9/7.0 on linux.
void GLDrawingCanvas::myRenderText(double x, double y, double z, const QString & str, const QFont & font,
                           int listBase)
{
#ifdef __linux__
  static bool use_alternate_rendering=false;
  static bool firsttime=true;
  if (firsttime) {
    firsttime=false;
    use_alternate_rendering=
      MIConfig::Instance()->GetProfileInt("View Parameters", "AlternateTextRendering", 0) != 0;
  }

  if (use_alternate_rendering) {
    makeCurrent();
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_CULL_FACE);

    glRasterPos3d(x, y, z);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glAlphaFunc(GL_GREATER, 0.0);
    glEnable(GL_ALPHA_TEST);
    glListBase(fontDisplayListBase(font, listBase));
    QByteArray cstr(str.toLatin1());
    glCallLists(cstr.size(), GL_UNSIGNED_BYTE, cstr.constData());

    glPopAttrib();
    return;
  }
#endif
  renderText(x, y, z, str, font, listBase);
}



void GLDrawingCanvas::paintGL() {

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();


  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColor3f(0.0f, 0.0f, 0.0f);
  glLineWidth(1.0);

  for (size_t i=0; i<_items.size(); ++i) {
    draw_item &item=_items[i];

    //FIXME: push name on stack for picking
    glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_ENABLE_BIT);

    switch (item.type) {

      case draw_item::Atom:
        glDisable(GL_DEPTH_TEST);
        glColor3ub((GLubyte)item.color.red,
                   (GLubyte)item.color.green,
                   (GLubyte)item.color.blue);
        drawCircle(item.x1, item.y1, item.r);
        break;

      case draw_item::SingleBond:
          glLineWidth(item.width);
          glBegin(GL_LINES);
          glVertex2f(item.x1, item.y1);
          glVertex2f(item.x2, item.y2);
          glEnd();
        break;
      case draw_item::DoubleBond:
        //FIXME: render two lines here?
        glLineWidth(item.width);
        glBegin(GL_LINES);
        glVertex2f(item.x1, item.y1);
        glVertex2f(item.x2, item.y2);
        glEnd();
        break;

      case draw_item::HydrogenBond:
         glLineStipple(4,0xAAAA);
         glEnable(GL_LINE_STIPPLE);
         glLineWidth(item.width);
         glBegin(GL_LINES);
         glVertex2f(item.x1, item.y1);
         glVertex2f(item.x2, item.y2);
         glEnd();
        break;

      case draw_item::Sun:
      {
        QFont f=QFont();
        if (item.font.size()) {
          f.fromString(item.font.c_str());
        } else {
          //FIXME: set from prefs?
          f.setPointSizeF(item.font_size);
        }

        myRenderText(item.x1, item.y1, 0.0f, item.text.c_str(),f);
        break;
      }

      case draw_item::Spokes:
      {
        float exact_direction = RAD2DEG * atan2((double) item.y2 - item.y1,
                                                (double) item.x2 - item.x1);

        int direction = SPOKESPACING * (floor(0.5 + exact_direction / SPOKESPACING));

        DrawSpokesGL(item.x1, item.y1, item.r, item.width, direction);
        DrawSpokesGL(item.x2, item.y2, item.r2, item.width2, direction + 180);

        break;
      }

      case draw_item::Label:
      {
        if (i==_selected_item) {
          glColor3f(1.0f, 0.0f, 0.0f);
        }

        QFont f=QFont();
        if (item.font.size()) {
          f.fromString(item.font.c_str());
        } else {
          //FIXME: set from prefs?
          f.setPointSizeF(item.font_size);
        }

        myRenderText(item.x1 + item.off * item.off_x,
                   item.y1 + item.off * item.off_y,
                   0.0f, item.text.c_str(),f);
        break;
      }
    }
    glPopAttrib();
  }
}

void GLDrawingCanvas::FitToPage(std::vector<float> origin,
                                std::vector<float> dimensions) {
  orig_x=origin[0];
  orig_y=origin[1];
  dim_x=dimensions[0];
  dim_y=dimensions[1];
  makeCurrent();
  resizeGL(width(), height());
}


void GLDrawingCanvas::Finish() {
  updateGL();
}



static float dsq(float x1, float y1, float x2, float y2) {
  return ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}


void GLDrawingCanvas::mousePressEvent(QMouseEvent *e) {
  e->accept();
  _selected_item=UINT_MAX;

  float x,y;
  toWorldCoords(e->pos(),x,y);

  unsigned int closest=UINT_MAX;
  float d, cx, cy, cd=999.0;
  for (unsigned int i=0; i < _items.size(); ++i) {
    if (_items[i].type==draw_item::Label) {
      if (closest == UINT_MAX) {
        cx=_items[i].x1;
        cy=_items[i].y1;
        cd=dsq(_items[i].x1,_items[i].y1, x, y);
        closest=i;
      }
      d=dsq(_items[i].x1,_items[i].y1, x,y);
      if (d < cd) {
        closest=i;
        cx=_items[i].x1;
        cy=_items[i].y1;
        cd=d;
      }
    }
  }

  if (cd < 0.25) { // approx 1/3 C-C bond length, squared
    _selected_item=closest;
    updateGL();
  }
}

void GLDrawingCanvas::mouseReleaseEvent(QMouseEvent *e) {
  e->accept();
  updateGL();
}

void GLDrawingCanvas::mouseMoveEvent(QMouseEvent *e) {
  e->accept();
  if (_selected_item == UINT_MAX)
    return;
  float x,y;
  toWorldCoords(e->pos(),x,y);
  _items[_selected_item].x1=x;
  _items[_selected_item].y1=y;
  updateGL();
}


void GLDrawingCanvas::OnEditCopy() {
  QApplication::clipboard()->setPixmap(renderPixmap());
}

void GLDrawingCanvas::OnPrint() {
  QPixmap imageToPrint = renderPixmap();   //TODO: could alter canvas size here
  QPrinter prn;
  QPrintDialog printDialog( &prn, this );
  if (printDialog.exec()) {
    QPainter painter(&prn);
    QRect rect = painter.viewport( );
    QSize size = imageToPrint.size( );

    size.scale( rect.size( ), Qt::KeepAspectRatio );
    painter.setViewport( rect.x( ), rect.y( ), size.width( ), size.height( ) );
    painter.setWindow( imageToPrint.rect( ) );
    painter.drawImage( 0, 0, imageToPrint.toImage( ) );
  }
}

void GLDrawingCanvas::OnExportImage() {
  std::string filter = "JPEG Image format (*.jpg)|*.jpg"
                       "|PNG Image format (*.png)|*.png"
                       "|TIFF Image format (*.tif)|*.tif";

  MIFileDialog dlg(this, "Export Image As",
                   QDir::currentPath().toStdString(), "export", filter, MI_SAVE_MODE);
  MIData data;
  data["path"].str = "./export";
  data["filterIndex"].radio = 0;
  data["filterIndex"].radio_count = 5;

  if (!dlg.GetResults(data)) {
    return;
  }

  std::string defaultExt="";
  std::string ext(file_extension(data["path"].str.c_str()));
  switch (data["filterIndex"].radio) {
    default:
    case 0:
      defaultExt = ".jpg";
      break;
    case 1:
      defaultExt = ".png";
      break;
    case 2:
      defaultExt = ".tif";
      break;
  }
  if (ext.size()==0) {
    data["path"].str += defaultExt;
  }

  QPixmap image = renderPixmap();   //TODO: could alter canvas size here
  image.save(data["path"].str.c_str());
}

void GLDrawingCanvas::OnSelectFont() {
  std::string f=_current_font;
  if (MIGetFontFromUser(f)) {
    _current_font=f;
  }
}

void GLDrawingCanvas::OnAddAnnotation() {
  MIGetStringDialog dlg(0,"Get annotation","Enter annotation");
  MIData dat;
  dat["val"].str="Custom annotation";
  if (!dlg.GetResults(dat))
    return;

  float x,y;
  toWorldCoords(QPoint(int(0.3*width()),int(0.2*height())), x, y);
  DrawLabel(x, y, 0.0f, 0.0f, 0.0f, Drawing::ANNOTATION_LABEL_SIZE, dat["val"].str);
  updateGL();
}


void GLDrawingCanvas::OnChangeFonts(const MIActionEvent &evt) {
  std::string f=_current_font;
  if (!MIGetFontFromUser(f)) {
    return;
  }

  if (evt.GetId() == ID_ASP_CURRENT_ITEM_FONT && _selected_item != UINT_MAX) {
    _items[_selected_item].font=f;
    updateGL();
    return;
  }

  float fsize=0.0f;
  switch (evt.GetId()) {
    case ID_ASP_ATOM_FONTS:          fsize=Drawing::ATOM_LABEL_SIZE;  break;
    case ID_ASP_RESIDUE_FONTS:       fsize=Drawing::RESIDUE_LABEL_SIZE;  break;
    case ID_ASP_HYDROGEN_BOND_FONTS: fsize=Drawing::HBOND_LABEL_SIZE;  break;
    case ID_ASP_ANNOTATION_FONTS:    fsize=Drawing::ANNOTATION_LABEL_SIZE; break;
    default: break;
  }

  for (size_t i=0; i<_items.size(); ++i) {
    draw_item &item=_items[i];
    if (item.text.size() &&
        (evt.GetId() == ID_ASP_ALL_FONTS || item.font_size == fsize)) {
      item.font=f;
    }
  }
  updateGL();
}
