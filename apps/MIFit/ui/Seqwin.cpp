#ifdef USE_SEQ_WINDOW

#include <vector>
#include <algorithm>
#include <boost/bind.hpp>

#include "chemlib.h"
#include "corelib.h"

#include "molw.h"
#include "Seqwin.h"
#include "id.h"
#include "MainFrame.h"
#include "wxCMolwView.h"
#include "find_control.h"
#include "MainFrame.h"
#include "MIHistory.h"

using namespace chemlib;

BEGIN_EVENT_TABLE(SequenceWindow, wxPanel)
EVT_MOUSE_EVENTS(SequenceWindow::OnMouseEvent)
EVT_KEY_DOWN(SequenceWindow::OnKeyEvent)
EVT_PAINT(SequenceWindow::OnPaint)
EVT_SIZE(SequenceWindow::OnSize)
EVT_TIMER(ID_SEQWIN_TIMER, SequenceWindow::OnTimer)
//EVT_COMMAND_SCROLL(ID_SEQWINSCROLLBAR, SequenceWindow::OnScroll)
END_EVENT_TABLE()


SequenceWindow::SequenceWindow(wxView* v, wxWindow* parent, const wxPoint& pos, const wxSize& size, long style)
  : wxPanel(parent, -1, pos, size, style), rowheight(0), Resizing(false) {
  Models = NULL;
  MouseStillTime = 0;
  MouseInWindow = false;
  ToolTipInterval = MIConfig::Instance()->GetProfileInt("View Parameters", "ToolTipInterval", 190);
  ClockTick = 100;
  m_timer = new wxTimer(this, ID_SEQWIN_TIMER);
  m_timer->Start(ClockTick);
  Models = NULL;
#ifdef _WIN32
  bigFont = wxFont(8, wxDEFAULT, wxNORMAL, wxLIGHT);
  medFont = wxFont(7, wxDEFAULT, wxNORMAL, wxLIGHT);
  smallFont = wxFont(6, wxDEFAULT, wxNORMAL, wxLIGHT);
#else
  bigFont = wxFont(10, wxDEFAULT, wxNORMAL, wxLIGHT);
  medFont = wxFont(8, wxDEFAULT, wxNORMAL, wxLIGHT);
  smallFont = wxFont(6, wxDEFAULT, wxNORMAL, wxLIGHT);
#endif
  theView = (CMolwView*)v;
  Models = theView->GetDocument()->GetDisplaylist();
  currentChangedConnection=Models->currentMoleculeChanged.connect(boost::bind(&SequenceWindow::currentModelChanged, this, _1, _2));
  AtomStack = theView->AtomStack;
  m_palette = Application::instance()->GetLPpal();
  pDC = NULL;
  popup_menu = new wxMenu(wxMENU_TEAROFF);
  popup_menu->Append(ID_COLOR_RESIDUE, "Color residue");
  popup_menu->Append(ID_CPK_RESIDUE, "CPK Radius residue");
  popup_menu->AppendSeparator();
  popup_menu->Append(ID_HIDE_RESIDUE, "Show/Hide residue");
  maxy = 200;
  scroll = 0;
}

void SequenceWindow::currentModelChanged(Molecule*, Molecule*) {
  currentChangedConnection.disconnect();
  Refresh();
}

void SequenceWindow::OnPaint(wxPaintEvent& WXUNUSED (event)) {
  if (!theView->IsCreated()) {
    return;
  }
  wxPaintDC dc(this);
  PrepareDC(dc);
  OnDraw(&dc);
}

SequenceWindow::~SequenceWindow() {
  delete popup_menu;
  delete m_timer;
}

void SequenceWindow::OnKeyEvent(wxKeyEvent& event) {
  if (!theView) {

    event.Skip();
  } else {
    theView->handleKeyEvent(event);
    if (theView->ViewChanged()) {
      Refresh(false);
    }
  }
}

void SequenceWindow::OnDraw(wxDC* theDC) {
  if (!theView->IsCreated()) {
    return;
  }
  pDC = theDC;
  wxScrollBar* scrollbar;  findControl(scrollbar, this, ID_SEQWINSCROLLBAR);
  scroll = scrollbar->GetThumbPosition()*rowheight;
  DrawSelf();
  pDC = NULL;
  scrollbar->SetScrollbar(scroll/rowheight, (frame.bottom-frame.top)/rowheight, maxy/rowheight+1, rowheight);
}

void SequenceWindow::OnMouseEvent(wxMouseEvent& event) {
  if (!theView->IsCreated()) {
    return;
  }
  unsigned short flags = 0;
  if (!theView) {
    return;
  }

  wxClientDC dc(this);
  PrepareDC(dc);
  pDC = &dc;

  dc.SetPen(*wxBLACK_PEN);

  wxPoint pt(event.GetLogicalPosition(dc));

  // flags for the key states
  if (event.ControlDown()) {
    flags += MK_CONTROL;
  }
  if (event.ShiftDown()) {
    flags += MK_SHIFT;
  }
  if (event.LeftIsDown()) {
    flags += MK_LBUTTON;
  }
  if (event.RightIsDown()) {
    flags += MK_RBUTTON;
  }

  //handle the type of event
  if (event.Entering()) {
    MouseInWindow = true;
  } else if (event.Leaving()) {
    MouseInWindow = false;
  } else {
    if (event.Dragging() || event.Moving()) {
      OnMouseMove(flags, CPoint((int)pt.x, (int)pt.y));
    }
    if (event.LeftUp()) {
      OnLButtonUp(flags, CPoint((int)pt.x, (int)pt.y));
    }
    if (event.LeftDown()) {
      OnLButtonDown(flags, CPoint((int)pt.x, (int)pt.y));
    }
    if (event.LeftDClick()) {
      OnLButtonDblClk(flags, CPoint((int)pt.x, (int)pt.y));
    }
    if (event.RightUp()) {
      OnRButtonUp(flags, CPoint((int)pt.x, (int)pt.y));
    }
    if (event.RightDown()) {
      OnRButtonDown(flags, CPoint((int)pt.x, (int)pt.y));
    }
  }

  if (PaletteChanged) {
    DrawSelf();
    PaletteChanged = false;
  }
  pDC = NULL;
}

bool SequenceWindow::CharPosition(int seqpos, int modelno, int& x, int& y) {
  x = frame.left + start + cwidth*(seqpos - (seqpos/nrow)*nrow);
  x += (seqpos - (seqpos/nrow)*nrow)/10*cwidth;
  y = top + rowheight*(seqpos/nrow);
  y += modelno * 12;
  // keep track of maximum y for scrollbar purposes
  if (y > maxy) {
    maxy = y;
  }
  //y += frame.top;
  if (y < frame.top /*+scroll*/ || y > frame.bottom+rowheight /*+scroll*/) {
    return false;
  }
  y -= scroll;
  return true;
}

void TextDraw(int x, int y, const char* s);

//@{
// Imitation of MFC CSize for portability purposes.
//@}
struct CSize
{
  int cx;
  int cy;
};

void SequenceWindow::DrawSelf() {
  if (!theView->IsCreated()) {
    return;
  }
  if (pDC == NULL) {
    return;
  }
  wxString buf;
  char id[50];
  RESIDUE* prev_res = NULL;
  int x, y;
  if (!theView) {
    return;
  }
  Models = theView->GetDocument()->GetDisplaylist();
  AtomStack = theView->AtomStack;
  if (AtomStack == NULL) {
    return;
  }
  m_palette = Application::instance()->GetLPpal();
  char c;
  CSize charsize;
  CSize charsize2;
  wxCoord descent;
  wxCoord descent2;
  //pDC->SelectObject(&bigFont);
  pDC->SetFont(bigFont);
  pDC->SetPen(*wxBLACK_PEN);
  pDC->SetTextForeground(*wxBLACK);
  pDC->SetBackground(*wxWHITE_BRUSH);
  pDC->SetBackgroundMode(wxSOLID);
  pDC->Clear();

  frame.left = 0;
  frame.top = scroll;
  GetClientSize(&x, &y);
  frame.right = x;
  frame.bottom = y+scroll;

  pDC->GetTextExtent("W", &charsize.cx, &charsize.cy, &descent);
  pDC->GetTextExtent("Q", &charsize2.cx, &charsize2.cy, &descent2);
  if (charsize2.cx > charsize.cx) {
    charsize.cx = charsize2.cx;
  }
  if (charsize2.cy > charsize.cy) {
    charsize.cy = charsize2.cy;
  }
  if (descent2 > descent) {
    descent = descent2;
  }
  cwidth = charsize.cx-1;
  bigY = (charsize.cy + descent)*8/10;
  //pDC->SelectObject(&medFont);
  pDC->SetFont(medFont);
  pDC->GetTextExtent("W", &charsize.cx, &charsize.cy, &descent);
  pDC->GetTextExtent("Q", &charsize2.cx, &charsize2.cy, &descent2);
  if (charsize2.cx > charsize.cx) {
    charsize.cx = charsize2.cx;
  }
  if (charsize2.cy > charsize.cy) {
    charsize.cy = charsize2.cy;
  }
  if (descent2 > descent) {
    descent = descent2;
  }
  medY = (charsize.cy+descent)*8/10;
  //pDC->SelectObject(&smallFont);
  pDC->SetFont(smallFont);
  pDC->GetTextExtent("W", &charsize.cx, &charsize.cy, &descent);
  pDC->GetTextExtent("Q", &charsize2.cx, &charsize2.cy, &descent2);
  if (charsize2.cx > charsize.cx) {
    charsize.cx = charsize2.cx;
  }
  if (charsize2.cy > charsize.cy) {
    charsize.cy = charsize2.cy;
  }
  if (descent2 > descent) {
    descent = descent2;
  }
  smallY = (charsize.cy+descent)*8/10;
  pDC->SetFont(bigFont);
  //start = 40;
  start = cwidth/2;
  top = 54;
  rowheight = 52;
  maxy = 0;
  nrow = (((frame.right-frame.left)-start -15 /* scrolbar width*/)/cwidth)/11*10;
  if (nrow < 1) {
    nrow = 1;
  }
  if (!Models) {
    return;
  }
  int modelno = 0;
  nmodels = 1;
  Molecule* node = Models->CurrentItem();
  if (!node) {
    return;
  }
  int color = Colors::sec_colors[Colors::HELIX];
  color = PaletteIndex(color)+4;
  HelixRGB.SetColour(wxColour(m_palette->colors[color].red,
      m_palette->colors[color].green,
      m_palette->colors[color].blue));
  color = Colors::sec_colors[Colors::SHEET];
  color = PaletteIndex(color)+4;
  SheetRGB.SetColour(wxColour(m_palette->colors[color].red,
      m_palette->colors[color].green,
      m_palette->colors[color].blue));
  color = Colors::sec_colors[Colors::COIL];
  color = PaletteIndex(color)+4;
  CoilRGB.SetColour(wxColour(m_palette->colors[color].red,
      m_palette->colors[color].green,
      m_palette->colors[color].blue));
  color = Colors::sec_colors[Colors::TURN];
  color = PaletteIndex(color)+4;
  TurnRGB.SetColour(wxColour(m_palette->colors[color].red,
      m_palette->colors[color].green,
      m_palette->colors[color].blue));

  pDC->SetFont(bigFont);
  buf = "Model: ";
  buf += node->compound.c_str();
  wxCoord tw, th;
  wxCoord fw = (wxCoord)(frame.right-frame.left-25);
  pDC->GetTextExtent(buf, &tw, &th);
  //int ncpy = (frame.right-frame.left-25)/cwidth;
  if (fw < tw) {
    int ncpy = buf.Length()*fw/tw;
    wxString trim = buf.Right(ncpy);
    buf = trim;
  }
  //pDC->TextOut(frame.left+4, frame.top + top -scroll-bigY, buf, strlen(buf));
  pDC->DrawText(buf, frame.left+3, 1-scroll);
  pDC->SetBackgroundMode(wxSOLID);
  char alt_name1;
  int in_stack;

  for (MIIter<RESIDUE> res = node->GetResidues(); res; ++res) {
    if (CharPosition(res->seqpos, modelno, x, y)) {
      if ((in_stack = AtomStack->InStack(res)) != 0) {
        if (in_stack == 1) {
          pDC->SetTextBackground(wxColour(0, 255, 255));
        } else {
          pDC->SetTextBackground(wxColour(180, 255, 255));
        }
      } else {
        pDC->SetTextBackground(wxColour(255, 255, 255));
      }
      buf = res->name1();
      pDC->DrawText(buf, x, y-bigY);
      pDC->SetTextBackground(wxColour(255, 255, 255));
      //alternate sequence
      alt_name1 = node->GetSeq(res->seqpos);
      buf = alt_name1;
      pDC->DrawText(buf, x, y+1);
      DrawColor(res, x, y);
    }
    prev_res = res;
  }
  pDC->SetTextBackground(wxColour(255, 255, 255));
  for (int i = 0; i < node->SeqMax(); i++) {
    if (CharPosition(i, modelno, x, y)) {
      for (MIIter<RESIDUE> res = node->GetResidues(); res; ++res) {
        if (res->seqpos == (unsigned int)i) {
          goto found;
        }
      }
      pDC->DrawText(_T("."), x, y-bigY);
      //alternate sequence
      buf = alt_name1 = node->GetSeq(i);
      pDC->DrawText(buf, x, y+1);
    }
found: ;
  }
  pDC->SetBackgroundMode(wxTRANSPARENT);

  for (MIIter<RESIDUE> res = node->GetResidues(); res; ++res) {
    if (CharPosition(res->seqpos, modelno, x, y)) {
      if (IsPeptide(*res)) {
        DrawSecStr(res, x, y-20, prev_res);
        if (chargetype(res->name1()) != ' ') {
          buf = c = chargetype(res->name1());
          pDC->SetFont(smallFont);
          pDC->DrawText(buf, x, y-30-smallY);
        }
      }
      if ((res->seqpos)%10 == 0) {
        sprintf(id, "%s %c", res->name().c_str(), (char)(res->chain_id()&255));
        pDC->SetFont(medFont);
        buf = id;
        pDC->DrawText(buf, x, y-10-medY);
      }
    }
    prev_res = res;
  }

  pDC->SetFont(wxNullFont);
  pDC->SetPen(wxNullPen);
  pDC->SetBrush(wxNullBrush);
}

void SequenceWindow::DrawColor(RESIDUE* res, int x, int y) {
  MIAtom* a;
  int color;
  if ((a = atom_from_name("CA", *res)) == NULL) {
    a = res->atom(0);
  }
  if (a) {
    color = abs(a->color);
    color = PaletteIndex(color)+2;
    wxColour col(m_palette->colors[color].red,
                 m_palette->colors[color].green,
                 m_palette->colors[color].blue);
    wxPen pen(col, 1, wxSOLID);
    pDC->SetPen(pen);
    pDC->DrawLine(x, y+1, x+cwidth, y+1);
    pDC->DrawLine(x, y+2, x+cwidth, y+2);
    pDC->SetPen(wxNullPen);
  }
}

//#define FORECOLOR(c) pDC->SetColour(c);
//#define MoveTo pDC->MoveTo
//#define LineTo pDC->LineTo

void SequenceWindow::DrawSecStr(RESIDUE* res, int x, int y, RESIDUE* prev_res) {
  if (!pDC) {
    return;
  }
  if (res->secstr == 'C' || res->secstr == 'U') {
    pDC->SetPen(CoilRGB);
    pDC->DrawLine(x, y-4, x+cwidth, y-4);
  } else if (res->secstr == 'T') {
    pDC->SetPen(TurnRGB);
    pDC->DrawLine(x, y-4, x+cwidth/2, y-7);
    pDC->DrawLine(x+cwidth/2, y-7, x+cwidth, y-4);
  } else if (res->secstr == 'H') {
    pDC->SetPen(HelixRGB);
    pDC->SetBrush(*wxTRANSPARENT_BRUSH);
    pDC->DrawArc(x+cwidth/2, y-4, x, y-3, x+cwidth/4, y-4);
    pDC->DrawArc(x+cwidth/2, y-3, x+cwidth, y-4, x+3*cwidth/4, y-4);
  } else if (res->secstr == 'S') {
    pDC->SetPen(SheetRGB);
    if (res->next() == NULL || res->next()->secstr != 'S') {
      pDC->DrawLine(x, y-8, x+cwidth, y-4);
      pDC->DrawLine(x+cwidth, y-4, x, y);
      pDC->DrawLine(x, y, x, y-8);
    } else {
      pDC->DrawLine(x, y-2, x+cwidth, y-2);
      pDC->DrawLine(x, y-6, x+cwidth, y-6);
    }
    if (prev_res == NULL || prev_res->secstr != 'S') {
      pDC->DrawLine(x, y-2, x, y-6);
    }
  }
  pDC->SetPen(*wxBLACK_PEN);
}

void SequenceWindow::OnMouseMove(unsigned short WXUNUSED (nFlags), CPoint point) {
  if (!theView->IsCreated()) {
    return;
  }
  if (Resizing) {
    if (lastdraw.y != point.y) {
      // don't let draw go below pane or above top
      if (point.y < frame.bottom-5 && point.y > 1) {
      }
    }
  }
  MouseStillTime = 0;
  mouse = point;
}

void SequenceWindow::OnLButtonUp(unsigned short WXUNUSED (nFlags), CPoint WXUNUSED (point)) {
  // this is taken care of by the framework with MFC
}

void SequenceWindow::OnRButtonUp(unsigned short WXUNUSED (nFlags), CPoint point) {
  if (point.x > mousestart.x+1 || point.x < mousestart.x-1
      || point.y > mousestart.y+1 || point.y < mousestart.y-1) {
    return;
  }
  int x, y;
  //MIAtom *a;
  if (!Models) {
    return;
  }
  int modelno = 0;
  nmodels = 1;
  Molecule* node = Models->CurrentItem();
  if (!node) {
    return;
  }
  for (MIIter<RESIDUE> res = node->GetResidues(); res; ++res) {
    if (CharPosition(res->seqpos, modelno, x, y)) {
      y -= scroll;
      if (point.x >= x && point.x < x+cwidth) {
        if (point.y <= y && point.y > y-10) {
          // todo - whatever we want to do when a r-pick occurs...
          break;
        }
      }

    }
  }
}

static bool SendToHistory(const std::string& type, Displaylist* dl,
                          Molecule* model, RESIDUE* res, MIAtom* atom) {
  MIHistory* h = MIGetHistory();
  unsigned int anum = 0, rnum = 0, mnum = 0;

  if (!h->PickSerialize(dl, atom, res, model, anum, rnum, mnum)) {
    return false;
  }
  MIData data;
  data["command"].str = "seqw";
  data["type"].str = type;
  data["atom"].u = anum;
  data["res"].u = rnum;
  data["mol"].u = mnum;
  data["ss"].str = res->secstr;
  return h->AddCommand(data);
}

void SequenceWindow::OnLButtonDown(unsigned short WXUNUSED (nFlags), CPoint point) {
  int x, y;
  mousestart = point;
  if (!Models) {
    return;
  }
  int modelno = 0;
  nmodels = 1;
  Molecule* node = Models->CurrentItem();
  if (!node) {
    return;
  }
  ResidueClicked = NULL;
  GetClientSize(&x, &y);

  for (MIIter<RESIDUE> res = node->GetResidues(); res; ++res) {
    if (CharPosition(res->seqpos, modelno, x, y)) {
      //y-=scroll;
      if (point.x >= x && point.x < x+cwidth) {
        if (point.y <= y && point.y > y-10) {
          pDC->SetPen(*wxBLACK_PEN);
          pDC->DrawLine(x, y, x+cwidth, y);
          pDC->SetPen(wxNullPen);
          theView->select(NULL, res, NULL);
          SendToHistory("select", theView->GetDocument()->GetDisplaylist(), 0, res, 0);
          PaletteChanged = true;
          ResidueClicked = res;
          theView->ReDraw();
          break;
        } else if (point.y < y-20 && point.y > y-29) {
          // rotate through the secondary structure types
          if (res->secstr == 'C') {
            res->secstr = 'H';
          } else if (res->secstr == 'H') {
            res->secstr = 'S';
          } else if (res->secstr == 'S') {
            res->secstr = 'T';
          } else if (res->secstr == 'T') {
            res->secstr = 'C';
          }
          SendToHistory("change", theView->GetDocument()->GetDisplaylist(), 0, res, 0);
          PaletteChanged = true;
          theView->ReDraw();
          break;
        }
      }
    }
  }
  theView->canvas->getWindow()->SetFocus();
}

void SequenceWindow::OnRButtonDown(unsigned short WXUNUSED (nFlags), CPoint WXUNUSED (point)) {
  //PopupMenu(popup_menu, point.x, point.y);
}

void SequenceWindow::OnLButtonDblClk(unsigned short WXUNUSED (nFlags), CPoint WXUNUSED (point)) {
  if (ResidueClicked != NULL) {
    theView->setFocusResidue(ResidueClicked);
    SendToHistory("focus", theView->GetDocument()->GetDisplaylist(), 0, ResidueClicked, 0);
  }
}

void SequenceWindow::OnKeyDown(unsigned short WXUNUSED (nChar), unsigned short WXUNUSED (nRepCnt), unsigned short WXUNUSED (nFlags)) {
}

void SequenceWindow::OnSize(wxSizeEvent&) {
  if (!theView->IsCreated()) {
    return;
  }
  int x, y;
  GetClientSize(&x, &y);
  SeqScroller* scrollbar;  findControl(scrollbar, this, ID_SEQWINSCROLLBAR);
  scrollbar->SetSize(x-18, 0, 18, y);
  Refresh();
}

BEGIN_EVENT_TABLE(SeqScroller, wxScrollBar)
EVT_KEY_DOWN(SeqScroller::OnKeyEvent)
EVT_COMMAND_SCROLL(ID_SEQWINSCROLLBAR, SeqScroller::OnScroll)
END_EVENT_TABLE()

void SeqScroller::OnScroll(wxScrollEvent& event) {
  if (event.GetEventType() == wxEVT_SCROLL_CHANGED) {
    return;
  }
  seqwin->Refresh();
  seqwin->theView->canvas->getWindow()->SetFocus();
  event.Skip();
}

void SeqScroller::OnKeyEvent(wxKeyEvent& event) {
  seqwin->OnKeyEvent(event);
}

void SequenceWindow::OnTimer(wxTimerEvent&) {
  if (!theView->IsCreated()) {
    return;
  }
  if (MouseInWindow) {
    MouseStillTime += ClockTick;
    if (MouseStillTime > ToolTipInterval) {
      OnToolTip();
    }
  }
}

void SequenceWindow::OnToolTip() {
  if (!theView->IsCreated()) {
    return;
  }
  int x, y;
  static int oldx = -1, oldy = -1;
  wxString s;
  if (mouse.x != oldy || mouse.y != oldy) {
    if (!Models) {
      return;
    }
    int modelno = 0;
    nmodels = 1;
    Molecule* node = Models->CurrentItem();
    if (!node) {
      return;
    }
    if (mouse.y + scroll < 20) {
      // in title - display full title on bar
      s = node->compound.c_str();
      oldx = mouse.x;
      oldy = mouse.y;
    } else {

      for (MIIter<RESIDUE> res = node->GetResidues(); res; ++res) {
        if (CharPosition(res->seqpos, modelno, x, y)) {
          if (mouse.x >= x && mouse.x < x + cwidth
              && mouse.y <= y && mouse.y > y-10) {

            oldx = mouse.x;
            oldy = mouse.y;
            s = resid(res).c_str();
            break;
          }
        }
      }
    }
    MIMainWindowRightFooter(s);
  }
  MouseStillTime = 0;
}

void SequenceWindow::ReDraw() {
  if (theView) {
    if (!theView->IsCreated()) {
      return;
    }
    if (theView->IsDrawing()) {
      return;
    }
  }
  wxClientDC dc(this);
  OnDraw(&dc);
}

#endif //USE_SEQ_WINDOW
