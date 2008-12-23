#ifndef SequenceWIND
#define SequenceWIND

#ifdef USE_SEQ_WINDOW


class MIPalette;
class CMolwView;
class Displaylist;
class Stack;

class SequenceWindow : public wxPanel {
public:
  CMolwView* theView;
  int start, top;
  int rowheight;
  int cwidth;
  int nmodels;
  int nrow;
  int maxy, scroll;
  short font;
  RECT frame;
  wxPen CoilRGB, HelixRGB, SheetRGB, TurnRGB;
  wxFont bigFont, smallFont, medFont;
  int bigY, medY, smallY;
private:
  void DrawSecStr(chemlib::RESIDUE *res, int, int, chemlib::RESIDUE *prev_res);
  void DrawColor(chemlib::RESIDUE* res, int x, int y);
  chemlib::RESIDUE* ResidueClicked;
  CPoint mousestart;
  bool Resizing;
  CPoint lastdraw;
  long MouseStillTime;
  long ToolTipInterval;
  bool MouseInWindow;
  int ClockTick;
  wxTimer* m_timer;
  wxMenu* popup_menu;   // popup menu for right click
  //LGWorld * mGWorld;
  boost::signals::connection currentChangedConnection;


public:
  void ReDraw();
  CPoint mouse;
  void OnToolTip();
  void OnTimer(wxTimerEvent& evt);
  void OnSize(wxSizeEvent& evt);
  Displaylist* Models;
  Stack* AtomStack;
  MIPalette* m_palette;
  wxDC* pDC;
  SequenceWindow(wxView* v, wxWindow* parent, const wxPoint& pos, const wxSize& size, long style);
  ~SequenceWindow();
  bool CharPosition(int seqpos, int modelno, int& x, int& y);
  virtual void OnDraw(wxDC*);
  void OnPaint(wxPaintEvent& event);
  void DrawSelf();
  void OnLButtonUp(unsigned short nFlags, CPoint point);
  void OnRButtonUp(unsigned short nFlags, CPoint point);
  void OnLButtonDown(unsigned short nFlags, CPoint point);
  void OnRButtonDown(unsigned short nFlags, CPoint point);
  void OnLButtonDblClk(unsigned short nFlags, CPoint point);
  void OnKeyDown(unsigned short nChar, unsigned short nRepCnt, unsigned short nFlags);
  void OnMouseMove(unsigned short nFlags, CPoint point);
  void OnMouseEvent(wxMouseEvent& event);
  void OnKeyEvent(wxKeyEvent& event);
private:
  void currentModelChanged(Molecule* oldModel, Molecule* newModel);

  DECLARE_EVENT_TABLE();
};

class SeqScroller : public wxScrollBar {
private:
  SequenceWindow* seqwin;
public:
  SeqScroller(SequenceWindow* parent, int id, const wxPoint& pos, const wxSize& size, long style)
    : wxScrollBar(parent, id, pos, size, style) {
    seqwin = parent;
  }

  void OnScroll(wxScrollEvent& event);
  void OnKeyEvent(wxKeyEvent& event);
private:
  DECLARE_EVENT_TABLE();
};

#endif //USE_SEQ_WINDOW

#endif
