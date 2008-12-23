#ifndef mifit_ui_RenderStyle_h
#define mifit_ui_RenderStyle_h

class RenderStyle {

  /**
   * Whether to show atoms as balls.
   */
  bool atomBall;

  bool bondLine;
  float bondLineWidth;
  bool bondCylinder;
  float ballPercent;
  float stickPercent;

  static bool createDefaults();

  // These default values are unaffected by user preference settings.
  static const RenderStyle* defaultLine;
  static const RenderStyle* defaultBallAndLine;
  static const RenderStyle* defaultStick;
  static const RenderStyle* defaultBallAndStick;
  static const RenderStyle* defaultBall;

  static const RenderStyle* line;
  static const RenderStyle* ballAndLine;
  static const RenderStyle* stick;
  static const RenderStyle* ballAndStick;
  static const RenderStyle* ball;

public:

  // These default values are unaffected by user preference settings.
  static const RenderStyle& getDefaultLine();
  static const RenderStyle& getDefaultBallAndLine();
  static const RenderStyle& getDefaultStick();
  static const RenderStyle& getDefaultBallAndStick();
  static const RenderStyle& getDefaultBall();

  static const RenderStyle& getLine();
  static const RenderStyle& getBallAndLine();
  static const RenderStyle& getStick();
  static const RenderStyle& getBallAndStick();
  static const RenderStyle& getBall();

  RenderStyle();

  void set(const RenderStyle& style);

  bool isAtomBall();
  void setAtomBall(bool on);

  float getBallPercent();
  void setBallPercent(float percent);

  bool isBondLine();
  void setBondLine(bool on);

  float getBondLineWidth();
  void setBondLineWidth(float width);

  bool isBondCylinder();
  void setBondCylinder(bool on);

  float getStickPercent();
  void setStickPercent(float percent);

};

#endif
