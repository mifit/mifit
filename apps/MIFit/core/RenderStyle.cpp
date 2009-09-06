#include "core/corelib.h"

#include "RenderStyle.h"

const RenderStyle* RenderStyle::defaultLine = NULL;
const RenderStyle* RenderStyle::defaultBallAndLine = NULL;
const RenderStyle* RenderStyle::defaultStick = NULL;
const RenderStyle* RenderStyle::defaultBallAndStick = NULL;
const RenderStyle* RenderStyle::defaultBall = NULL;

const RenderStyle* RenderStyle::line = NULL;
const RenderStyle* RenderStyle::ballAndLine = NULL;
const RenderStyle* RenderStyle::stick = NULL;
const RenderStyle* RenderStyle::ballAndStick = NULL;
const RenderStyle* RenderStyle::ball = NULL;

bool RenderStyle::createDefaults() {
  float defaultBallSize = 10.0f / 100.0f;
  float defaultStickPercent = defaultBallSize * 33.0f / 100.0f;
  float ballPercent = (float) MIConfig::Instance()->GetProfileInt("View Parameters", "ballsize", 10) / 100.0f;
  float stickPercent = ballPercent
                       * (float) MIConfig::Instance()->GetProfileInt("View Parameters", "cylindersize", 33) / 100.0f;

  RenderStyle* style = new RenderStyle;
  style->setAtomBall(false);
  style->setBallPercent(defaultBallSize);
  style->setBondLine(true);
  style->setBondCylinder(false);
  style->setStickPercent(defaultStickPercent);
  delete defaultLine;
  defaultLine = style;

  style = new RenderStyle;
  style->setAtomBall(true);
  style->setBallPercent(defaultBallSize);
  style->setBondLine(true);
  style->setBondCylinder(false);
  style->setStickPercent(defaultStickPercent);
  delete defaultBallAndLine;
  defaultBallAndLine = style;

  style = new RenderStyle;
  style->setAtomBall(false);
  style->setBallPercent(defaultBallSize);
  style->setBondLine(false);
  style->setBondCylinder(true);
  style->setStickPercent(defaultStickPercent);
  delete defaultStick;
  defaultStick = style;

  style = new RenderStyle;
  style->setAtomBall(true);
  style->setBallPercent(defaultBallSize);
  style->setBondLine(false);
  style->setBondCylinder(true);
  style->setStickPercent(defaultStickPercent);
  delete defaultBallAndStick;
  defaultBallAndStick = style;

  style = new RenderStyle;
  style->setAtomBall(true);
  style->setBallPercent(1.0f);
  style->setBondLine(false);
  style->setBondCylinder(false);
  style->setStickPercent(defaultStickPercent);
  delete defaultBall;
  defaultBall = style;

  style = new RenderStyle;
  style->setAtomBall(false);
  style->setBallPercent(ballPercent);
  style->setBondLine(true);
  style->setBondCylinder(false);
  style->setStickPercent(stickPercent);
  delete line;
  line = style;

  style = new RenderStyle;
  style->setAtomBall(true);
  style->setBallPercent(ballPercent);
  style->setBondLine(true);
  style->setBondCylinder(false);
  style->setStickPercent(stickPercent);
  delete ballAndLine;
  ballAndLine = style;

  style = new RenderStyle;
  style->setAtomBall(false);
  style->setBallPercent(ballPercent);
  style->setBondLine(false);
  style->setBondCylinder(true);
  style->setStickPercent(stickPercent);
  delete stick;
  stick = style;

  style = new RenderStyle;
  style->setAtomBall(true);
  style->setBallPercent(ballPercent);
  style->setBondLine(false);
  style->setBondCylinder(true);
  style->setStickPercent(stickPercent);
  delete ballAndStick;
  ballAndStick = style;

  style = new RenderStyle;
  style->setAtomBall(true);
  style->setBallPercent(1.0f);
  style->setBondLine(false);
  style->setBondCylinder(false);
  style->setStickPercent(stickPercent);
  delete ball;
  ball = style;

  return true;
}

const RenderStyle& RenderStyle::getDefaultLine() {
  if (line == NULL) {
    createDefaults();
  }
  return *defaultLine;
}

const RenderStyle& RenderStyle::getDefaultBallAndLine() {
  if (ballAndLine == NULL) {
    createDefaults();
  }
  return *defaultBallAndLine;
}

const RenderStyle& RenderStyle::getDefaultStick() {
  if (stick == NULL) {
    createDefaults();
  }
  return *defaultStick;
}

const RenderStyle& RenderStyle::getDefaultBallAndStick() {
  if (ballAndStick == NULL) {
    createDefaults();
  }
  return *defaultBallAndStick;
}

const RenderStyle& RenderStyle::getDefaultBall() {
  if (ball == NULL) {
    createDefaults();
  }
  return *defaultBall;
}

const RenderStyle& RenderStyle::getLine() {
  if (line == NULL) {
    createDefaults();
  }
  return *line;
}

const RenderStyle& RenderStyle::getBallAndLine() {
  if (ballAndLine == NULL) {
    createDefaults();
  }
  return *ballAndLine;
}

const RenderStyle& RenderStyle::getStick() {
  if (stick == NULL) {
    createDefaults();
  }
  return *stick;
}

const RenderStyle& RenderStyle::getBallAndStick() {
  if (ballAndStick == NULL) {
    createDefaults();
  }
  return *ballAndStick;
}

const RenderStyle& RenderStyle::getBall() {
  if (ball == NULL) {
    createDefaults();
  }
  return *ball;
}

RenderStyle::RenderStyle() : bondLineWidth(2.0f) {
}

void RenderStyle::set(const RenderStyle& style) {
  atomBall = style.atomBall;
  bondLine = style.bondLine;
  bondLineWidth = style.bondLineWidth;
  bondCylinder = style.bondCylinder;
  ballPercent = style.ballPercent;
  stickPercent = style.stickPercent;
}

bool RenderStyle::isAtomBall() {
  return atomBall;
}

void RenderStyle::setAtomBall(bool on) {
  atomBall = on;
}

float RenderStyle::getBallPercent() {
  return ballPercent;
}

void RenderStyle::setBallPercent(float percent) {
  ballPercent = percent;
}

bool RenderStyle::isBondLine() {
  return bondLine;
}

void RenderStyle::setBondLine(bool on) {
  bondLine = on;
}

float RenderStyle::getBondLineWidth() {
  return bondLineWidth;
}

void RenderStyle::setBondLineWidth(float width) {
  bondLineWidth = width;
}

bool RenderStyle::isBondCylinder() {
  return bondCylinder;
}

void RenderStyle::setBondCylinder(bool on) {
  bondCylinder = on;
}

float RenderStyle::getStickPercent() {
  return stickPercent;
}

void RenderStyle::setStickPercent(float percent) {
  stickPercent = percent;
}

