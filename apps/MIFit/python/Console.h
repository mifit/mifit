#ifndef CONSOLE_h
#define CONSOLE_h

#include <QTextEdit>

/**
 * This Console class is a text editing widget which mimics the behavior of a
 * command line window. The prompt slot sets the prefix to the command line.
 * A return or enter key will cause a command signal with the text after the
 * prompt. Text before the prompt can't be edited, but can be selected and
 * copied using mouse selection. Keyboard cursor navigation and selection is
 * restricted to the command line.
 *
 * The message and prompt slots are for displaying text in the widget. The
 * message slot simply appends the text. The prompt slot appends the text and
 * sets a new start of command position. The start of command position is what
 * restricts the editing capabilities to the command line.
 * 
 * Note that a prompt is not automatically displayed after each command is
 * signaled. A call to the prompt slot must be performed for each prompting
 * of a new command.
 *
 * Although this class is generic, its design was primarily for operation with
 * an embedded Python interactive loop (the PyRun_InteractiveLoop function).
 *
 * The default font is "Lucida Console" with a style hint for TypeWriter. To
 * change the font, see the QTextEdit documentation.
 */
class Console : public QTextEdit {

  Q_OBJECT

  /**
   * Location marking the beginning of the command line.
   */
  int commandStart;

protected:

  /**
   * Filters and processes key events for restricted editing behavior.
   */
  virtual void keyPressEvent(QKeyEvent* event);

Q_SIGNALS:

  /**
   * Signals the text for a command when return or enter pressed.
   */
  void command(const QString& text);

public:

  /**
   * Creates an empty Console with given parent.
   * See QTextEdit for details.
   */
  Console(QWidget* parent = 0);
  
  /**
   * Creates an Console with given parent and text.
   * See QTextEdit for details.
   */
  Console(const QString& text, QWidget* parent = 0);

public Q_SLOTS:

  /**
   * Appends the prompt and sets a new start of command position.
   */
  void prompt(const QString& prompt);

  /**
   * Appends the text.
   */
  void message(const QString& text);
  
};

#endif
