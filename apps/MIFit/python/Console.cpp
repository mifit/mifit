#include "Console.h"

#include <QKeyEvent>
#include <QDebug>

Console::Console(QWidget *parent)
    : QTextEdit(parent),
      commandStart(0)
{

    QFont f("Lucida Console");
    f.setStyleHint(QFont::TypeWriter);
    setFont(f);
}

Console::Console(const QString &text, QWidget *parent)
    : QTextEdit(text, parent),
      commandStart(0)
{

    QFont f("Lucida Console");
    f.setStyleHint(QFont::TypeWriter);
    setFont(f);

    // Set the command start to after any text provided.
    QTextCursor cursor = textCursor();
    cursor.movePosition(QTextCursor::End);
    commandStart = cursor.position();
}

void Console::keyPressEvent(QKeyEvent *event)
{
    QTextCursor cursor = textCursor();

    // Always ignore these key sequences.
    if (event->matches(QKeySequence::MoveToNextLine)
        || event->matches(QKeySequence::MoveToPreviousLine)
        || event->matches(QKeySequence::MoveToNextPage)
        || event->matches(QKeySequence::MoveToPreviousPage)
        || event->matches(QKeySequence::SelectNextLine)
        || event->matches(QKeySequence::SelectPreviousLine)
        || event->matches(QKeySequence::SelectNextPage)
        || event->matches(QKeySequence::SelectPreviousPage))
    {

        return;
    }

    if (cursor.position() == commandStart)
    {

        // Ignore these key sequences when cursor at start of command line.
        if (event->matches(QKeySequence::Back)
            || event->matches(QKeySequence::MoveToPreviousChar)
            || event->matches(QKeySequence::MoveToPreviousWord)
            || event->matches(QKeySequence::MoveToStartOfLine)
            || event->matches(QKeySequence::MoveToStartOfBlock)
            || event->matches(QKeySequence::MoveToEndOfDocument)
            || event->matches(QKeySequence::SelectPreviousChar)
            || event->matches(QKeySequence::SelectPreviousWord)
            || event->matches(QKeySequence::SelectPreviousLine)
            || event->matches(QKeySequence::SelectPreviousPage)
            || event->matches(QKeySequence::SelectStartOfLine)
            || event->matches(QKeySequence::SelectEndOfLine)
            || event->matches(QKeySequence::SelectStartOfBlock)
            || event->matches(QKeySequence::SelectEndOfBlock)
            || event->matches(QKeySequence::SelectStartOfDocument)
            || event->matches(QKeySequence::SelectEndOfDocument))
        {

            return;
        }
    }
    else if (cursor.position() < commandStart)
    {

        // Ignore these key sequences when cursor before command line.
        if (event->matches(QKeySequence::Back)
            || event->matches(QKeySequence::Forward)
            || event->matches(QKeySequence::MoveToNextChar)
            || event->matches(QKeySequence::MoveToPreviousChar)
            || event->matches(QKeySequence::MoveToNextWord)
            || event->matches(QKeySequence::MoveToPreviousWord)
            || event->matches(QKeySequence::MoveToStartOfLine)
            || event->matches(QKeySequence::MoveToEndOfLine)
            || event->matches(QKeySequence::MoveToStartOfBlock)
            || event->matches(QKeySequence::MoveToEndOfBlock)
            || event->matches(QKeySequence::MoveToStartOfDocument)
            || event->matches(QKeySequence::MoveToEndOfDocument)
            || event->matches(QKeySequence::SelectNextChar)
            || event->matches(QKeySequence::SelectPreviousChar)
            || event->matches(QKeySequence::SelectNextWord)
            || event->matches(QKeySequence::SelectPreviousWord)
            || event->matches(QKeySequence::SelectStartOfLine)
            || event->matches(QKeySequence::SelectEndOfLine)
            || event->matches(QKeySequence::SelectStartOfBlock)
            || event->matches(QKeySequence::SelectEndOfBlock)
            || event->matches(QKeySequence::SelectStartOfDocument)
            || event->matches(QKeySequence::SelectEndOfDocument)
            || event->matches(QKeySequence::Cut))
        {

            return;
        }
        // When cursor before command line, reposition cursor at end if paste or regular text.
        if (event->matches(QKeySequence::Paste)
            || (!event->text().isNull() && ((event->modifiers() ^ Qt::ShiftModifier) == Qt::ShiftModifier)))
        {

            cursor.movePosition(QTextCursor::End);
            setTextCursor(cursor);
        }
    }
    // If enter or return, do command
    if (event->key() == Qt::Key_Return || event->key() == Qt::Key_Enter)
    {

        // Select command
        cursor.setPosition(commandStart);
        cursor.movePosition(QTextCursor::End, QTextCursor::KeepAnchor);

        // Convert breaks to \n
        QString line = cursor.selectedText();
        line.replace(QChar::ParagraphSeparator, '\n');

        // Move cursor to end of document before processing key
        cursor.movePosition(QTextCursor::End);
        setTextCursor(cursor);

        // Process event with default handling
        QTextEdit::keyPressEvent(event);

        // Send line in command signal
        command(line);

    }
    else if (cursor.position() > commandStart && event->matches(QKeySequence::MoveToStartOfLine))
    {

        // Home key sets position to start of command only if cursor located after prompt
        cursor.setPosition(commandStart);
        setTextCursor(cursor);

    }
    else
    {

        // Process event with default handling
        QTextEdit::keyPressEvent(event);

    }
}

void Console::message(const QString &text)
{
    // Move the cursor to the end for inserting and after insertion.
    QTextCursor cursor = textCursor();
    cursor.movePosition(QTextCursor::End);

    insertPlainText(text);

    cursor.movePosition(QTextCursor::End);
    setTextCursor(cursor);
    ensureCursorVisible();
}

void Console::prompt(const QString &prompt)
{
    message(prompt);
    commandStart = textCursor().position();
}
