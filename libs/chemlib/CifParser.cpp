#include <cstdio>
#include <cstring>

#ifdef _MSC_VER
//Disable annoying warning about max size of symbol being 255 char
#pragma warning(disable: 4786)

#endif

#include "CifParser.h"

/////////////////////////////////////////////////////////////////////////////
// Function:    CifParser (constructor)
// Purpose:		Constructs the parser and binds it to the given FILE ptr
// Input:       FILE pointer to mmCIF file
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
CifParser::CifParser(FILE* fp) {
  _file = fp;
  _start = true;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetNextBlock
// Purpose:
// Input:       None
// Output:      Stores a C string with the next block in the given CifDataBlock
// Requires:	The CIF file is delimited by lines starting with "data_"
/////////////////////////////////////////////////////////////////////////////
bool CifParser::GetNextBlock(CifDataBlock& block) {
  static const int MAX_CIF_LINE = 180;      //The spec is technically 80

  block.Clear();
  std::string buffer;
  char line[MAX_CIF_LINE];
  bool found = false;

  //Advance to the start of the data block, searching
  //for a line starting with "data_"
  if (_start) {
    while (fgets(line, MAX_CIF_LINE, _file) != 0) {
      if (!strncmp(line, "data_", 5)) {
        found = true;
        _start = false;
        _data = line;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }

  //Load block line-by-line into buffer, until we can't
  //read any more lines or we hit another line starting with "_data"
  while (fgets(line, MAX_CIF_LINE, _file) != 0) {
    if (!strncmp(line, "data_", 5)) {
      _data = line;
      block.SetInput(buffer.c_str());       //Pass the lines to the CifDataBlock
      return true;
    } else {
      buffer += line;                   //Load the line into the buffer
      found = true;
    }
  }

  //If we reach here, we're at the end of the file

  if (ferror(_file) || found == false) {
    return false;
  } else {
    block.SetInput(buffer.c_str());     //Pass the lines to the CifDataBlock
    return true;
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    CifTokenizer (constructor)
// Purpose:		Constructs the CifDataBlock tokenizer and binds to an
//				existing string containing a CifDataBlock
// Input:       std::string with the data block
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
CifTokenizer::CifTokenizer(std::string input) {
  _str = input;
  _length = _str.length();
  _cur_pos = 0;
}

//A set of state variables for parsing CIF files
#define PT_INTOKEN 1
#define PT_BETWEENTOKENS 2
#define PT_SINGLEQUOTE 3
#define PT_DOUBLEQUOTE 4
#define PT_MULTILINE 5                          //delimited by semi-colons

/////////////////////////////////////////////////////////////////////////////
// Function:    GetToken (constructor)
// Purpose:		Read a token from the cifDataBlock, if possible, and advance
//				the position to the end of the token
// Input:       None
// Output:      bool for whether the a valid token could be read
//				If returning true, a std::string containing the next token
// Requires:
/////////////////////////////////////////////////////////////////////////////
bool CifTokenizer::GetToken(std::string& token) {
  token.clear();
  int state = PT_BETWEENTOKENS;

  if (_length == 0 || _cur_pos >= _length) {
    return false;
  }

  for (size_t i = _cur_pos; i < _length; ++i) {
    switch (state) {
      case PT_INTOKEN:
        if (IsCifWhitespace(_str[i])) {
          token = _str.substr(_cur_pos, i - _cur_pos);
          _cur_pos = i;
          return true;
        }
        break;

      case PT_BETWEENTOKENS:
        if (IsCifWhitespace(_str[i])) {
          continue;
        } else if (_str[i] == '\'') {
          state = PT_SINGLEQUOTE;
          _cur_pos = i + 1;
        } else if (_str[i] == '"') {
          state = PT_DOUBLEQUOTE;
          _cur_pos = i + 1;
        } else if (_str[i] == ';'
                   && i > 0                     //ensure that i-1 is valid subscript
                   && _str[i-1] == '\n') {
          state = PT_MULTILINE;
          _cur_pos = i + 1;
        } else if (_str[i] == '#') {
          while (i < _length && _str[i] != '\n') {          //advance to end of comment
            i++;
          }
        } else {
          state = PT_INTOKEN;
          _cur_pos = i;
        }
        break;

      case PT_SINGLEQUOTE:
        if (_str[i] == '\''
            && (i+1 == _length
                || IsCifWhitespace(_str[i+1]))) {
          token = _str.substr(_cur_pos, i - _cur_pos);
          _cur_pos = i + 1;
          return true;
          //				token = std::string(_str + tok_start, i - tok_start);
          //				_str += i+1;
          //				return true;
        } else if (i+1 == _length) {
          state = PT_INTOKEN;
          _cur_pos--;
          i = _cur_pos;
        }
        break;

      case PT_DOUBLEQUOTE:
        if (_str[i] == '"'
            && (i+1 == _length
                || IsCifWhitespace(_str[i+1]))) {
          token = _str.substr(_cur_pos, i - _cur_pos);
          _cur_pos = i + 1;
          return true;
          //				(IsCifWhitespace(_str[i+1]) ||
          //				i + 1 == n)) {
          //				token = std::string(_str + tok_start, i - tok_start);
          //				_str += i+1;
          //				return true;
        } else if (i+1 == _length) {
          state = PT_INTOKEN;
          _cur_pos--;
          i = _cur_pos;
          //				tok_start--;
          //				i = tok_start - 1;
        }
        break;

      case PT_MULTILINE:
        if (i == _cur_pos && IsCifWhitespace(_str[i])) {
          _cur_pos++;
          continue;
        } else if (_str[i] == ';' && _str[i-1] == '\n') {
          token = _str.substr(_cur_pos, i - _cur_pos - 1);
          _cur_pos = i + 1;
          return true;
        }
        break;

    }     //switch(state)
  }  //for(i=0...)

  //We've reached the end of the block
  switch (state) {
    case PT_BETWEENTOKENS:
      _cur_pos = _length;
      return false;
    case PT_INTOKEN:
      token = _str.substr(_cur_pos);
      _cur_pos = _length;
      return true;
    default:
      return false;
  }
}

void CifTokenizer::SlurpValues(std::string& input) {
  //	int tok_start;
  //	char *values;
  //	int state = PT_BETWEENTOKENS;
  //	int n = strlen(_str);

  input.clear();
  int state = PT_BETWEENTOKENS;

  if (_length == 0 || _cur_pos >= _length) {
    return;
  }

  for (size_t i = _cur_pos; i < _length; ++i) {
    switch (state) {
      case PT_INTOKEN:
        if (IsCifWhitespace(_str[i])) {
          state = PT_BETWEENTOKENS;
        }
        //			else if (i == n) {
        //				input = std::string(_str, i);
        //				return;
        //			}
        break;

      case PT_BETWEENTOKENS:
        if (IsCifWhitespace(_str[i])) {
          continue;
        } else if (_str[i] == '_') {
          if (i > _cur_pos) {
            input = _str.substr(_cur_pos, i - _cur_pos - 1);
          } else {
            input = "";
          }
          _cur_pos = (i > 0) ? i-1 : 0;
          return;
          //				values = new char[i+1];
          //				strncpy(values, _str, i);
          //				values[i] = 0;
          //				_str += i;
          //				return values;
        } else if (_str[i] == 'l'
                   && _str[i+1] == 'o'
                   && _str[i+2] == 'o'
                   && _str[i+3] == 'p'
                   && _str[i+4] == '_') {
          input = _str.substr(_cur_pos, i - _cur_pos - 1);
          _cur_pos = (i > 0) ? i-1 : 0;
          return;
          //				input = std::string(_str,i);
        } else if (_str[i] == '\'') {
          state = PT_SINGLEQUOTE;
          //				_cur_pos = (i + 1 < _length) ? i + 1 : _length - 1;
        } else if (_str[i] == '"') {
          state = PT_DOUBLEQUOTE;
          //				_cur_pos = (i + 1 < _length) ? i + 1 : _length - 1;;
        } else if (_str[i] == '#') {                    //advance to the end of the comment
          while (i < _length && _str[i] != '\n') {
            i++;
          }
        } else if (_str[i] == ';'
                   && (i == 0                           //ensure that i-1 is a valid subscript
                       || _str[i-1] == '\n')) {
          state = PT_MULTILINE;
          //				_cur_pos = (i + 1 < _length) ? i + 1 : _length - 1;;
        }
        //			else if (i == n) {
        //				input = std::string(_str, i);
        //				return;
        //			}
        else {
          state = PT_INTOKEN;
        }
        break;

      case PT_SINGLEQUOTE:
        if (_str[i] == '\''
            && (i+1 == _length
                || IsCifWhitespace(_str[i+1]))) {
          state = PT_BETWEENTOKENS;
        }
        break;

      case PT_DOUBLEQUOTE:
        if (_str[i] == '"'
            && (i+1 == _length
                || IsCifWhitespace(_str[i+1]))) {
          state = PT_BETWEENTOKENS;
        }
        break;

      case PT_MULTILINE:
        //			if(i == _cur_pos && IsCifWhitespace(_str[i])) {
        //				_cur_pos++;
        //				continue;
        //			}
        if (_str[i] == ';' && _str[i-1] == '\n') {
          //				_cur_pos = (i + 1 < _length) ? i + 1 : _length - 1;;
          state = PT_BETWEENTOKENS;
        }
        break;

    }     //switch(case)
  }  //for(i=0...)

  input = _str.substr(_cur_pos);
  _cur_pos = _length;
  return;
}

void CifTokenizer::SlurpNames(std::vector<std::string>& names) {
  names.clear();
  int state = PT_BETWEENTOKENS;

  if (_length == 0 || _cur_pos >= _length) {
    return;
  }

  for (size_t i = _cur_pos; i < _length; ++i) {
    switch (state) {
      case PT_INTOKEN:
        if (IsCifWhitespace(_str[i])) {
          if (i == _cur_pos) {
            break;
          }
          //				names.Add( std::string( _str + tok_start, i - tok_start ) );
          names.push_back(_str.substr(_cur_pos, i - _cur_pos) );
          _cur_pos = i;
          state = PT_BETWEENTOKENS;
        }
        break;

      case PT_BETWEENTOKENS:
        if (IsCifWhitespace(_str[i])) {
          continue;
        } else if (_str[i] == '_') {
          _cur_pos = i+1;
          state = PT_INTOKEN;
        } else if (_str[i] == '#') {                        //Comments only end with newline
          while (i < _length && _str[i] != '\n') {
            i++;
          }
        } else {                        //Normal exit, we reach a token that's not
          _cur_pos = i;                 //a name
          return;
        }
        break;
    }
  }
  names.clear();                        //Shouldn't reach here with a proper cif...only if
  _cur_pos = _length;                   //the end of the block or file comes unexpectedly
                                        //in the middle of a loop header
}

