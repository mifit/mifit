#ifndef MMCIF_PARSER_H
#define MMCIF_PARSER_H

#include <string>

#include "CifData.h"

//Class to parse a file in mmCif format into individual
//data blocks (denoted by keywords of the form _data_comp_XX)
class CifParser
{
private:
    FILE *_file;
    bool _start;                                // Flag for when we're in the header
    std::string _data;                              // Line of input starting the block (e.g. "data_2SA")

public:
    CifParser(FILE *fp);
    bool GetNextBlock(CifDataBlock &block);
};

//Class to parse an individual data block in mmCif format into individual
//tokens (names, values, keywords)
class CifTokenizer
{
    std::string _str;
    size_t _cur_pos;
    size_t _length;
public:
    CifTokenizer(std::string input);
    //	CifTokenizer(const char *input);
    //	CifTokenizer::~CifTokenizer();

    bool  GetToken(std::string&);

    //Used when a parsing the list of values in a loop
    inline bool GetString(std::string &token)
    {
        return GetToken(token);
    }

    void SlurpNames(std::vector<std::string> &names);
    void SlurpValues(std::string &values);
};

//Macro to define whitespace characters that delimit CIF tokens
inline bool IsCifWhitespace(char c)
{
    return (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f');
}

#endif //MMCIF_PARSER_H
