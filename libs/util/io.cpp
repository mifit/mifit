#include "io.h"

io*io::defaultIo_ = NULL;

io::io()
{
}

io::~io()
{
}

void io::setDefaultIo(io *defaultIo)
{
    if (defaultIo_ != NULL)
    {
        delete defaultIo_;
    }
    defaultIo_ = defaultIo;
}

io*io::defaultIo()
{
    return defaultIo_->create();
}
