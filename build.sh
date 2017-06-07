#!env bash

set -e
if [ "$(uname)" == "MINGW64_NT-6.1" ]; then
echo Building for MINGW64
pushd configure >/dev/null
echo Building configure
ninja
popd >/dev/null
echo Configuring
./configure.exe
echo Building MIFit
ninja uics
ninja
fi
