#!/usr/bin/env bash

set -e
if [ "$(uname)" == "MINGW64_NT-6.1" ]; then
echo Building for MINGW64
pushd configure >/dev/null
echo Building configure
cp mingw.ninja platform.ninja
ninja
popd >/dev/null
echo Configuring
./configure.exe mingw
echo Building MIFit
ninja uics
ninja

elif [ "$(uname)" == "Linux" ]; then
echo Building for Ubuntu Trusty
pushd configure >/dev/null
echo Building configure
cp trusty.ninja platform.ninja
ninja
popd >/dev/null
echo Configuring
./configure.exe trusty
echo Building MIFit
ninja uics
ninja

fi
