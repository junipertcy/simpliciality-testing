#!/bin/bash

echo $PATH

if [[ "$TRAVIS_OS_NAME" == "linux" ]]
then
    echo "not doing anything..."
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]
then
    brew update
    brew upgrade
elif [[ "$TRAVIS_OS_NAME" == "windows" ]]
then
    choco install python --version 3.8.0
    python --version
fi

export INSTALL_TEST_REQUIREMENTS="true"
echo $INSTALL_TEST_REQUIREMENTS
