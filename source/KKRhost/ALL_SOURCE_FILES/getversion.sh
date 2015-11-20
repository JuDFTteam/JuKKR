#!/bin/bash
echo "      character(len=*), parameter :: version='" > tmpver
git describe >> tmpver
echo "'" >> tmpver
tr -d '\n' < tmpver > ALL_SOURCE_FILES/version
echo "" >> ALL_SOURCE_FILES/version
echo "Version: "
git describe
rm -f tmpver

