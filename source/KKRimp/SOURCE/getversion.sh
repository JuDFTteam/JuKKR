#!/bin/bash
echo "      character(len=*), parameter :: version='" > tmpver
git describe >> tmpver
echo ":" >> tmpver
cat compver >> tmpver
echo "'" >> tmpver
tr -d '\n' < tmpver > version
echo "" >> version
echo "Version: "
git describe
rm -f tmpver compver

