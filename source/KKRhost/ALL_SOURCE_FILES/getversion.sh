#!/bin/bash
echo "      character(len=*), dimension(4), parameter ::" > tmpver
echo "     +version=(/'" >> tmpver
echo -n "     +" >> tmpver;git describe >> tmpver
echo -n "     +', '" >> tmpver
cat compver >> tmpver
echo -n "     +', '" >> tmpver
cat compflag >> tmpver
echo -n "     +', '" >> tmpver
cat complib >> tmpver
echo "     +'/)" >> tmpver
cat tmpver > ALL_SOURCE_FILES/version
rm -f tmpver compver compflag complib

