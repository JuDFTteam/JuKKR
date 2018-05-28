#!/usr/bin/env bash

# first merge *.dyn files to get pgopti.dpi
profmerge
# then use static (*.spi) and dynamic (*.dpi) information to generate code coverage report
Codecov -spi pgopti.spi -dpi pgopti.dpi

