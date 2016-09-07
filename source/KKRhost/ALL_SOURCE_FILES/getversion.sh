#!/bin/bash

# print warning if .git was not found
function noversion {
    printf "WARNING_no_current_version_info:" >> tmpver
    printf "\n"
    printf "\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "!                                                                 !\n"
    printf "!                 WARNING NO VERSION INFO AVAILABLE               !\n"
    printf "!                                                                 !\n"
    printf "!             get access to git repository:                       !\n"
    printf "!             gitlab@iffgit.fz-juelich.de:kkr/kkrjm.git           !\n"
    printf "!                                                                 !\n"
    printf "!             https://iffgit.fz-juelich.de/kkr/kkrjm              !\n"
    printf "!                                                                 !\n"
    printf "!             I will use latest known version number:             !\n"
    cat ALL_SOURCE_FILES/last_version_number
    printf "!                                                                 !\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "\n"
    printf "\n"

    cat ALL_SOURCE_FILES/last_version_number >> tmpver
}

# write current version taken from .git
function write_version {
    git describe >> tmpver && git describe > ALL_SOURCE_FILES/last_version_number || noversion 
}


# gather version info
printf "character(len=*), dimension(4), parameter ::" > tmpver
printf " version=(/ '" >> tmpver
write_version
printf "', '" >> tmpver
cat compver >> tmpver
printf "', '" >> tmpver
cat compflag >> tmpver
printf "', '" >> tmpver
cat complib >> tmpver
printf "' /)" >> tmpver

# delete newlines
tr '\n ' ' ' < tmpver > tmpver_2

# write module
printf "module mod_version\n" > tmpver
printf "\n">> tmpver
printf "implicit none\n">> tmpver
printf "private\n">> tmpver
printf "public version\n">> tmpver
printf "\n">> tmpver
cat tmpver_2 >> tmpver
printf "\n">> tmpver
printf "\n">> tmpver
printf "end module mod_version\n">> tmpver


# write tmpver to version.f90
cp tmpver ALL_SOURCE_FILES/version.f90

# clean up
rm -f tmpver compver compflag complib tmpver_2

