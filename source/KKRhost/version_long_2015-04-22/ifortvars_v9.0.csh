# Run as "source ifortvars....."
#! /bin/tcsh

if !($?PATH) then
    setenv PATH /usr/local/intel/fc/9.0/bin
else
    setenv PATH /usr/local/intel/fc/9.0/bin:${PATH}
endif

if !($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH /usr/local/intel/fc/9.0/lib
else
    setenv LD_LIBRARY_PATH /usr/local/intel/fc/9.0/lib:${LD_LIBRARY_PATH}
endif

# DYLD_LIBRARY_PATH is used on MAC OS
if !($?DYLD_LIBRARY_PATH) then
    setenv DYLD_LIBRARY_PATH /usr/local/intel/fc/9.0/lib
else
    setenv DYLD_LIBRARY_PATH /usr/local/intel/fc/9.0/lib:${DYLD_LIBRARY_PATH}
endif

if !($?MANPATH) then
    setenv MANPATH /usr/local/intel/fc/9.0/man:`manpath`
else
    setenv MANPATH /usr/local/intel/fc/9.0/man:${MANPATH}
endif

if !($?INTEL_LICENSE_FILE) then
    setenv INTEL_LICENSE_FILE "/usr/local/intel/fc/9.0/licenses:/opt/intel/licenses:${HOME}/intel/licenses:/Users/Shared/Library/Application Support/Intel/Licenses"
else
    setenv INTEL_LICENSE_FILE "/usr/local/intel/fc/9.0/licenses:/opt/intel/licenses:${HOME}/intel/licenses:/Users/Shared/Library/Application Support/Intel/Licenses:${INTEL_LICENSE_FILE}"
endif


