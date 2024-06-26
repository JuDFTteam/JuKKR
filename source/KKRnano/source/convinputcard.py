#!/usr/bin/python

"""Read values of keys from Juelich KKR inputcards"""

import re

KEYS = (
"""
#new key name | old key name | num. values to read | offset
alat      ALATBASIS 1 0
bravais_a BRAVAIS 3 1
bravais_b BRAVAIS 3 2
bravais_c BRAVAIS 3 3
scfsteps  NSTEPS  1 1
imix      IMIX    1 1
tempr     TEMPR   1 1
fcm       FCM     1 1
gmax      GMAX    1 0
rmax      RMAX    1 0
icst      ICST    1 0
strmix    STRMIX  1 1
kte       KTE     1 0
kforce    KFORCE  1 0
kpre      KPRE    1 0
kxc       KEXCOR  1 0
qmrbound  QMRBOUND 1 1
ldau      LLDAU   1 0
kvrel     KVREL   1 0
rclust    RCLUSTZ 1 0
basisscale BASISCALE 3 0
emin      EMIN    1 1
emax      EMAX    1 1
jij       LJIJ    1 0
brymix    BRYMIX  1 1
rcutjij   RCUTJIJ 1 0
cartesian CARTESIAN 1 0
bzdivide  BZDIVIDE  3 0
npol      NPOL    1 1
npnt1     NPT1    1 1
npnt2     NPT2    1 1
npnt3     NPT3    1 1
""")

template = (
"""
# LATTICE
alat       = $alat
basisscale = $basisscale

#BRAVAIS
bravais_a = $bravais_a
bravais_b = $bravais_b
bravais_c = $bravais_c

cartesian = $cartesian

# number of k-points in each direction
bzdivide = $bzdivide

rclust   = $rclust  # reference cluster radius

# Energy contour
emin  = $emin
emax  = $emax
npnt1 = $npnt1
npnt2 = $npnt2
npnt3 = $npnt3
npol  = $npol           # Number of Matsubara poles, npol=0 triggers DOS calculation
tempr = $tempr

# Self-consistency options
scfsteps = $scfsteps
imix     = $imix
mixing   = $mixing
fcm      = $fcm

target_rms = $target_rms   # abort when target_rms error has been reached

# Parameters for Ewald sums
rmax = $rmax
gmax = $gmax

# Exchange correlation potential
kxc = $kxc

# Solver options
qmrbound = $qmrbound

icst    = $icst      # num. Born iterations for non-spherical potential
kpre    = $kpre
kforce  = $kforce
jij     = $jij
ldau    = $ldau
rcutjij = $rcutjij
nsra    = $nsra      # 1=non-scalar-relativistic 2=scalar-relativistic
kte     = $kte

#------------------------------------------------------------------------------
# Shape-function options
#------------------------------------------------------------------------------

rclust_voronoi = $rclust_voronoi  # radius of cluster used for Voronoi
nmin_panel     = $nmin_panel      # minimum number of points per panel

# number of points for 'muffin-tinization'
# create shape-function for num_MT_points in MT region
# used to restrict core wavefunction to MT region
# choose 0 for touching MT spheres
# suggested value = 10

num_MT_points  = $num_MT_points

# determines how to set new MT-radius, choose 0.0 to get new MT radius from
# atominfo file, otherwise choose 0.0 < MT_scale < 1.0 as factor to scale
# maximal MT radius to new MT radius

MT_scale       = $MT_scale

# Reference potential
# Radius of repulsive reference potential
# RMT_ref_scale <= 0.0: use same radius as for muffin-tin sphere
# otherwise use RMT_ref_scale * (max. possible muffin tin radius)
# Recommended: RMT_ref_scale = 0.98

RMT_ref_scale = $RMT_ref_scale

""")

f = open('inputcard', 'r')
lines = f.readlines()
f.close()
#print lines

class InputCardKeyError:
    pass

def extractNums(line, start, num):
    result = []
    splitted = line[start:].split()

    count = 0
    for entry in splitted:
        if count >= num:
            break
        if entry[0] == '=' and count == 0:
            entry = entry[1:]          
        if entry == '':
            continue
        
        result.append(entry)
        count = count + 1
        
    return result

def getNums(key, lines, num=1, offset=0):
    """
    offset=0: get value(s) from same line after key
    offset>0: get values from offset lines directly below key
    """
    result = []
    for nline, line in enumerate(lines):
        m = re.search(r'\b' + key + r'\b', line)
        if m:
           if offset == 0:
               result = extractNums(line, m.end(), num)
           else:
               result = extractNums(lines[nline+offset], m.start(), num)
    
    if (len(result) < num):
        print "Error with key " + key
        raise InputCardKeyError
    
    return result[:num]
           
           
def getKeyDict(keyfile, lines):
    result = {}
    for k in keyfile:
        if not k or k.find('#') == 0:
            continue
        splitted = k.split()
        
        try:
            key_values = getNums(key=splitted[1], lines=lines, 
                                 num=int(splitted[2]), offset=int(splitted[3]))
        except InputCardKeyError:
            key_values = ["???"]
        
        result[splitted[0]] = '  '.join(key_values)
    return result
    
keydict = getKeyDict(KEYS.split('\n'), lines)

def apply_kkrnano_rules(keydict):
    """Apply KKRnano specific transformation of input parameters
       Provide default values for new parameters not present in old
       'inputcard'"""
    newdict = keydict.copy()
    if int(keydict['imix']) > 3:
        mixing = keydict['brymix']
    else:
        mixing = keydict['strmix']
    newdict['mixing'] = mixing
    newdict['nsra'] = str(int(keydict['kvrel']) + 1)
    # provide default values for new shape-function parameters
    # to behave as in Voronoi-program
    
    # Voronoi cluster 1.6 times larger than ref. cluster
    newdict['rclust_voronoi'] = str(float(newdict['rclust'].lower().replace('d','e')) * 1.6)
    newdict['nmin_panel'] = '3'
    newdict['num_MT_points'] = '10'
    newdict['MT_scale'] = '0.0'
    newdict['RMT_ref_scale'] = '0.0'
    newdict['target_rms'] = '0.0'
    return newdict

keydict = apply_kkrnano_rules(keydict)


def substTemplate(template, keydict, fixed=False):
    """Substitute keys into template"""
    for k in keydict:
        reg = r'\$'  + k + r'\b'
        if fixed:
            subst = keydict[k].ljust(len(k) + 1)  # keep fixed form!!!
        else:
            subst = keydict[k]
        template = re.sub(reg, subst, template)
    return template

template = substTemplate(template, keydict, False)

if __name__ == "__main__":
    f = open("input.conf", "w")
    f.write(template)
    f.close()

