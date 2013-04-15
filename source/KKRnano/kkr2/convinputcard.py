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
        result[splitted[0]] = '  '.join(getNums(key=splitted[1], lines=lines, 
                                      num=int(splitted[2]), offset=int(splitted[3])))
    return result
    
keydict = getKeyDict(KEYS.split('\n'), lines)

def apply_kkrnano_rules(keydict):
    newdict = keydict.copy()
    if int(keydict['imix']) > 3:
        mixing = keydict['brymix']
    else:
        mixing = keydict['strmix']
    newdict['mixing'] = mixing
    #del newdict['strmix']
    #del newdict['brymix']
    newdict['nsra'] = str(int(keydict['kvrel']) + 1)
    #del newdict['kvrel']
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

f = open("input.conf", "w")
f.write(template)
f.close()
    
#f = open("inputcard.template", "r")
#template = f.read()
#f.close()
#print substTemplate(template, keydict, True)
