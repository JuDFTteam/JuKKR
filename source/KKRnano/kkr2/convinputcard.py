#!/usr/bin/python

"""Read values of keys from Juelich KKR inputcards"""

import re

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
           
           
def getKeyDict(keyfilename, lines):
    result = {}
    keyfile = open(keyfilename, 'r')
    for k in keyfile:
        if not k or k.find('#') == 0:
            continue
        splitted = k.split()
        result[splitted[0]] = '  '.join(getNums(key=splitted[1], lines=lines, 
                                      num=int(splitted[2]), offset=int(splitted[3])))
    keyfile.close()
    return result
    
keydict = getKeyDict('keys.dat', lines)

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

f = open("input.template", "r")
template = f.read()
f.close()
template = substTemplate(template, keydict, False)

f = open("input.conf", "w")
f.write(template)
f.close()
    
#f = open("inputcard.template", "r")
#template = f.read()
#f.close()
#print substTemplate(template, keydict, True)
