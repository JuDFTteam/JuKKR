import numpy as np
import glob
import os, sys

def find_files(src_dir):
    ext=["*.f90","*.F90","*.f","*.F"]

    files=[]
    for ii in range(len(ext)):
        files=files+glob.glob(src_dir+ext[ii])

    out_files=[]
    for ii in range(len(files)):
        out_files.append(os.path.basename(files[ii]))
    return out_files
#####################################################################################
# Function to extract the categories and subroutines from each file
#####################################################################################
def extract_from_file(filename,code='kkrhost'):
    mod_file=open(filename)
    cat=[]
    sub=[]
    module=[]
    # Go through the file line by line finding the occurrences of the keywords
    # "category","subroutine","module" and apppend them to the appropriate arrays
    for line in mod_file:
        data = str.split(line)
        if len(data)>1 and data[0].lower()=='module':
            module.append(data[1].lower())
        if len(data)>2 and (data[1].lower()=='category:' or data[1].lower()=='kkrtags:'):
            tmp=[]
            for ii in range(2,len(data)):
                # Eliminate any commas in the categories
                tmp_str=data[ii].replace(",","")
                tmp.append(tmp_str.lower())
            cat.append(tmp)
        if len(data)>1 and (data[0].lower()=='subroutine' or data[0].lower()=='function' or data[0].lower()=='program'):
            # Get rid of any "(" in the subroutines
            try:
                indx=data[1].index('(')
                sub.append(data[1][0:indx].lower())
            except:
                sub.append(data[1].lower())
        if len(data)>1 and (data[0].lower()!='end' and data[1].lower()=='function'):
            if data[0]!='!>':
                try:
                    indx=data[2].index('(')
                    sub.append(data[2][0:indx].lower())
                except:
                    sub.append(data[2].lower())
        if len(data)>2 and (data[0].lower()!='end' and data[2].lower()=='function'):
            if data[0]!='!>' and data[0]!='!':
                try:
                    indx=data[3].index('(')
                    sub.append(data[3][0:indx].lower())
                except:
                    sub.append(data[3].lower())
        if len(data)>1 and (data[0].lower()=='type' and data[1].lower()=='::'):
            sub.append(data[2])
    # If there is no module call it no module
    if len(module)==0:
        module.append('no module only file '+filename)
    if len(sub)==0:
        sub.append('no subroutine only file: '+filename)
    if len(cat)==0 or len(cat)!=len(sub):
        for ii in range(len(sub)-len(cat)):
           tmp=[code,'undefined']
           cat.append(tmp)
    return sub, cat,module

#####################################################################################
# Create a dictionary that prints the category for each subroutine
#####################################################################################
def subroutine_dict(mod,sub,cat,code='kkrhost'):
    # Generate a dictionary for each file
    from ruamel.yaml import YAML
    mod_dict=dict()
    mod_dict['module']=dict()
    for jj in range(0,len(mod)):
        mod_dict['module']['name']=mod[jj]
        mod_dict['module'][mod[jj]]=dict()
        for ii in range(0,len(sub)):
            mod_dict['module'][mod[jj]]['subroutine '+str(ii)]=sub[ii]
            mod_dict['module'][mod[jj]][sub[ii]]=dict()
            mod_dict['module'][mod[jj]][sub[ii]]=cat[ii]

    yaml = YAML()
    yaml.indent(mapping=4, sequence=6, offset=3)
    with open('KKR_dictionary_'+code+'_.yml', 'a') as outfile:
        yaml.dump(mod_dict, outfile)
    return mod_dict

#####################################################################################
# Takes the unique categories and then assigns the subroutine that belongs to each 
# category
#####################################################################################
def category_dict(cat,subroutine,code='kkrhost'):
    from ruamel.yaml import YAML
    # Flatten the list of the categories
    temp=[j for sub in cat for j in sub]
    # Find the unique entries
    unq=list(set(temp))
    # Loop over each category
    cat_to_sub=[]
    sub_num=[]
    cat_dict=dict()
    cat_dict['category']=dict()
    for ii in range(0,len(unq)):
        sub_cat=[]
        for jj in range(0,len(subroutine)):
            if unq[ii] in cat[jj]:
                sub_cat.append(subroutine[jj])
        sub_num.append(len(sub_cat))
        cat_to_sub.append(sub_cat)
        cat_dict['category'][unq[ii]]=dict()
        cat_dict['category'][unq[ii]]['num_subroutines']=sub_num[ii]
        cat_dict['category'][unq[ii]]['subroutines']=cat_to_sub[ii]
    yaml = YAML()
    yaml.indent(mapping=4, sequence=6, offset=3)
    with open('KKR_cat_dictionary_'+code+'_.yml', 'w') as outfile:
        yaml.dump(cat_dict, outfile)
    return cat_dict

#####################################################################################
# 
#####################################################################################
def main():
    src_dir='./'
    code='PKKprime'
    files=find_files(src_dir)

    total_sub=[]
    total_mod=[]
    total_cat=[]

    for filename in files:
        complete_path=src_dir+filename
        sub,cat,mod=extract_from_file(complete_path,code)
        print filename
        print 'sub'
        print sub
        print 'cat'
        print cat
        print 'mod'
        print mod
        subroutine_dict(mod,sub,cat,code)
        total_sub=total_sub+sub
        total_mod=total_mod+mod
        total_cat=total_cat+cat

    category_dict(total_cat,total_sub,code)
    return

if __name__ == '__main__':
    main()

