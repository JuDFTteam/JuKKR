import numpy as np
import glob
import os, sys

#####################################################################################
# Create a bar plot showing the categories
#####################################################################################
def cat_bar_plot(code,cat,sub_num):
    import matplotlib.pyplot as plt
    from matplotlib import cm as cm
    import numpy as np
    plt.rc('text', usetex=False)#True)
    plt.rc('font', family='serif')

    font_size=14 #28
    fig_dpi=200 #800
    tol=1e-5

    prune_key=[]
    percent=[]
    total_num=0
    for key in cat:
        if key!=code and key!='undefined':
            prune_key.append(key)
            percent.append(sub_num[cat.index(key)])
            total_num=total_num+sub_num[cat.index(key)] 
    indx=np.argsort(percent)
    sort_per=[]
    sort_key=[]
    for ii in range(len(indx)-1,-1,-1):
        sort_key.append(prune_key[indx[ii]])
        sort_per.append(percent[indx[ii]]*100./total_num)
        #sort_per.append(percent[indx[ii]]) #*100./total_num)
    print(code, sum(sort_per))
    x_axis=np.arange(len(sort_per))
    colors=cm.gnuplot(np.linspace(0,1,len(sort_per)))

    fig = plt.figure(figsize=(16,10))
    ax_per= fig.add_subplot(111)
    ax_per.bar(x_axis,sort_per,color=colors)
    yl = ax_per.get_ylim()
    #ax_per.set_ylim(yl[0]-1, yl[1])
    #plt.axhline(0, color='grey', lw=1)
    #ax_per.set_yscale('log')
    plt.xticks(x_axis,sort_key,rotation='vertical',fontsize=font_size*0.75)
    plt.yticks(fontsize=font_size*0.75)
    plt.ylabel('percentage of calls',fontsize=font_size)
    for axis in ['bottom','left','right','top']:
        ax_per.spines[axis].set_linewidth(3)
    #plt.show()
    FigName='call_pct_'+code+'_.pdf'
    fig.savefig(FigName,transparent=False,dpi=fig_dpi,bbox_inches='tight')
    return

#####################################################################################
# Create a list of files over which one will be looping about
#####################################################################################
def find_files(src_dir):
    ext=["*.f90","*.F90","*.f","*.F"]

    files=[]
    for dirs in os.walk(src_dir):
        print src_dir, dirs[0]
        for ii in range(len(ext)):
            files=files+glob.glob(dirs[0]+ext[ii])

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
    call_type=[]
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
            call_type.append(data[0].lower())
        if len(data)>1 and (data[0].lower()!='end' and data[1].lower()=='function'):
            if data[0]!='!>':
                try:
                    indx=data[2].index('(')
                    sub.append(data[2][0:indx].lower())
                except:
                    sub.append(data[2].lower())
                call_type.append(data[1].lower())
        if len(data)>2 and (data[0].lower()!='end' and data[2].lower()=='function'):
            if data[0]!='!>' and data[0]!='!':
                try:
                    indx=data[3].index('(')
                    sub.append(data[3][0:indx].lower())
                except:
                    sub.append(data[3].lower())
                call_type.append(data[2].lower())
        if len(data)>1 and (data[0].lower()=='type' and data[1].lower()=='::'):
            sub.append(data[2])
            call_type.append(data[0].lower())
    # If there is no module call it no module
    if len(module)==0:
        module.append('no module only file '+filename)
    if len(sub)==0:
        sub.append('no subroutine only file: '+filename)
        call_type.append('misc')
    if len(cat)==0 or len(cat)!=len(sub):
        for ii in range(len(sub)-len(cat)):
           tmp=[code,'undefined']
           cat.append(tmp)
    return sub, cat,module,call_type

#####################################################################################
# Create a dictionary that prints the category for each subroutine
#####################################################################################
def subroutine_dict(mod,sub,cat,call_type,code='kkrhost'):
    # Generate a dictionary for each file
    from ruamel.yaml import YAML
    mod_dict=dict()
    mod_dict['module']=dict()
    for jj in range(0,len(mod)):
        mod_dict['module']['name']=mod[jj]
        mod_dict['module'][mod[jj]]=dict()
        for ii in range(0,len(sub)):
            mod_dict['module'][mod[jj]]['call_type '+str(ii)]=call_type[ii]
            mod_dict['module'][mod[jj]]['call '+str(ii)]=sub[ii]
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
def category_dict(cat,subroutine,call_type,code='kkrhost'):
    from ruamel.yaml import YAML
    # Flatten the list of the categories
    temp=[j for sub in cat for j in sub]
    # Find the unique entries
    unq=list(set(temp))
    unq_call=list(set(call_type))
    num_call=[]
    for key in unq_call:
        num_call.append(call_type.count(key))
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
    return cat_dict,unq,sub_num,num_call,unq_call

#####################################################################################
# 
#####################################################################################
def main():
    name_dirs = [['voronoi', 'voronoi'], ['kkrhost', 'KKRhost'], ['kkrimp', 'KKRimp'], ['pkkprime', 'PKKprime'], ['kkrsusc', 'KKRsusc/solver_module_v2'], ['rhoq', 'rhoq'], ['common', 'common']]
    for code, src_dir in name_dirs:
        code = code.lower()
        src_dir = '../source/'+src_dir+'/'
        files=find_files(src_dir)
 
        total_sub=[]
        total_mod=[]
        total_cat=[]
        total_call=[]
        for filename in files:
            complete_path=src_dir+filename
            sub,cat,mod,call_type=extract_from_file(complete_path,code)
            subroutine_dict(mod,sub,cat,call_type,code)
            total_sub=total_sub+sub
            total_mod=total_mod+mod
            total_cat=total_cat+cat
            total_call=total_call+call_type
        cat_dict,unq,sub_num,num_call,unq_call=category_dict(total_cat,total_sub,total_call,code)
 
        cat_bar_plot(code,unq,sub_num)
        cat_bar_plot('types_'+code,unq_call,num_call)
    return

if __name__ == '__main__':
    main()

