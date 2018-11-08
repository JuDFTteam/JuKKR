def get_source_files(folder):
	from os import listdir
	all_source_files = listdir(folder)
	all_source_files.remove('test')
	all_source_files.remove('runopt_update_tool')
	return all_source_files

def construct_keywords_dict(filename='runoptions_list.txt'):
	listtxt = open(filename).readlines()
	keywords =dict()

	for line in listtxt:
		new_keyword=line.split('::')[1].split('=')[0].strip()
		description=line.split('!!')[1].strip()
		old_keywords = [word.strip().split("'")[1] for word in line.split('former:')[1].strip()[:-1].split(' and ')]
		keywords[new_keyword] = [old_keywords, description]

	return keywords

def write_src(filename,txtlines):
	open(filename,'w').write(''.join(txtlines))

def get_regex_compiled(oldkey):
	from re import compile, IGNORECASE
	oldkey8='%-8s'%(oldkey)
	regexp_string = r' {0,10}\( {0,10}["'+"'"+']'+oldkey8+'["'+"'"+r'] {0,10}\)'
	compiled_regexp=[ compile(modestr+regexp_string, flags=IGNORECASE) for modestr in ['opt','test']]
	return compiled_regexp

def replace_line(line,compiled_regexp,newkey):
	counter = 0
	for c in compiled_regexp:
		line, icount = c.subn(newkey,line)
		counter += icount
	return line, counter

def automatic_replacement_all_sources():
	from os.path import sep

	#manual settings
	srcfolder = '.'
	outfolder = '.'

	#initialize variables
	WARNINGS=[]
	total_replacements = 0
	keywords = construct_keywords_dict()
	all_source_files = get_source_files(srcfolder)


	for srcfile in all_source_files:

		#open source and write backup
		srctxt = open(srcfolder+sep+srcfile,'r').readlines()
		write_src(outfolder+sep+srcfile+'.bak',srctxt)

		for newkey in keywords:
			for oldkey in keywords[newkey][0]:
				newsrctxt=[]
				compiled_regexp = get_regex_compiled(oldkey)
				for iline, line in enumerate(srctxt):
					newline, replacements = replace_line(line,compiled_regexp,newkey)
					newsrctxt.append(newline)
					total_replacements += replacements
					if (replacements>0):
						print('%-20s: %i replacements on line %i with keyword %s'%(srcfile,replacements,iline+1,newkey))
					if (replacements>1):
						WARNINGS.append('%-20s: %i replacements on line %i with keyword %s'%(srcfile,replacements,iline+1,newkey))

				srctxt = newsrctxt

		write_src(outfolder+sep+srcfile,newsrctxt)

	open(outfolder+sep+'WARNINGS.txt','w').write('\n'.join(WARNINGS))
	print(30*'#')
	print('In total %i replacements in %i files with %i warnings.'%(total_replacements,len(all_source_files),len(WARNINGS)))


def find_subroutine_blocks(srctxt):
	from re import compile, IGNORECASE
	c_start = compile(r'^\s{0,100}subroutine\s.*\(', flags=IGNORECASE)
	c_end   = compile(r'^\s{0,100}end\s{1,100}subroutine[\s!]', flags=IGNORECASE)

	istart = [iline for iline, line in enumerate(srctxt) if c_start.match(line)!=None]
	iend   = [iline for iline, line in enumerate(srctxt) if c_end.match(line)!=None]

	if(len(istart)!=len(iend)):
		raise ValueError('unequal starts and ends of subroutines are found')

	return zip(istart,iend)




if __name__ == '__main__':
	#automatic_replacement_all_sources()

	from os.path import sep
	from re import compile, IGNORECASE

	#manual settings
	srcfolder='.'
	outfolder = '.'

	#initialize variables
	WARNINGS=[]
	keywords = construct_keywords_dict(filename='runoptions_list.txt')
	all_source_files = get_source_files(srcfolder)

	sources_modified = []


	for srcfile in all_source_files:

		#open source and write backup
		srctxt = open(srcfolder+sep+srcfile,'r').readlines()
		write_src(outfolder+sep+srcfile+'.bak',srctxt)

		#find subroutine blocks
		try:
			sub_blocks = find_subroutine_blocks(srctxt)
		except ValueError:
			WARNINGS.append('Inspect %s'%(srcfile))
			continue

		#scan through subroutines and filter runoption keywords
		sub_blocks_keys_found = []
		for sub in sub_blocks:
			tmp_keys_found = []
			for line in srctxt[sub[0]:sub[1]]:
				for key in keywords.keys():
					if(key in line): tmp_keys_found.append(key)
			sub_blocks_keys_found.append(set(tmp_keys_found))

		#construct the module-load statement and add after subroutine definition (needs manual adapration afterwards)
		counter = 0
		for isub, sub in enumerate(sub_blocks):
			keys = sub_blocks_keys_found[isub]
			if(len(keys)>0):
				#indenttxt = srctxt[sub[0]+counter].split('subroutine')[0]
				#addline =  indenttxt+'  use :: mod_runoptions, only: '+', '.join(sorted(list(keys)))+' --manopt-- \n'
				addline = 'use :: mod_runoptions, only: '+', '.join(sorted(list(keys)))+' --manopt-- \n'
				srctxt.insert(sub[0]+counter+1,addline)
				counter += 1

		if(counter>0): sources_modified.append(srcfile)


		#write the source file
		write_src(outfolder+sep+srcfile,srctxt)


	open(outfolder+sep+'WARNINGS.txt','w').write('\n'.join(WARNINGS))
	print 'open '+' '.join(sources_modified)