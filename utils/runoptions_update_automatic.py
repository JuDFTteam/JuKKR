def get_source_files(folder):
	from os import listdir
	all_source_files = listdir(folder)
	remove_files = ['test','runopt_update_tool','runoptions_update_automatic.py','out','runoptions_list.txt','vinterface.f90','runoptions.F90','rinput13.F90']
	for file in remove_files:
		try:
			all_source_files.remove(file)
		except:
			pass

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

	sources_modified = []

	for srcfile in all_source_files:

		#open source and write backup
		srctxt = open(srcfolder+sep+srcfile,'r').readlines()
		write_src(outfolder+sep+srcfile+'.bak',srctxt)

		counter = 0

		for newkey in keywords:
			for oldkey in keywords[newkey][0]:
				newsrctxt=[]
				compiled_regexp = get_regex_compiled(oldkey)
				for iline, line in enumerate(srctxt):
					newline, replacements = replace_line(line,compiled_regexp,newkey)
					newsrctxt.append(newline)
					total_replacements += replacements
					counter += replacements
					if (replacements>0):
						print('%-20s: %i replacements on line %i with keyword %s'%(srcfile,replacements,iline+1,newkey))
					if (replacements>1):
						WARNINGS.append('%-20s: %i replacements on line %i with keyword %s'%(srcfile,replacements,iline+1,newkey))

				srctxt = newsrctxt

		write_src(outfolder+sep+srcfile,newsrctxt)
		if(counter>0): sources_modified.append(srcfile)

	open(outfolder+sep+'WARNINGS.txt','w').write('\n'.join(WARNINGS))
	print(30*'#')
	print('In total %i replacements in %i files with %i warnings.'%(total_replacements,len(sources_modified),len(WARNINGS)))
	print 'open '+' '.join(sources_modified)


def find_subroutine_blocks(srctxt):
	from re import compile, IGNORECASE
	c_start = compile(r'^\s{0,100}subroutine\s.*\(', flags=IGNORECASE)
	c_end   = compile(r'^\s{0,100}end\s{1,100}subroutine[\s!]', flags=IGNORECASE)

	istart = [iline for iline, line in enumerate(srctxt) if c_start.match(line)!=None]
	iend   = [iline for iline, line in enumerate(srctxt) if c_end.match(line)!=None]

	if(len(istart)!=len(iend)):
		raise ValueError('unequal starts and ends of subroutines are found')

	return zip(istart,iend)

def append_use_statement_to_subroutine():
	from os.path import sep

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


def filter_opt_test_calls():
	from os.path import sep
	from re import compile, IGNORECASE

	#manual settings
	srcfolder='..'
	outfolder = 'out'

	#initialize variables
	#WARNINGS=[]
	#keywords = construct_keywords_dict(filename='runoptions_list.txt')
	all_source_files = get_source_files(srcfolder)

	sources_modified = []

	c_8_opt_call  = compile(r'opt {0,10}\( {0,10}["'+"'"+'].{8,8}["'+"'"+r'] {0,10}\)', flags=IGNORECASE)
	c_8_test_call = compile(r'test {0,10}\( {0,10}["'+"'"+'].{8,8}["'+"'"+r'] {0,10}\)', flags=IGNORECASE)
	c_any_opt_call  = compile(r'opt {0,10}\( {0,10}["'+"'"+'].*["'+"'"+r'] {0,10}\)', flags=IGNORECASE)
	c_any_test_call = compile(r'test {0,10}\( {0,10}["'+"'"+'].*["'+"'"+r'] {0,10}\)', flags=IGNORECASE)

	o_8_opt_call = open(outfolder+sep+'c_8_opt_call.txt','w')
	o_8_test_call = open(outfolder+sep+'c_8_test_call.txt','w')
	o_any_opt_call = open(outfolder+sep+'c_any_opt_call.txt','w')
	o_any_test_call = open(outfolder+sep+'c_any_test_call.txt','w')

	for srcfile in all_source_files:

		#open source and write backup:q	
		srctxt = open(srcfolder+sep+srcfile,'r').readlines()

		for iline, line in enumerate(srctxt):
			if(len(c_8_opt_call.findall(line))>0): o_8_opt_call.write('%-30s: %5i : %s\n'%(srcfile,iline+1,line.strip()))
			if(len(c_8_test_call.findall(line))>0): o_8_test_call.write('%-30s: %5i : %s\n'%(srcfile,iline+1,line.strip()))
			if(len(c_any_opt_call.findall(line))>0): o_any_opt_call.write('%-30s: %5i : %s\n'%(srcfile,iline+1,line.strip()))
			if(len(c_any_test_call.findall(line))>0): o_any_test_call.write('%-30s: %5i : %s\n'%(srcfile,iline+1,line.strip()))


if __name__ == '__main__':
	#automatic_replacement_all_sources()
	#append_use_statement_to_subroutine()
	#filter_opt_test_calls()

	from os.path import sep
	from re import compile, IGNORECASE

	cs=[]

	cs.append(compile(r'logical.*opt', flags=IGNORECASE))
	cs.append(compile(r'external.*opt', flags=IGNORECASE))
	cs.append(compile(r'logical.*test', flags=IGNORECASE))
	cs.append(compile(r'external.*test', flags=IGNORECASE))


	#manual settings
	srcfolder='..'
	outfolder = '..'

	#initialize variables
	WARNINGS=[]
	all_source_files = get_source_files(srcfolder)

	sources_modified = []
	outfile = open(outfolder+sep+'testopt_all.txt','w')
	for srcfile in all_source_files:

		#load source and write backup
		srctxt = open(srcfolder+sep+srcfile,'r').readlines()
		write_src(outfolder+sep+srcfile+'.bak',srctxt)

		delete_lines = []

		#scan lines which shall be removel
		for iline, line in enumerate(srctxt):
			ls = any([c.search(line)!=None for c in cs])
			if(ls):
				outfile.write('%-30s: %5i : %s\n'%(srcfile,iline+1,line.strip()))
				delete_lines.append(iline)

		#remove lines
		for idel, iline in enumerate(delete_lines):
			srctxt.pop(iline-idel)

		#write the source file
		write_src(outfolder+sep+srcfile,srctxt)

	outfile.close()