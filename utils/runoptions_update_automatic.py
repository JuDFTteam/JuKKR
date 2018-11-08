def get_source_files(folder):
	from os import listdir
	all_source_files = listdir(folder)
	return all_source_files

def construct_keywords_dict():
	listtxt = open('options_list.txt').readlines()
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



if __name__ == '__main__':
	from os.path import sep

	#manual settings
	srcfolder = 'sourcefiles'
	outfolder = 'sourcefiles_new'

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