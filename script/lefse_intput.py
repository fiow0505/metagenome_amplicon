# Author : Kyeongeui Yun (fiow3250@gmail.com)
# 2023.06.02 

import os
import pandas as pd
import sys
from optparse import OptionParser as OP

desc = '''
Program description:
	Generating LEfSe inpput file
'''

debug = 0
def cmdparameter(argv):
	if len(argv) == 1:
		global desc
		#print(desc, file=sys.stderr)
		cmd = 'python ' + argv[0] + ' -h'
		os.system(cmd)
		sys.exit(1)

	usages = "%prog -i taxonomy-levelAll.RelAbundance.tsv -m metadata.txt"
	parser = OP(usage=usages)
	parser.add_option("-i", "--input-file", dest="input_file",
                      help="'taxonomy-levelAll.RelAbundance.tsv' all level file")
	parser.add_option("-m", "--metadata", dest="metadata", 
                      help="metadata file")
	parser.add_option("-n", "--group", dest="group", 
                      help="Group column in metdata")
	parser.add_option("-s", "--sortvalue", dest="sortvalue", 
                      help="sort by group value.  ex.Group1,Group2,Group3")
	parser.add_option("-o", "--output-file", dest="output_file", 
                      help="Output filename [default: lefse_in.txt]")
	parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
     
	(options, args) = parser.parse_args(argv[1:])
	assert options.input_file != None, "A filename needed for -i"
	return (options, args)

def main():
	global debug
	options, args = cmdparameter(sys.argv)
	debug = options.debug
	relevfile = options.input_file
	metafile = options.metadata
	groupID = options.group
	sortvalue = options.sortvalue
	output_file = "lefse_in.txt" if options.output_file is None else options.output_file
	output_path = os.path.dirname(output_file)
	if os.path.exists(output_path)==False:
		os.makedirs(output_path)
    
	allrelev = pd.read_csv(relevfile, sep="\t", index_col=0)
	success_samples = allrelev.columns.tolist()
	allrelev.index = allrelev.index.str.replace(";", "|").str.replace(" ", "_")
	metadata = pd.read_csv(metafile, sep="\t", index_col=0, dtype={"sample-id":str}).rename(columns = {groupID:"group_type"})
	metadata.index = metadata.index.astype(str)
	groupvalue = metadata["group_type"].unique().tolist() if sortvalue==None else sortvalue.split(",")
	metadata = metadata[metadata["group_type"].isin(groupvalue)][["group_type"]]
	samples = metadata["group_type"].index.tolist()
	samples = list(set(samples)&set(success_samples))
	samples_allrelev = allrelev[samples]
	samples_allrelev = samples_allrelev.loc[samples_allrelev.sum(axis=1)>0,]
	lefseinput = pd.concat([pd.DataFrame(metadata[["group_type"]].T)[samples], samples_allrelev]).reset_index().rename(columns = {"index":"subject_id"})
	lefseinput.to_csv(output_file, sep="\t", index=False)
    
	if os.path.exists(relevfile)==False:
		sys.exit(1)

if __name__ == '__main__':
    main()
