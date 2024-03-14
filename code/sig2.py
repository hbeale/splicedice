
import numpy as np
from scipy.stats import ranksums

# Suppress Warnings
import warnings
warnings.simplefilter("ignore")

# Arguments and config parsing
from config import get_config

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",nargs="?",default="compare",choices=["compare","fit_beta","query"])
    parser.add_argument("-m","--manifest",default=None,
                        help="TSV file with list of samples (first column) and group labels (second column).")  
    parser.add_argument("-p","--ps_table",default=None,
                        help="Filename and path for .ps.tsv file, output from MESA.")
    parser.add_argument("-s","--sig_file",default=None,
                        help="Filename and path for .sig.tsv file, previously output from splicedice.")
    parser.add_argument("-a","--annotation",default=None,
                        help="GTF or splice_annotation.tsv file with gene annotation (optional for labeling/filtering)")
    parser.add_argument("-o","--output_prefix",
                        help = "Path and file prefix for the output file. '.sig.tsv' or '.match.tsv' will be appended to prefix.")
    parser.add_argument("-c","--config_file",default=None,
                        help = "Optional. For adjusting parameters of splicing analysis.")
    parser.add_argument("-ctrl","--control_name",default=None,
                        help="Sample group label that represents control for comparative analysis (default is first group in manifest).")
    parser.add_argument("-n","--n_threads",default=1,
                        help="Maximum number of processes to use at the same time.")
    return parser.parse_args()
     
def check_args_and_config(args,config):
    if args.mode == "compare":
        if not args.ps_table or not args.manifest:
            exit()
    elif args.mode == "fit_beta":
        return True
    elif args.mode == "compare":
        return True
    return True

#### Main ####
def main():
    # Arguments and configuration specifications
    args = get_args()
    config = get_config(args.config_file)
    check_args_and_config(args=args,config=config)

    print("Reading manifest")
    if args.manifest:
        manifest = Manifest(filename=args.manifest,control_name=args.control_name)
    else:
        manifest = Manifest()

    if args.sig_file:
        manifest.read_sig(args.sig_file)
        
    if args.ps_table:
        ps_table = Table(filename=args.ps_table,store=None)

    if args.mode == "compare":
        print("comparing...")
        manifest.compare(ps_table)
        print("writing...")
        manifest.write_sig(args.output_prefix)

    elif args.mode == "fit_beta":
        manifest.fit_betas(args)
        manifest.write_sig(args.output_prefix)

    elif args.mode == "query":
        manifest.query_and_write(ps_table)


#### Manifest Class ####    
class Manifest:
    def __init__(self,filename=None,control_name=None):
        self.groups = {}
        self.data = {}
        self.columns = []
        self.get_group = {}
        if filename:
            with open(filename) as manifest_file:
                self.groups = {}
                self.samples = []
                for line in manifest_file:
                    row = line.rstrip().split("\t")
                    sample_name,group_name = row[0:2]
                    self.samples.append(sample_name)
                    self.get_group[sample_name] = group_name
                    if group_name not in self.groups:
                        if not control_name:
                            control_name = group_name
                        self.groups[group_name] = Signature(name=group_name,manifest=self)
                    self.groups[group_name].add_sample(sample_name)

            self.sample_group_names = [name for name in self.groups.keys() if name != control_name]
            for group_name in self.sample_group_names:
                self.groups[group_name].add_control(control_name)
            self.index = {sample:i for i,sample in enumerate(self.samples)}

    def read_sig(self,sig_file):
        data = {}
        with open(sig_file) as tsv:
            columns = tsv.readline().rstrip().split("\t")[1:]
            for i,column in enumerate(columns):
                group_name = column.split("_")[-1]
                if group_name not in self.groups:
                    self.groups[group_name] = Signature(name=group_name,manifest=self)
                    self.groups[group_name].add_column(column,i)
            for line in tsv:
                row = line.rstrip().split()
                data[row[0]] = row[1:]
        self.data = data
        self.columns = columns

    def write_sig(self,output_prefix,intervals=None):
        header = ["splice_interval"]
        indices = []
        i = 0
        for name,signature in self.groups.items():
            if signature.control:
                header.extend([f"mean_{name}",f"median_{name}",f"delta_{name}",f"ranksums_{name}"])
                indices.extend([i+x for x in (0,1,2,3)])
            else:
                header.extend([f"mean_{name}",f"median_{name}"])
                indices.extend([i+x for x in (0,1)])
            i += 4
        if not intervals:
            intervals = self.data.keys()
        with open(f"{output_prefix}.sig.tsv",'w') as tsv:
            tab = '\t'
            tsv.write(f"{tab.join(header)}\n")
            for interval in intervals:
                tsv.write(f"{interval}\t{tab.join([str(self.data[interval][i]) for i in indices])}\n")

    def compare(self,ps_table):
        data = {}
        group_indices = self.get_group_indices(ps_table.get_samples())
        for interval,row in ps_table.get_rows():
            nan_check = [np.isnan(x) for x in row]
            group_values = {g:[row[i] for i in index if not nan_check[i]] for g,index in group_indices.items()}
            medians = {g:np.median(values) for g,values in group_values.items()}
            data[interval] = []
            for group_name,values in group_values.items():
                control_name = self.groups[group_name].control
                if control_name:
                    if len(values)>2 and len(group_values[control_name])>2:
                        D,pval = ranksums(values, group_values[control_name])
                    else:
                        pval = None
                    delta = medians[group_name] - medians[control_name]
                else:
                    pval,delta = None,None
                data[interval].extend([np.mean(values),medians[group_name],delta,pval])
        self.data = data
        return None
    
    def get_group_indices(self,samples):
        group_indices = {k:[] for k in self.groups.keys()}
        for i,s in enumerate(samples):
            if s in self.get_group:
                group_indices[self.get_group[s]].append(i)
        return group_indices

     #### Table class

#### Table Class ####    
class Table:
    def __init__(self,filename=None,samples=None,intervals=None,data=None,store=None):
        self.store = store
        if intervals and samples and data:
            self.samples = samples
            self.intervals = intervals
            self.data = data
        else:
            self.filename = filename
            self.samples = None
            self.intervals = None
            self.data = None

    def get_samples(self):
        if self.samples:
            return self.samples
        else:
            with open(self.filename) as tsv:
                return tsv.readline().rstrip().split('\t')[1:]
            
    def get_rows(self,interval_set=None):
        if self.store == None:
            with open(self.filename) as data_file:
                header = data_file.readline().rstrip().split('\t')[1:]
                if interval_set:
                    for line in data_file:
                        row = line.rstrip().split("\t")
                        if row[0] in interval_set:
                            yield (row[0],[float(x) for x in row[1:]])
                else:
                    for line in data_file:
                        row = line.rstrip().split("\t")
                        yield (row[0],[float(x) for x in row[1:]])
                        
#### Signature Class ####    
class Signature:
    def __init__(self,name=None,manifest=None,samples=[],control_group=None):
        self.name = name
        self.manifest = manifest
        self.samples = samples
        self.control = control_group
        
    def add_sample(self,sample):
        self.samples.append(sample)
        
    def add_control(self,control_group):
        self.control = control_group

    def add_column(self,label,i):
        self.columns[label] = i
        
    def sample_set(self):
        return set(self.samples)
        
# Run main
if __name__ == "__main__":
    main()


