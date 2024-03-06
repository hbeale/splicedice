

## Python package imports
import numpy as np
from statsmodels.stats.multitest import multipletests
## splicedice code imports
from annotation import Annotation

#### ####               #### ####
     #### Table class

class Table:
    def __init__(self,filename=None,intervals=None,samples=None,data=None,store=None):

        if intervals and samples and data:
            self.samples = samples
            self.intervals = intervals
            self.data = data
        # elif filename and store:
        #     self.samples,self.intervals,self.data = load_from_file.read_table(filename)
        else:
            self.filename = filename
            self.samples = None
            self.intervals = None
            self.data = None

    # def load_from_file(self,filename)
    #        pass

    # def write_table(self,out_filename)
    #     with open(out_filename,'w') as tsv:
    #         tab = "\t"
    #         tsv.write(f"splice_interval\t{tab.join(self.samples.samples)}\n")
    #         for interval,row in zip(self.intervals.intervals,self.data):
    #             tsv.write(f"{interval}\t{tab.join(str(val) for val in row)}\n")
    #     return None
     
    # def get_sample(self,name)
    #     return [self.data[i][self.sample_index[name]] for i in range(self.intervals.n)]
    
    # def combine(self,other,keep="both"):
    #     if keep == "union":
    #         intervals = self.intervals.union(other.intervals)
    #     elif keep == "intersection":
    #         intervals = self.intervals.intersection(other.intervals)
    #     elif keep == "left":
    #         intervals = self.intervals
    #     elif keep == "right":
    #         intervals = other.intervals
    #     samples = self.samples.combine(other.samples)
    #     data = [[] for i in range(intervals.n)]
    #     for i,interval in enumerate(intervals):
    #         data[i] = self.get_row(interval) + other.get_row(interval)
    #     return Table(intervals=intervals,samples=samples,data=data)
    
    # def reset_na(self,splice_interval,counts)
    # name,left,right,strand = Intervals.parse_interval(splice_interval)
    # contig = (name,strand)
    # span_i = (left,right,self.intervals.index[splice_interval])
    # ex_count = np.zeros(range(self.samples.m))
    # for exclusion in self.exclusions[contig][span_i]:
    #     ex_count += counts[exclusion]
    #     splice_interval =  self.intervals.intervals[span_i[2]]

    def get_row(self,interval,get_exclusion=False):
        try:
            return self.data[self.intervals.index[interval]]
        except KeyError:
            if get_exclusion:
                row = ["nan" for i in range(self.samples.m)]
                # for other in self.exclusions[interval]:
                #     other>0             
                #return ...
                return row
            else:
                return ["nan" for i in range(self.samples.m)]
            
    def get_rows(self):
        if self.store == None:
            with open(self.filename) as data_file:
                header = data_file.readline().rstrip().split('\t')[1:]
                yield header[1:]
                for line in data_file:
                    row = line.rstrip().split("\t")
                    yield (row[0],[float(x) for x in row[1:]])

                    
    def get_rows_by_group(self):
        for i,interval in enumerate(self.intervals):
            for name,group in self.samples.groups.items():
                yield (name,interval,[self.data[i][j] for j in group])

#### ####               #### ####
    #### Samples class
                
class Samples(list):

    def __init__(self,manifest=None,control_name=None,man1=None,man2=None,samples=[],groups={}):
        self.control_name = control_name
        if manifest:
            samples,self.groups = self.parse_manifest(manifest)
        elif man1 and man2:
            samples,self.groups = self.parse_manifests(man1,man2)
            self.control_name = "2"
        elif groups or samples:
            self.groups = groups
        else:
            return None
        list.init(self,samples)
        self.index = {sample:i for i,sample in enumerate(self.samples)}

        self.sample_groups = {k:v for k,v in self.groups.items() if k!=self.control_name}

    def parse_manifest(self,manifest):
        with open(manifest) as file:
            samples = []
            groups = {}
            for i,line in enumerate(file):
                row = line.rstrip().split("\t")
                samples.append(row[0])
                if row[1] in groups:
                    groups[row[1]].append(i)
                else:
                    groups[row[1]] = [i]
                    if not self.control_name:
                        self.control_name = row[1]
        return samples,groups
    
    def parse_manifests(self,man1,man2):
        samples = []
        groups = {"1":[],"2":[]}
        with open(man1) as file:
            for i,line in enumerate(file):
                row = line.rstrip().split("\t")
                samples.append(row[0])
                groups["1"].append(i)
        with open(man2) as file:
            for j,line in enumerate(file):
                row = line.rstrip().split("\t")
                samples.append(row[0])
                groups["2"].append(i+1+j) #index in samples list
        return samples,groups
 
    # def combine(self,other):
    #     for sample in other.samples:
    #         if sample in self.index:
    #             raise ValueError('Overlapping sample names, cannot combine.')
    #     new_groups = self.groups.items() ^ other.groups.items()
    #     for group in self.groups.keys() & other.groups.keys():
    #         new_groups[group] = self.groups[group] + other.groups[group]
    #     new_samples = self.samples + other.samples
    #     return Samples(samples=new_samples,groups=new_groups)
    
#### ####               #### ####
    #### Intervals class

class Intervals(list):

    @staticmethod
    def parse_interval(string):
        name,span,strand = string.split(":")
        left,right = (int(s) for s in span.split("-"))
        return (name,left,right,strand)

    @staticmethod
    def get_string(name,left,right,strand):
        return f"{name}:{left}-{right}:{strand}"
    
    def __init__(self,intervals=[],contigs={}):
        if contigs:
            self.contigs = contigs
            if not intervals:
                intervals = self.get_intervals(contigs)
        list.__init__(self,intervals)
        # Store properties
        self.index = {interval:i for i,interval in enumerate(self.intervals)}
        self.n = len(self.intervals)
        self.exclusions = {}

    def get_contigs(self,intervals):
        contigs = {}
        for i,interval in enumerate(intervals):
            name,left,right,strand = self.parse_interval(interval)
            try:
                contigs[(name,strand)].append((left,right,i))
            except KeyError:
                contigs[(name,strand)] = [(left,right,i)]
        for span_list in contigs.values():
            span_list.sort()
        return contigs

    def get_intervals(self,contigs):
        intervals = []
        for contig,spans in contigs.items():
            name,strand = contig
            for left,right in spans:
                intervals.append(f"{name}:{left}-{right}:{strand}")
        return intervals


    def read_exclusion(self,exclusion_file):
        exclusions = {}
        with open(exclusion_file) as tsv:
            for line in tsv:
                interval,exclusive = line.rstrip().split("\t")
                exclusions[interval] = exclusive.split(",") # <<<--------- CHECK THIS CHARACTER
        return exclusions

    def find_exclusion(self,contigs):
        exclusions = {}
        for contig,span_list in contigs.items():
            span_list.sort()
            exclusions[contig] = {}
            active = []
            for new_span in span_list:
                left = new_span[0]
                exclusions[contig][new_span] = []
                new_active = []
                for span in active:
                    if left <= span[1]:
                        exclusions[contig][span].append(new_span)
                        exclusions[contig][new_span].append(span)
                        new_active.append(span)
                new_active.append(new_span)
                active = new_active
                exclusions[contig][new_span].append(new_span)
        return exclusions
    
    def combine(self,other,keep="union"):
        if keep == "union":
            contigs = self.contigs.items() ^ other.contigs.items()
            for contig in self.contigs.keys():
                if contig in other.contigs:
                    contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        elif keep == "intersection":
            contigs = {}
            for contig in self.contigs.keys():
                if contig in other.contigs:
                    contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)
    
    def intersection(self,other):
        contigs = {}
        for contig in self.contigs.keys():
            if contig in other.contigs:
                contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)
    
    def union(self,other):
        contigs = self.contigs.items() ^ other.contigs.items()
        for contig in self.contigs.keys():
            if contig in other.contigs:
                contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)

#### ####               #### ####
#### Signature class
from scipy.stats import ranksums

class Signature:

    def __init__(self,vs_file=None,ps_table=None):
        # Lists of samples in each group
        if vs_file:
            self.vs_label,self.vs_data = self.read_vsfile(vs_file)
        elif ps_table:
            self.ps_table = ps_table
            self.vs_label,self.vs_data = None,None
        self.beta = Beta()
        self.params = None

    def read_vsfile(self,vs_filename):
        vs_data = {}
        with open(vs_filename) as tsv:
            vs_label = tsv.readline().rstrip().split("\t")[1:]
            for line in tsv:
                row = line.rstrip().split()
                vs_data[row[0]] = row[1:]
        return vs_label,vs_data
        
    def write_vsfile(self,output_prefix):
        tab = "\t"
        with open (f"{output_prefix}.sig.tsv","w") as tsv:
            tsv.write(f"splice_interval\t{tab.join(self.vs_label)}\n")
            for interval,data in self.vs_data.items():
                tsv.write(f"{interval}\t{tab.join(str(x) for x in data)}\n")
    
    def compare(self,ps_table,samples=None):
        if not samples:
            samples = ps_table.samples
        vs_data = {}
        vs_label = [f"median_{samples.control_name}",f"mean_{samples.control_name}"]
        for group in samples.sample_groups.keys():
            vs_label.extend([f"median_{group}",f"mean_{group}",f"delta_med_{group}",f"ranksums_p_{group}"])
        for interval,row in ps_table.get_rows(dtype="arrays"):
            nan_check = [np.isnan(x) for x in row]
            control_values = [row[i] for i in samples.control_group if not nan_check[i]]
            if len(control_values) < 3:
                continue
            control_median = np.median(control_values)
            vs_data[interval] = [control_median,np.mean(control_values)]
            for group in samples.sample_groups.values():
                ps_values = [row[i] for i in group if not nan_check[i]]
                if len(ps_values) < 3:
                    vs_data[interval].extend([np.median(ps_values),np.mean(ps_values),None])
                    continue
                D,pval = ranksums(ps_values,control_values)
                median = np.median(ps_values)
                vs_data[interval].extend([median,np.mean(ps_values),median-control_median,pval])
        self.vs_label = vs_label
        self.vs_data = vs_data
        return None
    
    ## Beta fit helper functions for multiprocessing
    def collect_params(self,item,params):
        group,interval,mab = item
        try:
            params[group][interval] = mab
        except KeyError:
            params[group] = {interval:mab}
            
    def save_params(self,params):
        self.params = params
        return params

    def fit_betas(self):
        mr = MultiReader(self.ps_table.get_rows_by_group, 
                    self.beta.fit_beta,
                    self.collect_params,
                    self.save_params)
        
    # ## Set of querying functions for use by MultiReader
    def query_row(self,row):
        interval,values = row
        probabilities = {}
        for group in self.groups:
            probabilities[group] = []
            m,a,b = self.mab_params(group,interval)
            for x in values:
                probabilities[group].append(self.beta.cdf(x,m,a,b,self.beta.loc,self.beta.scale))
        return probabilities

    def gather_probabilities(self,item,out_list):
        for i,values in enumerate(item):
            for j,x in enumerate(values):
                out_list[i][j].append(x)

    def write_pvals(self,out_list):
        control_probabilities = out_list[0]
        for values in out_list[1:]:
            for probabilities in values:
                D,pval = ranksums(probabilities,control_probabilities)

    def fit_betas(self,ps_table):
        mr = MultiReader(ps_table.get_rows, 
                    self.query_row,
                    self.gather_probabilities,
                    self.write_pvals)

#### ####               #### ####
#### #### Beta distribution
from scipy.stats import beta as stats_beta

class Beta:

    @staticmethod
    def cdf(x,m,a,b,loc=-0.001,scale=1.002):
        if x == m:
            return 1
        else:
            if x > m:
                return 1 - stats_beta.cdf(x,a,b,loc=loc,scale=scale)
            else:
                return stats_beta.cdf(x,a,b,loc=loc,scale=scale)

    


    def __init__(self,ps_table,samples,intervals,exclude_zero_ones=False):

        self.ps_table = ps_table  
        self.samples = samples
        
        self.exclude = exclude_zero_ones
        if self.exclude:
            self.loc = 0
            self.scale = 1
        else:
            self.loc = -0.001
            self.scale = 1.002

        
        #self.alphas,self.betas,self.medians = self.fit_betas(self.sig.groups)
        


    def fit_beta(self,row):
        group,interval,values = row
        values = [x for x in values if not np.isna(x)]
        if not values:
            return (None,None,None)
        if self.exclude:
            values = [x for x in values if x != 0 and x != 1]
        try:
            a,b,l,s = stats_beta.fit(values,floc=self.loc,fscale=self.scale)
        except:
            a,b = None,None 
        return group,interval,(np.median(values),a,b)

    def query(self,group,event):
        return (self.medians[group][event],self.alphas[group][event],self.betas[group][event])
        
# functions for read_function,row_function,out_function
def yield_rows(filename):
    with open(filename) as tsv:
        for line in tsv:
            row = line.rstrip().split("\t")
            yield (row[0],row[1:])

def gz_yield_rows(filename):
    import gzip
    with gzip.open(filename) as tsv:
        for line in tsv:
            row = line.decode().rstrip().split("\t")
            yield (row[0],row[1:])

#### ####               #### ####
#### Multiprocessing
import multiprocessing

class MultiReader:
    def __init__(self,read_function,row_function,collect_function,out_function,n_threads=3):
        if n_threads < 3:
            return None
        import multiprocessing
        self.n = n_threads
        buffer_ratio = 10
        with multiprocessing.Manager() as manager:
            q1 = manager.Queue(maxsize = self.n * buffer_ratio)
            q2 = manager.Queue()
            o = manager.list()
            self.read_and_do(read_function,q1,row_function,q2,collect_function,o,out_function)
            self.out = out_function(o)
        return None
            
    def get(self):
        return self.out
            
    def reader(self,read_function,q):
        for item in read_function():
                q.put(item)
        for i in range(self.n):
            q.put("DONE")
        return None

    def do_rows(self,q,f,o):
        while True:
            item = q.get()
            if item == "DONE":
                break
            o.put(f(item))
        return None
    
    def collect_rows(self,q,f,o):
        while True:
            item = q.get()
            if item == "DONE":
                break
            f(item,o)
        return None

    def read_and_do(self,read_function,q1,row_function,q2,collect_function,o):
        read_process = multiprocessing.Process(target=self.reader,args=(read_function,q1))
        read_process.start()
        pool = [multiprocessing.Process(target=self.do_rows,args=(q1,row_function,q2)) for n in range(self.n-2)] 
        for p in pool:
            p.start()
        out_process = multiprocessing.Process(target=self.collect_rows,args=(q2,collect_function,o))
        read_process.join()
        for p in pool:
            p.join()
        q2.put("DONE")
        out_process.join()
        return None
    

#### Arg Parse
def get_parser():
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",nargs="?",default="compare",choices=["compare","fit_beta","query"])

    parser.add_argument("-p","--ps_table",default=None,
                        help="Filename and path for .ps.tsv file, output from MESA.")
    parser.add_argument("-v","--vs_table",default=None,
                        help="Filename and path for .sig.tsv file, previously output from splicedice.")
    
    parser.add_argument("-a","--annotation",default=None,
                        help="GTF or splice_annotation.tsv file with gene annotation (optional for labeling/filtering)")
    parser.add_argument("-m","--manifest",
                        help="TSV file with list of samples and labels.")
    parser.add_argument("-c","--control_label",default="control",type=str,
                        help="Label for samples in control group (when using single manifest).")
    parser.add_argument("-m1","--manifest1",
                        help="File with list of samples for signature group.")
    parser.add_argument("-m2","--manifest2",
                        help="File with list of samples for control group.")    
    parser.add_argument("-o","--output_prefix",
                        help = "Path and file prefix for the output file. '.sig.tsv' or '.match.tsv' will be appended to prefix.")
    parser.add_argument("-c","--config_file",default=None,
                        help = "Optional. For adjusting parameters of splicing analysis.")
    parser.add_argument("-ctrl","--control_name",default=None,
                        help="Sample group label that represents control for comparative analysis.")
    parser.add_argument("-n","--n_threads",default=1,
                        help="Maximum number of processes to use at the same time.")
    
    return parser
     
# Import config reader and default configs
from .config import get_config

# Function to check for appropriate input specs
def check_args_and_config(args,config):
    if args.mode == "compare":
        if not args.ps_table:
            exit()
        if not ((args.manifest and args.control) or (args.manifest1 and args.manifest2)):
            exit()
    elif args.mode == "fit_beta":
        return True
    elif args.mode == "compare":
        return True
    return True

# Main
def main():
    # Arguments and configuration specifications
    args = get_parser().parse_args()
    config = get_config(args.config_file)
    check_args_and_config(args=args,config=config)

    #### ####
    #### Initiating objects with files
    if args.manifest:
        samples = Samples(manifest=args.manifest,control=args.control_name)
    elif args.manifest1 and args.manifest2:
        samples = Samples(man1=args.manifest1, man2=args.manifest2)
    if args.ps_table:
        ps_table = Table(filename=ps_table,store=None)
    if args.vs_table:
        signature = Signature(vs_table=args.vs_table)
    else:
        signature = Signature(ps_table=ps_table,samples=samples)
    if args.annotation:
        annotation = Annotation(filename=args.annotation)

    # Execute the specified mode
    if args.mode == "compare":
        signature.compare(samples)
        signature.write_vsfile(args.output_prefix)

    elif args.mode == "fit_beta":
        signature.fit_betas(samples)
        signature.write_vsfile()

    elif args.mode == "query":
        signature.query_and_write(ps_table)

if __name__=="__main__":
    main()
    

    
                

    
    
        
            
        
