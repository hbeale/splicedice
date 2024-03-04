

## Python package imports
import numpy as np
from statsmodels.stats.multitest import multipletests
## splicedice code imports
from annotation import Annotation


class Table:
    def __init__(self,filename=None,intervals=None,samples=None,data=None,store=None):

        if intervals and samples and data:
            self.samples = samples
            self.intervals = intervals
            self.data = data
        elif filename and store:
            self.samples,self.intervals,self.data = self.read_table(filename)
        else:
            self.filename = filename
            self.samples = None
            self.intervals = None
            self.data = None

    def read_table(self,filename):
        with open(filename) as tsv:
            samples = tsv.readline().rstrip().split("\t")[1:]
            intervals = []
            data = []
            for line in tsv:
                row = line.rstrip().split("\t")
                intervals.append(row[0])
                data.append([float(x) for x in row[1:]])
            intervals = Intervals(intervals)
        return samples,intervals,data

    def write_table(self,out_filename):
        with open(out_filename,'w') as tsv:
            tab = "\t"
            tsv.write(f"splice_interval\t{tab.join(self.samples.samples)}\n")
            for interval,row in zip(self.intervals.intervals,self.data):
                tsv.write(f"{interval}\t{tab.join(str(val) for val in row)}\n")
        return None
    
    # def read_table_filter(self,table_file,intervals):
    #     interval_set = set(intervals)
    #     intervals = []
    #     with open(table_file) as table:
    #         samples = table.readline().rstrip().split("\t")[1:]
    #         data = []
    #         for line in table:
    #             row = line.rstrip().split("\t")
    #             if row[0] in interval_set:
    #                 intervals.append(row[0])
    #                 data.append([float(x) for x in row[1:]])
    #     return samples,intervals,data
    
    # def get_sample(self,name):
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
    
    # def reset_na(self,splice_interval,counts):
        
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
        # elif self.store == "byrow":
        #     for interval_row in zip(self.intervals,self.data):
        #         yield interval_row
        # elif self.store == "array":
        #     for i in range(self.m):
        #         yield self.intervals[i],self.data[i,:]
        # # elif other storage cases...
                    
    def get_rows_by_group(self):
        for i,interval in enumerate(self.intervals):
            for name,group in self.samples.groups.items():
                yield (name,interval,[self.data[i][j] for j in group])


#### ####
    #### Samples class
class Samples(list):

    def __init__(self,manifest=None,control_group=None,man1=None,man2=None,samples=[],groups={}):
        if manifest:
            samples,self.groups = self.parse_manifest(manifest)
            self.control_group = control_group
        elif man1 and man2:
            samples,self.groups = self.parse_manifests(man1,man2)
            self.control_group = "2"
        elif groups or samples:
            self.groups = groups
            self.control_group = control_group
        else:
            return None
        list.init(self,samples)
        self.index = {sample:i for i,sample in enumerate(self.samples)}

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
    ####

    
    # def combine(self,other):
    #     for sample in other.samples:
    #         if sample in self.index:
    #             raise ValueError('Overlapping sample names, cannot combine.')
    #     new_groups = self.groups.items() ^ other.groups.items()
    #     for group in self.groups.keys() & other.groups.keys():
    #         new_groups[group] = self.groups[group] + other.groups[group]
    #     new_samples = self.samples + other.samples
    #     return Samples(samples=new_samples,groups=new_groups)
    

            

#### #### 
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

#### ####
#### Signature class
from scipy.stats import ranksums
class Signature:

    def __init__(self,vs_file=None,ps_table=None):   #ps_file=None,sample_group=[],control_group=[]):

        # Lists of samples in each group
        # self.groups = [sample_group,control_group]
        # Group comparison data
        if vs_file:
            self.vs_label,self.vs_data, = self.read_vsfile(vs_file)
        elif ps_table:
            self.ps_table = ps_table
            self.vs_label,self.vs_data = None,None


#        self.beta = None
        
        # Select events with a cutoff
        #i_label = {x:i for i,x in self.vs_label}
        #self.events = self.get_events(self.vs_data,i=i_label["pvalue"])

    def compare(self,control=None):
        vs_data = {}
        vs_label = ["median1","median2","pvalue","mean1","mean2"]

        for interval,row in self.ps_table.get_rows_by_group():
            for group in groups:
                values = [x for x in ]



        """with open(ps_file) as tsv:
            names = tsv.readline().rstrip().split("\t")
            sample_indices = [i for i,x in enumerate(names) if x in sample_group]
            control_indices = [i for i,x in enumerate(names) if x in control_group]
            vs_data = {}
            vs_label = ["median1","median2","pvalue","mean1","mean2"]
            for line in tsv:
                row = line.rstrip().split("\t")
                interval = row[0]
                sig = [float(row[i]) for i in sample_indices if row[i] != "nan"]
                con = [float(row[i]) for i in control_indices if row[i] != "nan"]
                if len(sig) < 3 or len(con) < 3:
                    continue # default values????
                D,pval = ranksums(sig,con)
                median_sample = np.median(sig)
                median_control = np.median(con)
                mean_sample = np.mean(sig)
                mean_control = np.mean(con)
                vs_data[interval] = [median_sample,median_control,pval,
                                     mean_sample,mean_control]"""
        return vs_label,vs_data

    def add_annotation(self,annotation):
        if "gene" in self.vs_label:
            return None
        

    def read_vsfile(self,vs_filename):
        vs_data = {}
        with open(vs_filename) as tsv:
            vs_label = tsv.readline().rstrip().split("\t")[1:]
            for line in tsv:
                row = line.rstrip().split()
                vs_data[row[0]] = row[1:]
        return vs_label,vs_data
    
    def write_vsfile(self,output,threshold=0,i=2):
        tab = "\t"
        with open (output,"w") as tsv:
            tsv.write(f"interval\t{tab.join(self.vs_label)}\n")
            for interval,data in self.vs_data.items():
                if float(data[i]) <= threshold:
                    tsv.write(f"{interval}\t{tab.join(str(x) for x in data)}\n")
    
    def fit_beta(self,ps_filename):
        self.beta = Beta(self,ps_filename)
        self.vs_label.extend(["alpha1","beta1","alpha2","beta2"])
        for interval,data in self.vs_data.items():
            data([self.beta.alphas[0][interval],self.beta.betas[0][interval],
                         self.beta.alphas[1][interval],self.beta.betas[1][interval]])
    
    def mab_params(self,group,interval):
        return (self.beta.alphas[group][interval],
                self.beta.betas[group][interval])
   

    ## Set of functions for use by MultiReader
    def query_row(self,row):
        interval,values = row
        probabilities = {}
        for group in self.groups:
            probabilities[group] = []
            m,a,b = self.mab_params(group,interval)
            for x in values:
                probabilities[group].append(Betas.cdf(x,m,a,b,self.beta.loc,self.beta.scale))
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





class SignatureSet:

    def __init__(self,signatures=None,base_signature=None):
        if signatures and base_signature:
            self.signatures = [control] + signatures
        elif signatures:
            self.signatures = [None] + signatures
        elif base_signature:
            self.signatures = [base_signature]
        else:
            self.signatures = [None]


    def compare(self,ps_table,groups):
        vs_data = {}
        vs_label = []

        group_list = list(groups.keys())
        g0 = group_list[0]
        group_list = group_list[1:]

        vs_label = [f"median_{g0}",f"mean_{g0}"]
        for g in group_list:
            vs_label.extend((f"median_{g}",f"mean_{g}",f"ranksums_p_{g}"))


        for interval,row in ps_table.get_rows():
            nan_check = np.isnan(row)
            g0_values = [row[i] for i in groups[g0] if not nan_check[i]]
            if len(g0_values) < 3:
                continue
            vs_data[interval] = [np.median(g0_values),np.mean(g0_values)]
            for group in group_list:
                ps_values = [row[i] for i in positions if not nan_check[i]]
                if len(ps_values) < 3:
                    vs_data[interval].extend((None,None,None))
                    continue
                D,pval = ranksums(ps_values,g0_values)
                vs_data[interval].extend((np.median(ps_values),np.mean(ps_values),pval))
        return vs_label,vs_data


from scipy.stats import beta as stats_beta              
class Betas:

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

        
        self.alphas,self.betas,self.medians = self.fit_betas(self.sig.groups)
        
    def fit_betas(self,ps_table,samples,intervals):
        alphas = {group:{} for group in samples.groups}
        betas = {gp:{} for gp in (0,1)}
        medians = {gp:{} for gp in (0,1)}
        for interval,row_by_group in ps_table.get_rows_by_group(samples,intervals):
            for group,row in row_by_group:
                values = [float(ps) for ps in row if ps != "nan"]
                if len(values) == 0:
                            continue
                median = np.median(values)
                if self.exclude:
                    values = [x for x in values if x != 0 and x != 1]
                try:
                    a,b,l,s = stats_beta.fit(values,floc=self.loc,fscale=self.scale)
                except:
                    a,b,l,s = None,None,None,None 
                alphas[group][interval] = a
                betas[group][interval] = b
                medians[group][interval] = median
        return alphas,betas,medians

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

#### ####
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
    parser.add_argument("-a","--annotation",default=None
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
    parser.add_argument("-ctrl","--control",default=None,
                        help="Sample group label that represents control for comparative analysis.")
    parser.add_argument("-n","--n_threads",default=1,
                        help="Maximum number of processes to use at the same time.")
    
    return parser
     

from .config import get_config

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

def main():
    # Arguments and configuration specifications
    args = get_parser().parse_args()
    config = get_config(args.config_file)
    check_args_and_config(args=args,config=config)

    # Getting samples from manifest
    if args.manifest:
        samples = Samples(manifest=args.manifest,c_label=args.control)
    elif args.manifest1 and args.manifest2:
        samples = parse_manifests(man1=args.manifest1, man2=args.manifest2)

    # Loading annotation

    
    # Initiating objects with files
    if args.ps_table:
        ps_table = Table(filename=ps_table,store=None)
    if args.vs_table:
        signature = Signature(vs_table=args.vs_table)
    if args.annotation:
        annotation = Annotation(filename=args.annotation)


    
    if args.mode == "compare":
        signature = Signature(ps_table=ps_table,samples=samples)
        signature.compare()
        signature.write_vsfile()

    elif args.mode == "fit_beta":
        if args.vs_table:
            signature = Signature(ps_table=ps_table,vs_table=args.vs_table,samples=samples)
        else:
            signature = Signature(ps_table=ps_table,samples=samples)
            signature.compare()
        signature.fit_betas()
        signature.write_vsfile()

    elif args.mode == "query":
        signature = Signature(vs_table=args.vs_table)
        multireader = MultiReader(ps_table.get_rows,
                                  signature.query_row,
                                  signature.gather_probabilities,
                                  signature.write_pvals,
                                  args.n_threads)

        pass
    elif args.mode == None:
        pass



    # Loading vs_table or ps_table.
    if args.vs_table:
        signature = Signature(vs_file=args.vs_table)
    if args.ps_table:
        signature = Signature(ps_table=ps_table,samples=samples)
    else:
        print("PS table or previous comparison table is required. Exiting...")
        exit()

    if mode == "compare":
        signature.write_vsfile(args.output_prefix,threshold = args.threshold)
    elif mode == "fit_beta":
        signature.fit_beta(args.ps_table)
        signature.write_vsfile(args.output_prefix,threshold = args.threshold)
    elif mode == "query":
        if not signature.beta:
            print("VS Table does not contain beta distribution parameters. Run fit_beta method with original PS table. Exiting...")
            exit()
        signature.query(args.ps_table,samples)
        

if __name__=="__main__":
    main()
    

    
                

    
    
        
            
        
