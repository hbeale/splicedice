

import numpy as np
from scipy.stats import beta as stats_beta
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
             

class Signature:

    def __init__(self,ps_file=None,vs_file=None,sample_group=[],control_group=[]):

        # Lists of samples in each group
        self.groups = [sample_group,control_group]
        # Group comparison data
        if vs_file:
            self.vs_label, self.vs_data, = self.read_vsfile(vs_file)
        elif ps_file:
            self.vs_label, self.vs_data = self.compare(ps_file,sample_group,control_group)

        self.beta = None
        
        # Select events with a cutoff
        #i_label = {x:i for i,x in self.vs_label}
        #self.events = self.get_events(self.vs_data,i=i_label["pvalue"])

    def compare(self,ps_file,sample_group,control_g tyui yroup):
        with open(ps_file) as tsv:
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
                                     mean_sample,mean_control]
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
            data.extend([self.beta.alphas[0][interval],self.beta.betas[0][interval],
                         self.beta.alphas[1][interval],self.beta.betas[1][interval]])

    def query(self,ps_table):

        for interval,ps_vals in ps_table.get_rows():
            for val in ps_vals:
                
class Groups(dict):
    def __init__(self,groups):
        
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
            
            
        
        

class Betas:

    from scipy.stats import beta as stats_beta

    def __init__(self,ps_table,samples,intervals,exclude_zero_ones=False):

        self.ps_table = ps_table  
        self.samples = samples
        
        self.exclude = exclude_zero_ones
        if self.exclude:
            self.floc = 0
            self.fscale = 1
        else:
            self.floc = -0.001
            self.fscale = 1.002

        
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
                    a,b,l,s = stats_beta.fit(values,floc=self.floc,fscale=self.fscale)
                except:
                    a,b,l,s = None,None,None,None 
                alphas[group][interval] = a
                betas[group][interval] = b
                medians[group][interval] = median
        return alphas,betas,medians

    def probability(self,event,gp,x):
        
        
        



    
    def query(self,ps_values):
        
        pvals = {}    
        for interval,row in ps_table.get_rows():
            if interval not in self.intervals:
                continue
            for group in self.samples.groups:

                m = self.medians[group][event]
                if x == m:
                    p = 1
                else:
                    a = self.alphas[group][event]
                    b = self.betas[group][event]
                    if x > m:
                        p = 1 - stats_beta.cdf(val,a,b,loc=self.floc,scale=self.fscale)
                    else:
                        p = stats_beta.cdf(val,a,b,loc=self.floc,scale=self.fscale)
                P.append(p)
                
        
        return pvals
    


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
    def __init__(self,filename,read_function,row_function,collect_function,out_function,n_threads=3):
        if n_threads < 3:
            return None
        import multiprocessing
        self.n = n_threads
        buffer_ratio = 10
        with multiprocessing.Manager() as manager:
            q1 = manager.Queue(maxsize = self.n * buffer_ratio)
            q2 = manager.Queue()
            o = manager.list()
            self.read_and_do(read_function,filename,q1,row_function,q2,collect_function,o)
            self.out = out_function(o)
        return None
            
    def get(self):
        return self.out
            
    def reader(self,read_function,filename,q):
        for item in read_function(filename):
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

    def read_and_do(self,read_function,filename,q1,row_function,q2,collect_function,o):
        read_process = multiprocessing.Process(target=self.reader,args=(read_function,filename,q1))
        read_process.start()
        pool = [multiprocessing.Process(target=self.do_rows,args=(q1,row_function,q2)) for n in range(self.n-2)] 
        for p in pool:
            p.start()
        out_process = multiprocessing.Process(target=self.collect_rows,args=(q2,collect_function,o))
        for p in pool:
            p.join()
        read_process.join()
        q2.put("DONE")
        out_process.join()
        return None
    
#### Manifest parsing
def parse_manifest(manifest):
    if manifest:
        with open(manifest) as file:
            samples = []
            for line in file:
                row = line.rstrip().split("\t")
                sample = row[0]
                group = row[1]
                samples.append((sample,group))
    return samples

                
def parse_manifests(man1,man2):
    with open(man1) as file:
        signature_group = []
        for line in file:
            signature_group.append(line.rstrip().split("\t")[0])
    with open(man2) as file:
        control_group = []
        for line in file:
            control_group.append(line.rstrip().split("\t")[0])
    return signature_group,control_group


def get_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",nargs="?",default="compare",choices=["compare","fit_beta","query"])

    parser.add_argument("-p","--ps_table",default=None,
                        help="Filename and path for .ps.tsv file, output from MESA.")
    parser.add_argument("-v","--vs_table",default=None,
                        help="Filename and path for .sig.tsv file, previously output from mesa signature.")
    parser.add_argument("-a","--annotation_gtf",
                        help="GTF file with gene annotation (optional for labeling/filtering)")
    parser.add_argument("-m","--manifest",
                        help="TSV file with list of samples and labels.")
    parser.add_argument("-s","--sample_label",default="sample",type=str,
                        help="Label for samples in signature group (when using single manifest).")
    parser.add_argument("-c","--control_label",default="control",type=str,
                        help="Label for samples in control group (when using single manifest).")
    
    parser.add_argument("-m1","--manifest1",
                        help="File with list of samples for signature group.")
    parser.add_argument("-m2","--manifest2",
                        help="File with list of samples for control group.")
    
    parser.add_argument("-t","--threshold",default=1.0,type=float,
                        help="Significance threshold for inclusion in signature/output.")
    
    parser.add_argument("-ex","--exclude_zeros_ones",default=False,action="store_true")
    
    parser.add_argument("-o","--output_prefix",
                        help = "Path and file prefix for the output file. '.sig.tsv' or '.match.tsv' will be appended to prefix.")
    return parser
     
def main():

    from annotation import Annotation
    parser = get_parser()
    args = parser.parse_args()
    mode = args.mode
    if args.manifest:
        samples = parse_manifest(manifest=args.manifest)
        if args.sample_label:
            sample_group = [s[0] for s in samples if s[1]==args.sample_label]
            control_group = [s[0] for s in samples if s[1]==args.control_label]
    elif args.manifest1 and args.manifest2:
        sample_group,control_group = parse_manifests(man1=args.manifest1,
                                                 man2=args.manifest2)
    else:
        samples = None
    if args.vs_table:
        signature = Signature(vs_file=args.vs_table)
    elif args.ps_table:
        signature = Signature(ps_file=args.ps_table,
                          sample_group=sample_group,
                          control_group=control_group)
    else:
        print("PS table or previous comparison table is required. Exiting...")
        exit()
    if args.annotation_gtf:
        signature.add_annotation(Annotation(args.annotation_gtf))
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
    

    
                

    
    
        
            
        
