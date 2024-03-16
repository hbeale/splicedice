
import numpy as np
from scipy.stats import ranksums

from tools import Beta,Multi,Table

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

    if args.manifest:
        manifest = Manifest(filename=args.manifest)
    else:
        manifest = Manifest()

    if args.ps_table:
        ps_table = Table(filename=args.ps_table,store=None)

    if args.mode == "compare":
        groups,med_stats,compare_stats = manifest.compare(ps_table)
        manifest.write_sig(args.output_prefix,groups,med_stats,compare_stats)


    elif args.mode == "fit_beta":
        if args.sig_file:
            compare_stats = manifest.read_sig(args.sig_file)
            beta_stats = manifest.fit_betas(ps_table,compare_stats)
        else:
            compare_stats = manifest.compare(ps_table,only_significant=True)
            beta_stats = manifest.fit_betas(ps_table,compare_stats)
        manifest.write_sig(args.output_prefix, beta_stats)

    elif args.mode == "query":
        manifest.read_sig(args.sig_file)
        beta_stats = manifest.get_beta_stats()
        manifest.query(ps_table,beta_stats)

#### Manifest Class ####    
class Manifest:
    def __init__(self,filename=None,control_name=None):
        self.samples = []
        self.groups = {}
        self.data = {}
        self.columns = []
        self.get_group = {}
        self.beta = Beta()
        if filename:
            with open(filename) as manifest_file:
                for line in manifest_file:
                    row = line.rstrip().split("\t")
                    sample_name,group_name = row[0:2]
                    self.samples.append(sample_name)
                    self.get_group[sample_name] = group_name
                    if group_name not in self.groups:
                        if not control_name:
                            control_name = group_name
                        self.groups[group_name] = Signature(label=group_name,manifest=self)
                    self.groups[group_name].add_sample(sample_name)

            self.sample_group_names = [name for name in self.groups.keys() if name != control_name]
            for group_name in self.sample_group_names:
                self.groups[group_name].add_control(control_name)
            self.index = {sample:i for i,sample in enumerate(self.samples)}

    def get_group_indices(self,samples):
        group_indices = {k:[] for k in self.groups.keys()}
        for i,s in enumerate(samples):
            if s in self.get_group:
                group_indices[self.get_group[s]].append(i)
        return {k:v for k,v in group_indices.items() if len(v) > 0}
    
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
        return columns,data

    def write_sig(self,output_prefix,intervals=None):
        header = ["splice_interval"]        
        for name,signature in self.groups.items():
            if 'compare' in self.data:
                if signature.control:
                    header.extend([f"mean_{name}",f"median_{name}",f"delta_{name}",f"ranksums_{name}"])
                else:
                    header.extend([f"mean_{name}",f"median_{name}"])
                if 'fit' in self.data:
                    header.extend([f"alpha_{name}",f"beta_{name}"])
            elif 'fit' in self.data:
                header.extend([f"median_{name}",f"alpha_{name}",f"beta_{name}"])
        if 'fit' in self.data:
            intervals = self.data['fit'].keys()
        elif 'compare' in self.data:
            intervals = self.data['compare'].keys()
        with open(f"{output_prefix}.sig.tsv",'w') as tsv:
            tab = '\t'
            tsv.write(f"{tab.join(header)}\n")
            for interval in intervals:
                row = []
                if 'compare' in self.data and 'fit' in self.data:
                    for c,f in zip(self.data['compare'][interval],self.data['fit'][interval]):
                        row.extend([str(x) for x in c]+[str(x) for x in f[1:]])
                elif 'compare' in self.data:
                    for values in self.data['compare'][interval]:
                        row.extend([str(x) for x in values])
                elif 'fit' in self.data:
                    for values in self.data['fit'][interval]:
                        row.extend([str(x) for x in values])
                tsv.write(f"{interval}\t{tab.join(row)}\n")

    def compare(self,ps_table,threshold=1):
        med_stats = {}
        compare_stats = {}
        group_indices = self.get_group_indices(ps_table.get_samples())
        for interval,row in ps_table.get_rows():
            nan_check = [np.isnan(x) for x in row]
            values_by_group = {g_name:[row[i] for i in index if not nan_check[i]] for g_name,index in group_indices.items()}
            medians = {g:np.median(values) for g,values in values_by_group.items()}
            stats = [] 
            to_add = False
            for group_name,group_values in values_by_group.items():
                stats.append([medians[group_name],np.mean(group_values)])
                control_name = self.groups[group_name].control
                if control_name:
                    delta = medians[group_name] - medians[control_name]
                    if len(group_values)>2 and len(values_by_group[control_name])>2:
                        D,pval = ranksums(group_values, values_by_group[control_name])
                        if pval < threshold:
                            to_add = True
                    else:
                        pval = None
                else:
                    delta,pval = None,None
                stats.append([delta,pval])
            if to_add:
                compare_stats[interval] = stats
                med_stats[interval] = [[medians[group_name],np.mean(group_values)] for group_name,group_values in values_by_group.items()]
        return list(group_indices.keys()),med_stats,compare_stats
    
    def significant_intervals(self,compare_stats,threshold=0.05):
        significant = set()
        for interval,data in compare_stats.items():
            for d in data:
                if d[1] and d[1] < threshold:
                    significant.add(interval)
                    break
        return significant
    
    def row_fit_beta(self,row,group_indices):
        interval,row = row
        nan_check = [np.isnan(x) for x in row]
        mab_row = []
        for index in group_indices.values():
            mab_row.append(self.beta.fit_beta([row[i] for i in index if not nan_check[i]]))
        return interval,mab_row

    def fit_betas(self,ps_table,compare_stats):
        interval_set = self.significant_intervals(compare_stats)
        group_indices = self.get_group_indices(ps_table.get_samples())
        data = {}
        import multiprocessing
        n = 8
        buffer_ratio = 10
        with multiprocessing.Manager() as manager:
            q1 = manager.Queue(maxsize = n * buffer_ratio)
            q2 = manager.Queue()
            o = manager.dict()
            read_process = multiprocessing.Process(target=Multi.mp_reader,
                                                   args=(ps_table.get_rows,interval_set,q1,n))
            read_process.start()
            pool = [multiprocessing.Process(target=Multi.mp_do_rows,args=(q1,self.row_fit_beta,group_indices,q2)) for n in range(n-2)] 
            for p in pool:
                p.start()
            done_count = 0
            while True:
                item = q2.get()
                if item == "DONE":
                    done_count += 1
                    if done_count == n-2:
                        break
                    continue
                data[item[0]] = item[1]
            read_process.join()
            for p in pool:
                p.join()
        return data
    
    def get_beta_stats(self):

        pass

    def row_query_beta(self,row,beta_stats):
        interval,values = row
        probabilities = []
        for m,a,b in beta_stats[interval]:
            for x in values:
                probabilities.append(self.beta.cdf(x,m,a,b))
        return probabilities

    def query(self,ps_table,beta_stats):
        import multiprocessing
        interval_set = set(beta_stats.keys())
        n = 6
        buffer_ratio = 10
        with multiprocessing.Manager() as manager:
            q1 = manager.Queue(maxsize = n * buffer_ratio)
            q2 = manager.Queue()
            o = manager.dict()
            temp_interval,params = beta_stats.popitem()
            samples = ps_table.get_samples()
            probs_by_sample = [[[] for j in range(len(samples))] for i in range(len(params))]
            beta_stats[temp_interval] = params
            read_process = multiprocessing.Process(target=Multi.mp_reader,args=(ps_table.get_rows,interval_set,q1,n))
            read_process.start()
            pool = [multiprocessing.Process(target=Multi.mp_do_rows,args=(q1,self.row_query_beta,beta_stats,q2)) for n in range(n-2)] 
            for p in pool:
                p.start()
            done_count = 0
            while True:
                item = q2.get()
                if item == "DONE":
                    done_count += 1
                    if done_count == n-2:
                        break
                    continue
                for i,g_probs in enumerate(item):
                    for j,prob in enumerate(g_probs):
                        probs_by_sample[i][j].append(prob)
            read_process.join()
            for p in pool:
                p.join()
        a = len(probs_by_sample)
        pvals = [[1 for i in range(a)] for j in range(a)]
        labels = [[1 for i in range(a)] for j in range(a)]
        for i,sample_probs in enumerate(probs_by_sample):
            for j,comp_probs in enumerate(probs_by_sample):
                labels[i][j] = f"{samples[i]}_over_{samples[j]}"
                if i==j:
                    continue
                s,pval = ranksums(sample_probs,comp_probs,alternative="greater",nan_policy="omit")
                pvals[i][j] = pval
        return labels,pvals
    

#### Signature Class ####    
class Signature:
    def __init__(self,label=None,manifest=None,samples=[],control_group=None):
        
        # Saving or instantiating sample group information
        self.label = label
        self.manifest = manifest
        self.samples = samples
        self.control = control_group

        # Instantiating data
        self.data = []
        self.labels = {}

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


