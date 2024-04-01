
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

## 
from tools import Beta,Multi,Table

## Suppress Warnings
import warnings
warnings.simplefilter("ignore")

## Arguments and config parsing
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
    parser.add_argument("-b","--beta_file",default=None,
                        help="Filename and path for .beta.tsv file, previously output from fit_beta.")
    parser.add_argument("-a","--annotation",default=None,
                        help="GTF or splice_annotation.tsv file with gene annotation (optional for labeling/filtering)")
    parser.add_argument("-o","--output_prefix",
                        help = "Path and file prefix for the output file. '.sig.tsv' or '.match.tsv' will be appended to prefix.")
    parser.add_argument("-c","--config_file",default=None,
                        help = "Optional. For adjusting parameters of splicing analysis.")
    parser.add_argument("-ctrl","--control_name",default=None,
                        help="Sample group label that represents control for comparative analysis (default is first group in manifest).")
    parser.add_argument("-n","--n_threads",default=4,type=int,
                        help="Maximum number of processes to use at the same time.")
    parser.add_argument("-k","--threshold",default=0.05,type=float,
                       help="")
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

    manifest = Manifest(filename=args.manifest,n_threads=args.n_threads,threshold=args.threshold)

    ps_table = Table(filename=args.ps_table,store=None)

    if args.mode == "compare":
        print("Testing for differential splicing...")
        groups,med_stats,compare_stats = manifest.compare_multi(ps_table,threshold=args.threshold,delta_threshold=0.05)
        print("Writing...")
        manifest.write_sig(args.output_prefix,groups,med_stats,compare_stats)

    elif args.mode == "fit_beta":
        if args.sig_file:
            print("Reading...")
            groups,compare_stats = manifest.read_sig(args.sig_file)
            med_stats = None
        else:
            print("Testing for differential splicing...")
            groups,med_stats,compare_stats = manifest.compare_multi(ps_table,threshold=0.05,delta_threshold=0.05)

        print("Fitting beta distributions...")
        beta_stats = manifest.fit_betas(ps_table,compare_stats)
        print("Writing files...")
<<<<<<< HEAD
        if med_stats:
            manifest.write_sig(args.output_prefix,groups=groups,med_stats=med_stats,compare_stats=compare_stats)
        manifest.write_beta(args.output_prefix,groups=groups,beta_stats=beta_stats)
=======
        manifest.write_sig(args.output_prefix,groups=groups,med_stats=med_stats,compare_stats=compare_stats)
        manifest.write_beta(args.output_prefix,groups=manifest.groups,beta_stats=beta_stats)
>>>>>>> e71fc6bbf000f8f4856eef64637ffdab5a15cf37

    elif args.mode == "query":
        print("Reading...")
        groups,beta_stats = manifest.read_beta(args.beta_file)
        print("Querying...")
        samples,queries,pvals = manifest.query(ps_table,groups,beta_stats)
        print("Writing...")

        manifest.write_pvals(args.output_prefix,samples,queries,pvals)

#### Manifest Class ####    
class Manifest:
    def __init__(self,filename=None,control_name=None,n_threads=4,threshold=0.05):
        self.samples = []
        self.groups = {}
        self.get_group = {}
        self.beta = Beta()
        self.controls = {}
        self.n_threads = n_threads
        self.threshold = threshold
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
                        self.groups[group_name] = []
                    self.groups[group_name].append(sample_name)
            for group_name in self.groups.keys():
                self.controls[group_name] = control_name
            self.index = {sample:i for i,sample in enumerate(self.samples)}

    def get_group_indices(self,samples):
        group_indices = {k:[] for k in self.groups.keys()}
        for i,s in enumerate(samples):
            if s in self.get_group:
                group_indices[self.get_group[s]].append(i)
        return {k:v for k,v in group_indices.items() if len(v) > 0}
    

    def read_beta(self,beta_file):
        beta_stats = {}
        with open(beta_file) as tsv:
            header = tsv.readline().rstrip().split("\t")[1:]
            groups = {}
            for i,column in enumerate(header):
                info_type,group_name = column.split("_",1)
                if group_name not in groups:
                    groups[group_name] = {}
                groups[group_name][info_type] = i
            for line in tsv:
                row = line.rstrip('\n').split('\t')
                interval = row[0]
                beta_stats[interval] = []
                row = row[1:]
                for i,x in enumerate(row):
                    try:
                        row[i] = float(x)
                    except ValueError:
                        row[i] = float('nan')
                for group_name in groups:
                    median = row[groups[group_name]["median"]]
                    alpha = row[groups[group_name]["alpha"]]
                    beta = row[groups[group_name]["beta"]]
                    beta_stats[interval].append([median,alpha,beta])

        groups = list(groups.keys())
        return groups,beta_stats
    
    def read_sig(self,sig_file):
        compare_stats = {}
        with open(sig_file) as tsv:
            columns = tsv.readline().rstrip('\n').split("\t")[1:]
            groups = {}
            for i,column in enumerate(columns):
                info_type,group_name = column.split("_",1)
                if group_name not in groups:
                    groups[group_name] = {}
                groups[group_name][info_type] = i
            for line in tsv:
                row = line.rstrip('\n').split('\t')
                interval = row[0]
                compare_stats[interval] = []
                row = row[1:]
                for i,x in enumerate(row):
                    try:
                        row[i] = float(x)
                    except ValueError:
                        row[i] = float('nan')
                for group_name in groups:
                    delta = row[groups[group_name]["delta"]]
                    pval = row[groups[group_name]["pval"]]
                    compare_stats[interval].append([delta,pval])
        groups = list(groups.keys())
        return groups,compare_stats
          
    def write_sig(self,output_prefix,groups=None,med_stats=None,compare_stats=None):
        header = ["splice_interval"]        
        for name in groups:
            header.extend([f"median_{name}",f"mean_{name}",f"delta_{name}",f"pval_{name}"])
        with open(f"{output_prefix}.sig.tsv",'w') as tsv:
            tab = '\t'
            tsv.write(f"{tab.join(header)}\n")
            for interval in compare_stats.keys():
                row = [interval]
                for m,c in zip(med_stats[interval],compare_stats[interval]):
                    row.append(f'{m[0]}\t{m[1]}\t{c[0]}\t{c[1]}')
                tsv.write(f"{tab.join(row)}\n")

    def write_beta(self,output_prefix,groups=None,beta_stats=None,):
        header = ["splice_interval"]
        intervals = beta_stats.keys()
        for name in groups:
            header.extend([f"median_{name}",f"alpha_{name}",f"beta_{name}"])
        with open(f"{output_prefix}.beta.tsv",'w') as tsv:
            tab = '\t'
            tsv.write(f"{tab.join(header)}\n")
            for interval in intervals:
                mabs = []
                for mab in beta_stats[interval]:
                    mabs.extend(str(x) for x in mab)
                tsv.write(f"{interval}\t{tab.join(mabs)}\n")

    def write_pvals(self,output_prefix,samples,queries,pvals):
        with open(f"{output_prefix}.pvals.tsv",'w') as tsv:
            tab = "\t"
            tsv.write(f"query\t{tab.join(samples)}\n")
            for i in range(len(queries)):
                tsv.write(f"{queries[i]}\t{tab.join(str(x) for x in pvals[i])}\n")

    def compare(self,ps_table,threshold=1,delta_threshold=0):
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
                control_name = self.controls[group_name]
                if control_name:
                    delta = medians[group_name] - medians[control_name]                        
                    if len(group_values)>2 and len(values_by_group[control_name])>2:
                        D,pval = ranksums(group_values, values_by_group[control_name])

                        if pval < threshold and delta > delta_threshold:
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
    
    def compare_multi(self,ps_table,threshold=0.05,delta_threshold=0):

        med_stats = {}
        compare_stats = {}
        group_indices = self.get_group_indices(ps_table.get_samples())

        sizes = [f"{k} ({len(v)})" for k,v in group_indices.items()]
        print(f"Groups: {', '.join(sizes)}")
        import multiprocessing
        n = self.n_threads
        buffer_ratio = 10
        with multiprocessing.Manager() as manager:
            q1 = manager.Queue(maxsize = n * buffer_ratio)
            q2 = manager.Queue()
            o = manager.dict()
            read_process = multiprocessing.Process(target=Multi.mp_reader,
                                                   args=(ps_table.get_rows,None,q1,n))
            read_process.start()
            pool = [multiprocessing.Process(target=Multi.mp_do_rows,args=(q1,self.row_compare,group_indices,q2)) for n in range(n-2)] 
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
                interval,m_stats,c_stats = item
                for s in c_stats:
                    if s[0] and abs(s[0])>delta_threshold and s[1] and s[1] < threshold:
                        compare_stats[interval] = c_stats
                        med_stats[interval] = m_stats
                        break
            read_process.join()
            for p in pool:
                p.join()
        return list(group_indices.keys()),med_stats,compare_stats
    
    def row_compare(self,item,group_indices):
        interval,row = item
    
        nan_check = [np.isnan(x) for x in row]
        values_by_group = {g_name:[row[i] for i in index if not nan_check[i]] for g_name,index in group_indices.items()}
        medians = {g:np.median(values) for g,values in values_by_group.items()}
        c_stats = []
        m_stats = []
        for group_name,group_values in values_by_group.items():
            control_name = self.controls[group_name]
            if control_name:
                delta = medians[group_name] - medians[control_name]                        
                if len(group_values)>2 and len(values_by_group[control_name])>2:
                    D,pval = ranksums(group_values, values_by_group[control_name])
                else:
                    pval = None
            else:
                delta,pval = None,None
            c_stats.append([delta,pval])
            m_stats.append([medians[group_name],np.mean(group_values)])
        return interval,m_stats,c_stats
    
    def significant_intervals(self,compare_stats):
        significant = set()
        d_count = 0
        for interval,data in compare_stats.items():
            for d in data:
                if d[0] and abs(d[0]) > 0.1 and d[1] and d[1] < self.threshold :
                    significant.add(interval)
                    d_count += 1
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
        print("significant intervals:",len(interval_set))
        group_indices = self.get_group_indices(ps_table.get_samples())
        beta_stats = {}
        import multiprocessing
        n = self.n_threads
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
                beta_stats[item[0]] = item[1]
            read_process.join()
            for p in pool:
                p.join()
        return beta_stats
    
    def row_query_beta(self,row,beta_stats):
        interval,values = row
        probabilities = []
        for mab in beta_stats[interval]:
            if None in mab:
                probabilities.append([float('nan') for i in range(len(values))])
                continue
            m,a,b = mab
            sub_probs = []
            for x in values:
                sub_probs.append(self.beta.cdf(x,m,a,b))
            probabilities.append(sub_probs)
        return probabilities

    def query(self,ps_table,groups,beta_stats):
        import multiprocessing
        interval_set = set(beta_stats.keys())
        n = self.n_threads
        buffer_ratio = 8
        with multiprocessing.Manager() as manager:
            q1 = manager.Queue(maxsize = n * buffer_ratio)
            q2 = manager.Queue()
            samples = ps_table.get_samples()
            probs_by_sample = [[[] for j in range(len(samples))] for i in range(len(groups))]
            read_process = multiprocessing.Process(target=Multi.mp_reader,args=(ps_table.get_rows,interval_set,q1,n))
            read_process.start()
            pool = [multiprocessing.Process(target=Multi.mp_do_rows,args=(q1,self.row_query_beta,beta_stats,q2)) for n in range(n-2)] 
            for p in pool:
                p.start()
            done_count = 0
            loop_count = 0
            while True:
                item = q2.get()
                if item == "DONE":
                    done_count += 1
                    if done_count == n-2:
                        break
                    continue
                for i,sub_probs in enumerate(item):
                    for j,prob in enumerate(sub_probs):
                        probs_by_sample[i][j].append(prob)
                loop_count += 1
            read_process.join()
            for p in pool:
                p.join()
        pvals = []
        queries = []

        for i,group_of_probs in enumerate(probs_by_sample):
            for j in range(i+1,len(probs_by_sample)):
                comp_group_probs = probs_by_sample[j]
                queries.append(f"{groups[i]}_over_{groups[j]}")
                queries.append(f"{groups[j]}_over_{groups[i]}")
                first_pvals = []
                second_pvals = []
                for first,second in zip(group_of_probs,comp_group_probs):
                    s,pval = ranksums(first,second,alternative="greater",nan_policy="omit")
                    first_pvals.append(pval)
                    s,pval = ranksums(second,first,alternative="greater",nan_policy="omit")
                    second_pvals.append(pval)
                pvals.extend([first_pvals,second_pvals])

        return samples,queries,pvals
       
# Run main
if __name__ == "__main__":
    main()


