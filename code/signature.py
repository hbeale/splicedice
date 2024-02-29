

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

    def compare(self,ps_file,sample_group,control_group):
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
        

class Beta:

    from scipy.stats import beta as stats_beta

    def __init__(self,signature,ps_filename,exclude_zero_ones=False):

        self.sig = signature
        self.ps_filename = ps_filename   
        self.exclude = exclude_zero_ones
        if self.exclude:
            self.floc = 0
            self.fscale = 1
        else:
            self.floc = -0.001
            self.fscale = 1.002
        self.alphas,self.betas,self.medians = self.fit_betas(ps_filename,self.sig.groups)
        
    def fit_betas(self,ps_filename,groups):
        alphas = {gp:{} for gp in (0,1)}
        betas = {gp:{} for gp in (0,1)}
        medians = {gp:{} for gp in (0,1)}
        with open(ps_filename) as tsv:
            header = tsv.readline().rstrip().split("\t")
            pos = {}
            for gp in (0,1):
                pos[gp] = [i for i,name in enumerate(header) if name in groups[gp]]
            for line in tsv:
                row = line.rstrip().split("\t")
                interval = row[0]
                if True:
                    for gp in (0,1):
                        values = [float(row[i]) for i in pos[gp] if row[i] != "nan"]
                        if len(values) == 0:
                            continue
                        median = np.median(values)
                        if self.exclude:
                            values = [x for x in values if x != 0 and x != 1]
                        try:
                            a,b,l,s = stats_beta.fit(values,floc=self.floc,fscale=self.fscale)
                        except:
                            a,b,l,s = None,None,None,None 
                        alphas[gp][interval] = a
                        betas[gp][interval] = b
                        medians[gp][interval] = median
        return alphas,betas,medians

    def probability(self,event,gp,x):
        a = self.alphas[gp][event]
        b = self.betas[gp][event]
        med = self.medians[gp][event]
        
        if x == med:
            p = 1
        elif x > med:
            p = 1 - stats_beta.cdf(x,a,b,loc=self.floc,scale=self.fscale)
        elif x < med:
            p = stats_beta.cdf(x,a,b,loc=self.floc,scale=self.fscale)
        return p


    def get_pval(self,events,values):

        probabilities = {gp:[] for gp in (0,1)}
        for event,value in zip(events,values):
            for gp in (0,1):
                probabilities[gp].append(self.probability(event,gp,value))
                
        pval = ranksums(probabilities[0],probabilities[1],alternative="greater").pvalue
        return pval

    def get_all_pvals(self,ps_table):
        """
        zip(ps_table.values,ps_table.names)
                    
        pvals = {}
        for name,values in zip(header[1:],ps_values):
            pvals[name] = self.get_pval(events,values)


        from multiprocessing import Pool
    
        n_threads = 8
        with Pool(n_threads) as pool:
            for result in pool.imap(self.get_pval,ps_values):
                name,pval = result
                pvals[name] = pval
        """

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
    

    
                

    
    
        
            
        
