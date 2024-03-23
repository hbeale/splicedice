import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import time

class PvalMatrix:
    def __init__(self,manifest,pvals):
        groups = self.read_manifest(manifest)
        header,labels,group_counts = self.read_pvals(pvals,groups)

    def read_manifest(self,filename):
        with open(filename) as manifest:
            groups = {}
            for line in manifest:
                row = line.rstrip().split("\t")
                groups[row[0]] = row[1]
        return groups
    
    def read_pvals(self,pvals,groups,threshold=0.05):
        with open(pvals) as tsv:
            header = tsv.readline().split('\t')
            group_indices = {group:[] for group in groups.values()}
            for i,name in enumerate(header):
                if name in groups:
                    group_indices[groups[name]].append(i)
            group_counts = []
            labels = []
            for line in tsv:
                row = line.rstrip().split('\t')
                counts = []
                for index in group_indices.values():
                    counts.append(len([i for i in index if float(row[i])<threshold]))
                labels.append(row[0])
                group_counts.append(counts)
        return header[1:],labels,group_counts

class PS_distribution:
    def __init__(self,values):
        self.fig = plt.figure(figsize=(6,4))
        
        self.fig.hist(values)

    def add_beta(self,a,b):
        xs,ys = [],[]
        for x in np.arange(0.005,1,0.005):
            xs.append(x)
            ys.append(stats.beta.pdf(x,(a,b)))
        self.fig.plot(xs,ys)
        
    def save_fig(self,out_prefix,dpi=600):
        self.fig.save_fig(f"{out_prefix}.{time.time()}.png",dpi=dpi)

class PCA_plot:
    def __init__(self,xs,ys,xy_pairs=[]):
        self.fig = plt.figure(figsize=(6,4))
        self.fig.scatter(xs,ys)

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",
                        help="")
    parser.add_argument("-p","--plot_options",
                        help="")
    parser.add_argument("-o","--output_file_prefix",
                        help="")
    return parser.parse_args()


def main():
    args = get_args()



if __name__=="__main__":
    main()

