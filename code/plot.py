import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import stats
import numpy as np
import time


def get_color(x):
    if x < 0.5:
        return (1-2*x,1-2*x,1)
    else:
        return (0,0,1.5-x)
        
def get_text_color(x):
    if x>0.5:
        return 'white'
    else:
                return 'black'

class PvalMatrix:
    def __init__(self,manifest,pvals):
        groups = self.read_manifest(manifest)
        self.header,self.ylabels,self.group_indices,self.group_counts,self.group_props = self.read_pvals(pvals,groups)
        self.xlabels = []
        for label,index in self.group_indices.items():
            a=5
            if len(label) > a+1:
                new = label[:a]
                for i in range(a,len(label),a):
                    new = new+'\n'+label[i:i+a]
                label = new
            self.xlabels.append(f"{label}\n(n={len(index)})")

    def read_manifest(self,filename):
        with open(filename) as manifest:
            groups = {}
            for line in manifest:
                row = line.rstrip().split("\t")
                if row[1] == "Tumor":
                    continue
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
            group_props = []
            labels = []
            for line in tsv:
                row = line.rstrip().split('\t')
                if "Tumor" in row[0]:
                    continue
                counts = []
                props = []
                for index in group_indices.values():
                    counts.append(len([i for i in index if float(row[i])<threshold]))
                    props.append(counts[-1]/len(index))
                labels.append(row[0])
                group_counts.append(counts)
                group_props.append(props)
        return header[1:],labels,group_indices,group_counts,group_props

    def plot_table(self,filename,color=get_color,text=get_text_color):
        pw = .6 * len(self.xlabels)
        ph = .32 * len(self.ylabels)
        w = pw + 2.5
        h = ph + 1
        plt.figure(figsize=(w,h))
        panel = plt.axes([1.75/w,.9/h,pw/w,ph/h])
        colorbar = plt.axes([(1.85+pw)/w,.9/h,.15/w,.5/h]) 
        panel.set_yticks(range(1,len(self.ylabels)+1),self.ylabels)
        panel.set_xticks(range(1,len(self.xlabels)+1),self.xlabels)   
        colorbar.tick_params(left=False,right=True,bottom=False,top=False,
                             labelleft=False,labelright=True,labelbottom=False,labeltop=False)
        colorbar.set_yticks([0,0.5,1],['0%','50%','100%'])
        for i in range(100):
            p = i/100
            rectangle = patches.Rectangle([0,p],1,.01,lw=0,fc=color(p))
            colorbar.add_patch(rectangle)
        for i in range(len(self.ylabels)):
            for j in range(len(self.xlabels)):
                p = self.group_props[i][j]
                rectangle = patches.Rectangle([j+.5,i+.5],1,1,lw=0.5,ec='black',fc=color(p))
                panel.add_patch(rectangle)
                panel.text(j+1,i+1,self.group_counts[i][j],color=text(p),va='center',ha='center')
        panel.set_xlim(0.5,len(self.xlabels)+.5)
        panel.set_ylim(0.5,len(self.ylabels)+.5)
        plt.savefig(f"{filename}_match.png")


class PS_distribution:
    def __init__(self,values,label):
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
    parser.add_argument("-q","--query",default=None,
                        help="Table of p-values from signature query (pvals.tsv).")
    parser.add_argument("-m","--manifest",default=None,
                        help="Manifest of samples with groups for combining or labeling.")
    parser.add_argument("-o","--output_file_prefix",default="splicedice",
                        help="Output path and filename before extensions [Default: 'splicedice']")
    return parser.parse_args()




def main():
    args = get_args()
    if args.query and args.manifest:
        pmat = PvalMatrix(args.query,args.manifest)
        pmat.plot_table(args.out_prefix)


if __name__=="__main__":
    main()

