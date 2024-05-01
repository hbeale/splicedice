import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import stats
import numpy as np
import time
from tools import Table, Manifest

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
            group_indices = {v:[] for v in groups.values()}
            for i,name in enumerate(header):
                if name in groups:
                    group_indices[groups[name]].append(i)
            group_counts = []
            group_props = []
            ylabels = []
            for line in tsv:
                row = line.rstrip().split('\t')
                if "Tumor" in row[0]:
                    continue
                counts = []
                props = []
                for index in group_indices.values():
                    counts.append(len([i for i in index if float(row[i])<threshold]))
                    props.append(counts[-1]/len(index))
                ylabels.append(row[0])
                group_counts.append(counts)
                group_props.append(props)
        return header[1:],ylabels[::-1],group_indices,group_counts[::-1],group_props[::-1]

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

class ColorLoop:
    def __init__(self,colors=["darkblue","darkorange","green","gray","purple","black"]):
        self.colors = colors
        self.i = -1

    def next(self):
        self.i = (self.i + 1) % len(self.colors)
        return self.colors[self.i]


def read_betas(filename,interval_set):
    with open(filename) as tsv:
        header = tsv.readline().rstrip().split('\t')[1:]
        indices = {}
        for i,column in enumerate(header):
            stat,name = column.split("_",1)
            if stat == "median":
                indices[name] = [i]
            elif stat == "alpha":
                indices[name].append(i)
            elif stat == "beta":
                indices[name].append(i)
        betas = {}
        for line in tsv:
            interval,row = line.split("\t",1)
            if interval in interval_set:
                row = [float(x) for x in row.split('\t')]
                mabs = {}
                for name,index in indices.items():
                    mabs[name] = [row[i] for i in index]
                betas[interval] = mabs
        return betas


class PS_distribution:
    def __init__(self,interval,row,group_indices=None,betas={},width=0.05):
        self.width = width
        self.bins = np.arange(0,1+width,width)
        self.stack = [0 for bin in self.bins[:-1]]
        colors = ColorLoop()
        if group_indices:
            self.group_indices = group_indices
        else:
            self.group_indices = {"Values":[i for i in range(len(row))]}

        fw,fh = 6,3
        pw,ph = 4,2
        self.fig = plt.figure(figsize=(fw,fh))
        self.panel = self.fig.add_axes([0.5/fw,0.5/fh,pw/fw,ph/fh])
        self.labels = []
        for group,indices in self.group_indices.items():
            values = [row[i] for i in indices]
            self.add_hist(values,label=group,color=colors.next())
        for name,mab in betas.items():
            m,a,b = mab
            self.add_beta(a,b,label=name,color=colors.next())

        lw = 1 + (len(self.labels)//6)
        lh = min(6,1+(len(self.labels)%6)) / 2
        self.legend = self.fig.add_axes([4.6/fw,(3.5-lh)/fh,lw/fw,lh/fh])


    def add_hist(self,values,label,color,density=True):
        counts,bins = np.histogram(values,bins=self.bins,density=density)
        for i,floor in enumerate(self.stack):
            r = patches.Rectangle((bins[i],floor),self.width,counts[i],
                                  edgecolor="darkgray",facecolor=color)
            self.stack[i] += counts[i]
            self.panel.add_patch(r)
        self.panel.set_ylim(0,max(self.stack)*1.1)
        self.labels.append(["h",label,color])

    def add_beta(self,a,b,label,color):
        xs,ys = self.beta_points(a,b)
        self.panel.plot(xs,ys,color=color)
        self.labels.append(["b",label,color])

        
    def beta_points(self,a,b,xdist=0.005):
        xs,ys = [],[]
        for x in np.arange(xdist/2,1,xdist):
            xs.append(x)
            ys.append(stats.beta.pdf(x,a,b))
        return xs,ys
        
    def fill_legend(self):
        ys = [0,1,2,3,4,5]
        x = 0
        for which,label,color in self.labels:
            y = ys.pop()
            ys = [y] + ys
            if which == "h":
                r = patches.Rectangle((x+.1,y+.1),.8,.8,
                                  edgecolor="darkgray",facecolor=color)
                self.legend.add_patch(r)
            elif which == "b":
                self.legend.plot([0,1],[y+.5,y+.5],color=color)
            self.legend.text(x+.5,y+.5,label)
            

        self.legend.set_xticks([])
        self.legend.set_yticks([])  
        self.legend.set_xlim(0,x+1)
        if len(self.labels) > 5:
            self.legend.set_ylim(0,6)
        else:
            self.legend.set_ylim(y,6)         
        


    def save_fig(self,out_prefix,dpi=600):
        self.fill_legend()
        self.fig.savefig(f"{out_prefix}.{time.time()}.png",dpi=dpi)

class PCA_plot:
    def __init__(self,xs,ys,xy_pairs=[]):
        self.fig = plt.figure(figsize=(6,4))
        self.fig.scatter(xs,ys)

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--intervals",default=None,
                        help="List of intervals to plot (comma-separated).")
    parser.add_argument("-q","--query",default=None,
                        help="Table of p-values from signature query (pvals.tsv).")
    parser.add_argument("-m","--manifest",default=None,
                        help="Manifest of samples with groups for combining or labeling.")
    parser.add_argument("-p","--ps_table",default=None,
                        help="ps.tsv file with PS values for plotting.")
    parser.add_argument("-b","--beta",default=None,
                        help="beta.tsv file with parameters for fit beta distributions.")
    parser.add_argument("-o","--out_prefix",default="splicedice",
                        help="Output path and filename before extensions [Default: 'splicedice']")
    return parser.parse_args()




def main():
    args = get_args()
    if args.query and args.manifest:
        pmat = PvalMatrix(args.manifest,args.query)
        pmat.plot_table(args.out_prefix)

    if args.ps_table and args.intervals and args.beta and args.manifest:
        manifest = Manifest(args.manifest)
        intervals = set(args.intervals.split(","))
        betas = read_betas(args.beta,intervals)
        ps_table = Table(args.ps_table)
        group_indices = manifest.get_group_indices(ps_table.get_samples())
        for interval,row in ps_table.get_rows(interval_set=intervals):
            ps_plot = PS_distribution(interval,row,group_indices,betas[interval])
            ps_plot.save_fig(args.out_prefix)




if __name__=="__main__":
    main()

