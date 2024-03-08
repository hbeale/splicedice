import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import time

    
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

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input")
    parser.add_argument("-p","--plot_options")
    parser.add_argument("-o","--output_file_prefix")


    return parser.parse_args()


def main():
    pass

if __name__=="__main__":
    main()

