from scipy.stats import beta as stats_beta
import numpy as np

#### Beta Class ####        
class Beta:
    def __init__(self,exclude_zero_ones=False):
        self.exclude = exclude_zero_ones
        if self.exclude:
            self.loc = 0
            self.scale = 1
        else:
            self.loc = -0.001
            self.scale = 1.002

    def cdf(self,x,m,a,b):
        if x == m:
            return 1
        else:
            if x > m:
                return 1 - stats_beta.cdf(x,a,b,loc=self.loc,scale=self.scale)
            else:
                return stats_beta.cdf(x,a,b,loc=self.loc,scale=self.scale)

    def fit_beta(self,values):
        values = [x for x in values if not np.isnan(x)]
        if not values:
            return (None,None,(None,None,None))
        if self.exclude:
            values = [x for x in values if x != 0 and x != 1]
        try:
            a,b,l,s = stats_beta.fit(values,floc=self.loc,fscale=self.scale)
        except:
            a,b = None,None 
        return (np.median(values),a,b)
    
#### Multi Class ####        
class Multi:

    @staticmethod
    def mp_reader(read_function,filter_function,q,n): 
        for item in read_function(filter_function):
            q.put(item)

        for i in range(n):
            q.put("DONE")
        return None

    @staticmethod
    def mp_do_rows(q,f,info,filter_function,o):
        while True:
            item = q.get()
            if item == "DONE":
                o.put("DONE")
                break
            o.put(f(item,info,filter_function))
        return None

#### Table Class ####        
class Table:
    def __init__(self,filename=None,samples=None,intervals=None,data=None,store=None):
        self.store = store
        if intervals and samples and data:
            self.samples = samples
            self.intervals = intervals
            self.data = data
        else:
            self.filename = filename
            self.samples = None
            self.intervals = None
            self.data = None

    def get_samples(self):
        if self.samples:
            return self.samples
        else:
            with open(self.filename) as tsv:
                return tsv.readline().rstrip().split('\t')[1:]
            
    def get_rows(self,filter_function=None):
        if self.store == None:
            with open(self.filename) as data_file:
                header = data_file.readline().rstrip().split('\t')[1:]
                if filter_function:
                    for line in data_file:
                        row = line.rstrip().split("\t")
                        if filter_function(row):
                            yield (row[0],[float(x) for x in row[1:]])
                else:
                    for line in data_file:
                        row = line.rstrip().split("\t")
                        yield (row[0],[float(x) for x in row[1:]])
                        