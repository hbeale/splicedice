from scipy.stats import beta as stats_beta
import numpy as np
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
    

class Multi:

    @staticmethod
    def mp_reader(read_function,i,q,n): 
        for i,item in enumerate(read_function(i)):
                q.put(item)
                if i == 1000:
                    break
        for i in range(n):
            q.put("DONE")
        return None

    @staticmethod
    def mp_do_rows(q,f,o):
        while True:
            item = q.get()
            if item == "DONE":
                break
            o.put(f(item))
        return None