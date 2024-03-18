# For analyzing alternative splicing events
# Written by Dennis Mulligan

class AlternativeSplicing:

    def alternative_sites(self,contigs):
        five_prime_sites = {}
        three_prime_sites = {}
        for contig,spans in contigs:
            if contig[1] == "+":
                five,three = 0,1
            elif contig[1] == "-":
                five,three = 1,0
            for splice_sites in spans:
                if (contig,splice_sites[five]) in five_prime_sites:
                    five_prime_sites[(contig,splice_sites[five])].append(splice_sites[three])
                else:
                    five_prime_sites[(contig,splice_sites[five])] = [splice_sites[three]]
                if (contig,splice_sites[three]) in three_prime_sites:
                    three_prime_sites[(contig,splice_sites[three])].append(splice_sites[five])
                else:
                    three_prime_sites[(contig,splice_sites[three])] = [splice_sites[five]]
        return five_prime_sites,three_prime_sites
    
    def exon_skipping(self,intervals,exons):
        exon_skipping = {}
        for contig,interval_spans in intervals.contigs.items():
            exon_skipping[contig] = []
            for left,right in interval_spans:
                for exon in exons.contigs[contig]:
                    if left < exon[0] and right > exon[1]:
                        exon_skipping[contig].append(((left,right),exon))
                    elif left > exon[1]:
                        break
        return exon_skipping
    
    def other_events(self,intervals):
        sandwich = {}
        overlap = {}


        for contig,interval_spans in intervals.contigs.items():
            sandwich[contig] = []
            overlap[contig] = []
        for contig,span_list in contigs.items():
            span_list.sort()
            active = []
            for left,right in span_list:
                exclusions[(contig,(left,right))] = []
                new_active = []
                for span in active:
                    if left <= span[1]:
                        exclusions[(contig,span)].append((left,right))
                        exclusions[(contig,(left,right))].append(span)
                        new_active.append(span)
                new_active.append((left,right))
                active = new_active
                exclusions[(contig,(left,right))].append((left,right))
        return exclusions
    

class Gene:
    def __init__(self,name,bounds=[],splice_intervals=[],transcripts=[]):
        self.name = name
        self.bounds = bounds
        self.splice_intervals = splice_intervals
        self.transcripts = transcripts

    #### Intervals class

class Intervals(list):

    @staticmethod
    def parse_interval(string):
        name,span,strand = string.split(":")
        left,right = (int(s) for s in span.split("-"))
        return (name,left,right,strand)

    @staticmethod
    def get_string(name,left,right,strand):
        return f"{name}:{left}-{right}:{strand}"
    
    def __init__(self,intervals=[],contigs={}):
        if contigs:
            self.contigs = contigs
            if not intervals:
                intervals = self.get_intervals(contigs)
        list.__init__(self,intervals)
        # Store properties
        self.index = {interval:i for i,interval in enumerate(self.intervals)}
        self.n = len(self.intervals)
        self.exclusions = {}

    def get_contigs(self,intervals):
        contigs = {}
        for i,interval in enumerate(intervals):
            name,left,right,strand = self.parse_interval(interval)
            try:
                contigs[(name,strand)].append((left,right,i))
            except KeyError:
                contigs[(name,strand)] = [(left,right,i)]
        for span_list in contigs.values():
            span_list.sort()
        return contigs

    def get_intervals(self,contigs):
        intervals = []
        for contig,spans in contigs.items():
            name,strand = contig
            for left,right in spans:
                intervals.append(f"{name}:{left}-{right}:{strand}")
        return intervals


    def read_exclusion(self,exclusion_file):
        exclusions = {}
        with open(exclusion_file) as tsv:
            for line in tsv:
                interval,exclusive = line.rstrip().split("\t")
                exclusions[interval] = exclusive.split(",") # <<<--------- CHECK THIS CHARACTER
        return exclusions

    def find_exclusion(self,contigs):
        exclusions = {}
        for contig,span_list in contigs.items():
            span_list.sort()
            exclusions[contig] = {}
            active = []
            for new_span in span_list:
                left = new_span[0]
                exclusions[contig][new_span] = []
                new_active = []
                for span in active:
                    if left <= span[1]:
                        exclusions[contig][span].append(new_span)
                        exclusions[contig][new_span].append(span)
                        new_active.append(span)
                new_active.append(new_span)
                active = new_active
                exclusions[contig][new_span].append(new_span)
        return exclusions
    
    def combine(self,other,keep="union"):
        if keep == "union":
            contigs = self.contigs.items() ^ other.contigs.items()
            for contig in self.contigs.keys():
                if contig in other.contigs:
                    contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        elif keep == "intersection":
            contigs = {}
            for contig in self.contigs.keys():
                if contig in other.contigs:
                    contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)
    
    def intersection(self,other):
        contigs = {}
        for contig in self.contigs.keys():
            if contig in other.contigs:
                contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)
    
    def union(self,other):
        contigs = self.contigs.items() ^ other.contigs.items()
        for contig in self.contigs.keys():
            if contig in other.contigs:
                contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)

#### ####               #### ####
