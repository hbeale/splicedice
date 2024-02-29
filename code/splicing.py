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
