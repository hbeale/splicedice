
class Annotation:
        
    def __init__(self,gtf_filename):
        self.annotated,self.gene_coords,self.transcript_ids = self.getAnnotated(gtf_filename)

    def getAnnotated(self,gtf_filename):
        gid_coords = {}
        gene_coords = {}
        genes = {}
        transcripts = {}

        with open(gtf_filename) as gtf:
            for line in gtf:

                if line.startswith("#"):
                    continue

                row = line.rstrip().split('\t')
                info = [x.split('"') for x in row[8].split(';')]
                chrom = row[0]
                strand = row[6]

                start = int(row[3])
                stop = int(row[4])-1

                if row[2] == "transcript":
                    tid =  [x[1] for x in info if 'transcript_id' in x[0]][0]
                    gene_name = [x[1] for x in info if 'gene_name' in x[0]][0]
                    genes[tid] = gene_name
                    transcripts[(tid,chrom,strand)] = []

                elif row[2] == "exon":
                    tid =  [x[1] for x in info if 'transcript_id' in x[0]][0]
                    transcripts[(tid,chrom,strand)].append((start,stop))

                elif row[2] == "gene":
                    gene_name = [x[1] for x in info if 'gene_name' in x[0]][0]
                    gid = [x[1] for x in info if 'gene_id' in x[0]][0]
                    try:
                        try:
                            gene_coords[(chrom,strand)] [(start,stop)].append(gene_name)
                        except KeyError:    
                            gene_coords[(chrom,strand)] [(start,stop)] = [gene_name]

                        try:
                            gid_coords[(chrom,strand)] [(start,stop) ].append(gid)
                        except KeyError:
                            gid_coords[(chrom,strand)] [(start,stop) ] = [gid]
                    except KeyError:
                        gene_coords[(chrom,strand)] = {(start,stop) : [gene_name]}
                        gid_coords[(chrom,strand)] = {(start,stop) : [gid]}              
                    
        annotated = {}
        transcript_ids = {}
        for transcript,exons in transcripts.items():
            tid,chromosome,strand = transcript
            if strand == "+":
                for i in range(len(exons)-1):
                    junction = (chromosome,exons[i][1],exons[i+1][0],strand)
                    if junction in annotated:
                        if genes[tid] not in annotated[junction]:
                            annotated[junction].append(genes[tid])
                            transcript_ids[junction].append(tid)
                    else:
                        annotated[junction] = [genes[tid]]
                        transcript_ids[junction] = [tid]
            elif strand == "-":
                for i in range(len(exons)-1):
                    junction = (chromosome,exons[i+1][1],exons[i][0],strand)
                    if junction in annotated:
                        if genes[tid] not in annotated[junction]:
                            annotated[junction].append(genes[tid])
                            transcript_ids[junction].append(tid)
                    else:
                        annotated[junction] = [genes[tid]]
                        transcript_ids[junction] = [tid]

        return annotated,gene_coords,transcript_ids
