import sys
fi = sys.argv[1]
green= sys.argv[2]


def dsize(ggdata):
    return {i.split()[0]: i.split() for i in open(ggdata+".len")}

genes = {}
for i in open(fi+".sorted.bam.merged"):
    subtype, start, end, count = i.split()
    start = int(start)
    end = int(end)
    count = int(count)
    try:
        genes[subtype]['count'] += count
        genes[subtype]['length'] += abs(end-start)
    except:
        genes[subtype] = { 
            'count': count,
            'length': abs(end-start),
        }

fo = open(fi+".sorted.bam.merged.quant", 'w')
gene_len = dsize(green)
# print gene_len
for gene in genes:
    # print gene, gene_len[gene]
    cov = genes[gene]['length']/float(gene_len[gene][1])
    fo.write("\t".join([
        gene,
        str(genes[gene]['count']),
        str(genes[gene]['length']),
        gene_len[gene][1],
        str(round(cov, 3))
    ])+"\n")
fo.close()





