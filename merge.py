import sys

inputFile=sys.argv[1]
path_to_data=sys.argv[2]
genes = {}

def dsize(path_to_data):
    return {i.split()[0].split("|")[-1].upper(): i.split() for i in open('{}/database/v2/features.gene.length'.format(path_to_data))}

for i in open(inputFile+".merged"):
    subtype, start, end, count, Type = i.split()
    start = int(start)
    end = int(end)
    count = int(count)
    try:
        genes[subtype]['count'] += count
        genes[subtype]['length'] += abs(end-start)
        genes[subtype]['type'].append(Type)
        genes[subtype]['type'] = list(set(genes[subtype]['type']))
    except:
        genes[subtype] = {
            'count': count,
            'length': abs(end-start),
            'type': [Type]
        }
# print genes
fo = open(inputFile+".merged.quant", 'w')
gene_len = dsize(path_to_data)
# print gene_len
for gene in genes:
    # print gene, gene_len[gene]
    cov = genes[gene]['length']/float(gene_len[gene][1])
    fo.write("\t".join([
        gene,
        "/".join(genes[gene]['type']),
        str(genes[gene]['count']),
        str(genes[gene]['length']),
        gene_len[gene][1],
        str(round(cov, 3))
    ])+"\n")
fo.close()

