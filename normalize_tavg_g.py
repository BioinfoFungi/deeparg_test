import sys
import numpy as np
fiArg = sys.argv[1]

def convertMarker2Class(rat_nreads):
    wnreads = sorted([(float(count)/(np.absolute(geneLen-algLen)+1),(np.absolute(geneLen-algLen)+1) ,count) for count,algLen,geneLen in rat_nreads], key=lambda x:x[0])
    den,num = zip(*[v[1:] for v in wnreads])
    sum_count = sum(num)
    loc_ab = float(sum(num))/float(sum(den)) if any(den) else 0.0


    # sum_count = sum([count for count, algLen, geneLen in rat_nreads])
    # sum_geneLen = sum([geneLen for count, algLen, geneLen in rat_nreads])
    # loc_ab = float(sum_count)/float(sum_geneLen) 

    return (sum_count,loc_ab)
rat_nreads = []
def aggregate(fiArg):
    for i in open(fiArg):
        subtype, gtype, count, algLen, geneLen, cov = i.split()
        if float(cov) <= 0.01:
            continue
        rat_nreads.append((subtype,gtype,int(count),int(algLen),int(geneLen)))
    return rat_nreads
aggregate(fiArg)

grouped_data = {}
for subtype, gtype, count, algLen, geneLen in rat_nreads:
    if gtype in grouped_data:
        grouped_data[gtype].append((count, algLen, geneLen))
    else:
        grouped_data[gtype] = [(count, algLen, geneLen)]

category_ab = []
for key, value in grouped_data.items():
    num,r = convertMarker2Class(value)
    category_ab.append((key, num,r))

sum_ab = sum([r for key, num,r in  category_ab])



res = [(key,num, r/sum_ab) for key, num,r in  category_ab]




fo2 = open(fiArg+'.subtype.ttavg_g', 'w')
fo2.write("#ARG-group\tReadCount\ttavg_g\n")
for key,num,r in res:
    print(key,"\t",num,"\t",r)
    fo2.write("\t".join([
            key,
            str(num),
            str(r)
        ])+"\n")

