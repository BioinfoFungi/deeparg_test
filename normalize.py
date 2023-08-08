import sys
import numpy as np


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
aggregate(sys.argv[2])

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



fiArg = sys.argv[2]
fo2 = open(fiArg+'.subtype.ttavg_g', 'w')
fo2.write("#ARG-group\tReadCount\ttavg_g\n")
for key,num,r in res:
    print(key,"\t",num,"\t",r)
    fo2.write("\t".join([
            key,
            str(num),
            str(r)
        ])+"\n")


print(rat_nreads)
def normalize(fi16s, fiArg, covi, parameters={}):
    N16s = sum([int(i.split()[1]) for i in open(fi16s) if float(i.split()[2]) >= 100])
    print( "Total number of 16S Reads in the sample: {}".format(N16s) )
    L16s = 1432
    Atype = {}
    Asubtype = {}

    # print N16s
    fo1 = open(fiArg+'.subtype', 'w')
    fo1.write("#ARG-group\tReadCount\t16s-NormalizedReadCount\n")
    for i in open(fiArg):
        subtype, gtype, count, algLen, geneLen, cov = i.split()
        if float(cov) <= covi:
            continue
        if N16s > 0:
            Asubtype[subtype] = (int(count)/float(geneLen))/(float(N16s)/L16s)
        else:
            Asubtype[subtype] = 0

        fo1.write("\t".join([
            subtype,
            str(count),
            str(Asubtype[subtype])
        ])+"\n")

   
        


        try:
            Atype[gtype][0] += int(count)
        except:
            Atype[gtype] = [int(count), geneLen]


    

    fo2 = open(fiArg+'.type', 'w')
    fo2.write("#ARG-category\tReadCount\t16s-NormalizedReadCount\n")
    Xtype = {}
    for itype in Atype:
        if N16s > 0:
            Xtype[itype] = (Atype[itype][0]/float(Atype[itype][1])) / \
                (float(N16s)/L16s)
        else:
            Xtype[itype] = 0
        fo2.write("\t".join([
            itype,
            str(Atype[itype][0]),
            str(Xtype[itype])
        ])+"\n")

    return True
normalize(
    sys.argv[1],
    sys.argv[2],
    float(1)/100,
    {'coverage': 1, 'identity_16s_alignment': 100.0}
)


# normalize(
#     self.sample_name + '.clean.sorted.bam.merged.quant',
#     self.sample_name + '.clean.deeparg.mapping.ARG.merged.quant',
#     float(self.data['parameters']['coverage'])/100,
#     self.data['parameters']
# )