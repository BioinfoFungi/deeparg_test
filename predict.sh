diamond  blastx \
    -q ./output/clean/test.clean \
    -d /home/wangyang/workspace/gusphdproj-deeparg-ss-fbe063e24cf7/database/database/v2/features \
    -k 1000 \
    --id 80.0 \
    --sensitive \
    -e 1e-10 \
    -a ./tmp2/test.clean.deeparg.align


diamond view \
    -a tmp2/test.clean.deeparg.align.daa \
    -o tmp2/test.clean.deeparg.align.tsv
