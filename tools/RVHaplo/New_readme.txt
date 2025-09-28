For PacBio HiFi (recommended Racon path):

bash rvhaplo.sh \
  -i sample.hdv.sam \
  -r data/hdv_ref.fa \
  -o results \
  -p sample \
  -t 16 \
  --polisher racon


For ONT (if you have medaka=1.6.0 available):

--polisher medaka