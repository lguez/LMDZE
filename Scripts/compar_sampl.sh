# This is a script in Bash.

# This script produces vertical distributions for all sampling
# methods, for a given value of "llm".

trap 'echo -e \\a; exit 1' ERR

read -p "llm = ? " llm
# (only used to name the output files in this script, not used by the
# Fortran program)

set -x

for sigma_sampling in param LMD5 strato1 strato2
  do
  Compil_prod/test_disvert <<EOF
&disvert_nml SIGMA_SAMPLING="$sigma_sampling" /
EOF
  mv test_disvert.csv test_disvert_${llm}_$sigma_sampling.csv
  ln -sf test_disvert_${llm}_$sigma_sampling.csv \
      test_disvert_$sigma_sampling.csv
done
