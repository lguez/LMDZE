# This is a script in Bash.

# This script produces vertical distributions for all sampling
# methods, for a given value of "llm". Run this script from the
# directory containing "test_disvert".

trap 'echo -e \\a; exit 1' ERR

read -p "llm = ? " llm
# (only used to name the output files in this script, not used by the
# Fortran program)

set -x

for s_sampling in param LMD5 strato1 strato2
  do
  test_disvert <<EOF
&disvert_nml S_SAMPLING="$s_sampling" /
EOF
  mv test_disvert.csv test_disvert_${llm}_$s_sampling.csv
  ln -sf test_disvert_${llm}_$s_sampling.csv \
      test_disvert_$s_sampling.csv
done
