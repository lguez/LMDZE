# This is a script in Bash.

# This script is useful only to build into a directory other than the
# source directory.

trap 'exit 1' ERR
set -x

dest_dir=/usr/local/guez/LMDZ/LMDZE_work/Compil_prod_g95

make -C $dest_dir -f $PWD/GNUmakefile $* srcdir=$PWD
