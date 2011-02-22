# This is a script in Bash.

# This script allows to build LMDZE into a directory other than
# "libf". If you just want to build into "libf" with the current
# compiler file, you can invoke "make" directly.

# Usage:
# make.sh [options and arguments for make]

set -xe

compiler=gfortran

##dest_dir=~/Bureau/Compil_prod_$compiler
dest_dir=/save/workdir_Lionel/LMDZ_work/LMDZE/Compil_prod_$compiler
# (Do not just use the name of the compiler as a directory name, it
# confuses some compilers.)

cp --update --verbose ../Compilers/$compiler.mk $dest_dir/compiler.mk
make -C $dest_dir -f $PWD/GNUmakefile $* libf_dir=$PWD
