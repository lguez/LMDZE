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

if [[ ../Compilers/$compiler.mk -nt $dest_dir/compiler.mk ]]
then
    cp --force ../Compilers/$compiler.mk $dest_dir/compiler.mk
    chmod a-w $dest_dir/compiler.mk
fi

make -C $dest_dir -f $PWD/GNUmakefile -I$PWD $* libf_dir=$PWD
# "-I" option for "nag_rules.mk"
