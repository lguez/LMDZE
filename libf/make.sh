# This is a script in Bash.

# This script builds LMDZE. It is useful only to build into a
# directory other than "libf". If you want to build into "libf", you
# can just invoke "make" directly.

dest_dir=/home/guez/In_transit/LMDZE_work/Compil_prod
##dest_dir=$workdir/LMDZE/Compil_prod

set -xe
gmake -C $dest_dir -f $PWD/GNUmakefile -I$PWD $* libf_dir=$PWD
