# This is a script in Bash.

# This script allows to build LMDZE into a directory other than
# "libf". Optionally, a compiler may be chosen on the command line.

# -- If a compiler is chosen on the command line then the script links
# a compiler file to the directory "libf" and builds LMDZE into a
# directory associated to the chosen compiler.

# -- If the compiler is not chosen on the command line and a compiler
# file already exists in "libf" then the script builds LMDZE with that
# compiler into the directory associated to it.

# -- If the compiler is not chosen on the command line and there is no
# compiler file in "libf" then the script builds LMDZE with the
# default compiler macros of make into a default directory. The
# default directory is defined in the script.

# If you just want to build into "libf" with the current compiler
# file, you can invoke "make" directly.

# Usage:

# make.sh [-c <<<compiler>>>] [options and arguments for make]

# The file "<<<compiler>>>.mk" should exist in "../Compilers" and the
# directory "<<<compiler>>>" shoud exist in "$dest_dir".

# Note that "-c" was chosen because it is not an option of "make".

getopts :c: name
if ((($? == 0)) && [[ $name = c ]])
    then
    echo "Linking \"$OPTARG.mk\"..."
    ln -sf ../Compilers/$OPTARG.mk compiler.mk
    target=$OPTARG
    shift $((OPTIND - 1))
else
    # Find the current compiler, if any:
    if [[ -L compiler.mk ]]
	then
	target=`basename $(readlink compiler.mk) .mk`
    fi
fi

##dest_dir=$workdir/In_transit/LMDZE/Compil_prod
dest_dir=/usr/local/guez/LMDZ/LMDZE_work/Compil_prod${target:+_$target}
# (Do not just use the name of the compiler as a directory name, it
# confuses some compilers.)

set -xe
gmake -C $dest_dir -f $PWD/GNUmakefile -I$PWD $* libf_dir=$PWD
