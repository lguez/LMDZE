# This is a script in Bash.

# This script runs "gcm".

trap 'echo -e \\a; exit 1' ERR
set -x

# This is only useful if the program was compiled with G95:
G95_FPU_UNDERFLOW=${G95_FPU_UNDERFLOW:+No}
G95_FPU_INVALID=${G95_FPU_INVALID:+No}

date
time gcm <gcm_nml.txt >gcm_out.txt 2>gcm_err.txt
echo -e '\a' # beep
