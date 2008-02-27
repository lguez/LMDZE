# This is a script in Bash with directives for LoadLeveler.

# This script submits an execution to LoadLeveler on Zahir

##@ cpu_limit = 10:00
# (per process)
##@ data_limit = 1GB # (per process)

#@ job_name = gcm_79
##@ output = $(job_name).out
##@ error = $(output)

#@ queue

set -vx
trap 'exit 1' ERR
echo $SHELL
echo $PATH
echo $-
echo $PS1

cd $workdir/LMDZE
./gcm.sh <<EOF
igcm28
34
n

EOF
