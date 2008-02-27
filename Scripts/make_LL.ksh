# This is a script in KornShell 88 with directives for LoadLeveler.

#@ cpu_limit = 1:00:00

#@ job_name = make_gcm
#@ output = $(job_name).out
#@ error = $(output)

#@ queue

set -x
trap 'exit 1' ERR

cd ~/LMDZE_program/libf
time gmake gcm
