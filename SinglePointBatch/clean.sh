#!bin/bash

# Command line input
usage="Remove unnecessary Columbus7 files and directories for collecting in current directory

$(basename "$0") [-h]

optional arguments:
  -h             show this help message and exit"

while getopts ':h:' option; do
  case "$option" in
    h) echo "$usage"
       exit;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit;;
  esac
done

# Do the job
for entry in * ; do
    if [ -d $entry ]; then # This is a directory
        cd $entry
        if [ -f 'runc.log' ]; then # This directory has a started job
            success=`cat runc.log |grep timings`
            if [ -n "$success" ]; then # This job has succeeded
                rm *
                rm -r COSMO MOCOEFS MOLDEN RESTART
            fi
        fi
        cd ..
    fi
done