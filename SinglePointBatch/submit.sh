#!bin/bash

# Command line input
usage="Submit new Columbus7 job or resubmit failed job in current directory

$(basename "$0") [-h] [-s n]
                 JobScriptPath

positional arguments:
  JobScriptPath  location of the job script for your queuing system

optional arguments:
  -h             show this help message and exit
  -s             queuing system (default = slurm)"

queue='slurm' # We will read the queuing system log to determine whether a job is failed
while getopts ':hs:' option; do
  case "$option" in
    h) echo "$usage"
       exit;;
    s) queue=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit;;
  esac
done
shift $((OPTIND - 1))

JobScript=$(basename $1)

# Do the job
for entry in * ; do
    if [ -d $entry ]; then # This is a directory
        cd $entry
        if [ -f 'geom' ]; then # This is a job directory
            if [ ! -f 'runc.log' ]; then # This is a new job
                cp $1 .
                sbatch $JobScript
            else # This job has started
                success=`cat runc.log |grep timings`
                if [ ! -n "$success" ]; then # This job has not yet succeeded
                    error=`cat runc.log |grep 'Error occured!'`
                    timeout=`cat ${queue}* |grep TIMEOUT`
                    if [ -n "${error}${timeout}" ]; then # This job has failed
                        echo 'Failed job '$entry
                        if [ -d WORK ]; then # Remove Columbus7 scratch
                            rm -r WORK
                        fi
                        rm ${queue}* # Remove queuing system log
                        sbatch $JobScript
                    fi
                fi
            fi
        fi
        cd ..
    fi
done