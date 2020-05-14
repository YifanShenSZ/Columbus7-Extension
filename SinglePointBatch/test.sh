#!bin/bash

# Command line input
usage="Test whether Columbus7 job has not started or failed in current directory

$(basename "$0") [-h] [-q n]

optional arguments:
  -h             show this help message and exit
  -q             queuing system (default = slurm)"

queue='slurm' # We will read the queuing system log to determine whether a job is failed
while getopts ':hq:' option; do
  case "$option" in
    h) echo "$usage"
       exit;;
    q) queue=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit;;
  esac
done
shift $((OPTIND - 1))

# Do the job
for entry in * ; do
    if [ -d $entry ]; then # This is a directory
        cd $entry
        if [ -f 'geom' ] && [ -f 'control.run' ]; then # This is a job directory
            if [ ! -f 'runc.log' ]; then # This is a new job
                echo $entry' has not started'
            else # This job has started
                success=`cat runc.log |grep timings`
                if [ ! -n "$success" ]; then # This job has not yet succeeded
                    echo $entry' has not finished'
                    error=`cat runc.error |wc -l`
                    timeout=`cat ${queue}* |grep TIMEOUT`
                    if [ $error -gt 0 ] || [ -n "${timeout}" ]; then # This job has failed
                        echo $entry' has failed'
                    fi
                fi
            fi
        fi
        cd ..
    fi
done
