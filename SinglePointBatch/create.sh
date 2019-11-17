#!bin/bash

# Command REPLY input
usage="Create Columbus7 single point job directories naming from 1 to the number of geometries in current directory

$(basename "$0") [-h]
                 InputPath NAtoms GeomPath

positional arguments:
  InputPath      location of Columbus input template
  NAtoms         number of atoms in the molecule
  GeomPath       location of the appended geometry file

optional arguments:
  -h             show this help message and exit"

while getopts ':hs:' option; do
  case "$option" in
    h) echo "$usage"
       exit;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit;;
  esac
done

# Do the jub
i=1
count=1
while read; do
    if [ $count == 1 ]; then # Start a new geometry
        if [ ! -d $i ]; then
            mkdir $i
        fi
        cd $i
        cp $1/* .
        echo "$REPLY" > geom
        count=$[count+1]
    elif [ $count == $2 ]; then # End of current geometry
        echo "$REPLY" >> geom
        cd ..
        i=$[i+1]
        count=1
    else # Keep working on current geometry
        echo "$REPLY" >> geom
        count=$[count+1]
    fi
done < $3
