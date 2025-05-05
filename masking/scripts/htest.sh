#!/bin/bash

VIEWANGLES="0 10 20 30 40 50 60 70 80"
declare -A CONFIGURATIONS
CONFIGURATIONS=( ["normal"]="--excentricity=1"
                 ["excentric"]="--excentricity=3"
                 ["unexcentric"]="--excentricity=0.3" )
for config in "${!CONFIGURATIONS[@]}"
do
    for angle in ${VIEWANGLES}
    do
		../bin/htest --output=../bin/htest-${config}-${angle}.svg ${CONFIGURATIONS[$config]} --angle=${angle}
		inkscape --export-type="png" "../bin/htest-${config}-${angle}.svg"
    done
done 




