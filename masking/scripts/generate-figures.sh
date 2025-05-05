#!/bin/bash

PLOTSAMPLES=50
MCSAMPLES=100
MCELEMENTS=20
VIEWANGLES="0 10 20 30 40 50 60 70 80"

declare -A BARBULE_CONFIGURATIONS
BARBULE_CONFIGURATIONS=( ["opaque"]="--excentricity=1 --separation=0"
                      ["normal"]="--excentricity=1 --separation=1"
                      ["transmittance"]="--excentricity=1 --separation=10" 
                      ["excentric"]="--excentricity=2 --separation=1")

for config in "${!BARBULE_CONFIGURATIONS[@]}"
do
    echo Barbule - ${config}
    echo Analytical
    time ../bin/plot-barbules --output=../bin/barbule-${config}-plot.svg ${BARBULE_CONFIGURATIONS[$config]} --plot-samples=${PLOTSAMPLES}
    echo
    echo Monte-Carlo ${MCSAMPLES}spp
    time ../bin/plot-barbules --output=../bin/barbule-${config}-montecarlo-plot.svg ${BARBULE_CONFIGURATIONS[$config]} --plot-samples=${PLOTSAMPLES} --montecarlo-barbules=${MCELEMENTS} --montecarlo-samples=${MCSAMPLES}
    for angle in ${VIEWANGLES}
    do
        ../bin/barbule-masking --output=../bin/barbule-${config}-view-${angle}.svg --angle=${angle}.01 ${BARBULE_CONFIGURATIONS[$config]}
    done
    echo
    echo
done 

declare -A BARB_CONFIGURATIONS
BARB_CONFIGURATIONS=( ["planar"]="--excentricity=1 --barbule-length=2 --barbule-inclination=0"
                      ["planartransparent"]="--excentricity=1 --barbule-length=2 --barbule-inclination=0 --left-transparency=0.5 --right-transparency=0.5"
                      ["planarhalftransparent"]="--excentricity=1 --barbule-length=2 --barbule-inclination=0 --left-transparency=0.5 --right-transparency=0"
                      ["upwards"]="--excentricity=1 --barbule-length=1 --barbule-inclination=45"
                      ["upwardstransparent"]="--excentricity=1 --barbule-length=1 --barbule-inclination=45 --left-transparency=0.5 --right-transparency=0.5"
                      ["upwardshalftransparent"]="--excentricity=1 --barbule-length=1 --barbule-inclination=45 --left-transparency=0.5 --right-transparency=0" 
                      ["downwards"]="--excentricity=1 --barbule-length=1 --barbule-inclination=-45"
                      ["excentric"]="--excentricity=3 --barbule-length=1 --barbule-inclination=45"
                      ["unexcentric"]="--excentricity=0.3 --barbule-length=1 --barbule-inclination=45" 
                      ["unexcentrictransparent"]="--excentricity=0.3 --barbule-length=1 --barbule-inclination=45 --left-transparency=0.5 --right-transparency=0.5" 
                      ["unexcentrichalftransparent"]="--excentricity=0.3 --barbule-length=1 --barbule-inclination=45 --left-transparency=0.5 --right-transparency=0"
                      ["longbarbules"]="--excentricity=1.5 --barbule-length=4 --barbule-inclination=45 --left-transparency=0.5 --right-transparency=0.5"  )
for config in "${!BARB_CONFIGURATIONS[@]}"
do
    echo Barb - ${config}
    echo Analytical
    time ../bin/plot-barbs --output=../bin/barb-${config}-plot.svg ${BARB_CONFIGURATIONS[$config]}  --plot-samples=${PLOTSAMPLES}
    echo
    echo Monte-Carlo ${MCSAMPLES}spp
    time ../bin/plot-barbs --output=../bin/barb-${config}-montecarlo-plot.svg ${BARB_CONFIGURATIONS[$config]} --plot-samples=${PLOTSAMPLES} --montecarlo-barbs=${MCELEMENTS} --montecarlo-samples=${MCSAMPLES}
    for angle in ${VIEWANGLES}
    do
        ../bin/barb-masking --output=../bin/barb-${config}-view-${angle}.svg --angle=${angle}.01 ${BARB_CONFIGURATIONS[$config]}
    done
    echo
    echo
done 

declare -A FEATHER_CONFIGURATIONS
FEATHER_CONFIGURATIONS=( ["standard"]="--barb-excentricity=1 --barbule-excentricity=1 --barbule-length=1 --barbule-rotation=45 --barbule-inclination=0 --barbule-separation=1"
                      ["strongbarb"]="--barb-excentricity=3 --barbule-excentricity=1 --barbule-length=1 --barbule-rotation=45 --barbule-inclination=0 --barbule-separation=1"
                      ["strongbarbule"]="--barb-excentricity=1 --barbule-excentricity=1 --barbule-length=3 --barbule-rotation=45 --barbule-inclination=45 --barbule-separation=1"
)
for config in "${!FEATHER_CONFIGURATIONS[@]}"
do
    echo FEATHER - ${config}
    echo Analytical
    time ../bin/plot-feathers --output=../bin/feather-${config}-plot.svg ${FEATHER_CONFIGURATIONS[$config]}  --plot-samples=${PLOTSAMPLES}
    echo
    echo Monte-Carlo ${MCSAMPLES}spp
    time ../bin/plot-feathers --output=../bin/feather-${config}-montecarlo-plot.svg ${FEATHER_CONFIGURATIONS[$config]}  --plot-samples=${PLOTSAMPLES} --montecarlo-samples=${MCSAMPLES} 
    echo
    echo
done 



