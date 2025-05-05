#!/bin/bash

PLOTSAMPLES=50
MCSAMPLES=100
MCELEMENTS=20

VIEWANGLES="25 70"
BARBULE_CONFIG="--excentricity=2 --separation=1"

echo Analytical
time ../bin/plot-barbules --output=../bin/barbule-${config}-plot.svg ${BARBULE_CONFIG} --plot-samples=${PLOTSAMPLES}
echo
echo Monte-Carlo ${MCSAMPLES}spp
time ../bin/plot-barbules --output=../bin/barbule-${config}-montecarlo-plot.svg ${BARBULE_CONFIG} --plot-samples=${PLOTSAMPLES} --montecarlo-barbules=${MCELEMENTS} --montecarlo-samples=${MCSAMPLES}

for angle in ${VIEWANGLES}
do
    ../bin/barbule-masking --output=../bin/barbule-masking-${angle}.svg --angle=${angle}.01 ${BARBULE_CONFIG}
done

BARB_CONFIG="--excentricity=1.5 --barbule-length=2.75 --barbule-inclination=45 --left-transparency=0.5 --right-transparency=0.5"

echo Analytical
time ../bin/plot-barbs --output=../bin/barb-${config}-plot.svg ${BARB_CONFIG}  --plot-samples=${PLOTSAMPLES}
echo
echo Monte-Carlo ${MCSAMPLES}spp
time ../bin/plot-barbs --output=../bin/barb-${config}-montecarlo-plot.svg ${BARB_CONFIG} --plot-samples=${PLOTSAMPLES} --montecarlo-barbs=${MCELEMENTS} --montecarlo-samples=${MCSAMPLES}

for angle in ${VIEWANGLES}
do
    ../bin/barb-masking --output=../bin/barb-masking-${angle}.svg --angle=${angle}.01 ${BARB_CONFIG}
done


declare -A FEATHER_CONFIGURATIONS
FEATHER_CONFIGURATIONS=(["strongbarb"]="--barb-excentricity=3 --barbule-excentricity=1 --barbule-length=1 --barbule-rotation=45 --barbule-inclination=0 --barbule-separation=1"
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



