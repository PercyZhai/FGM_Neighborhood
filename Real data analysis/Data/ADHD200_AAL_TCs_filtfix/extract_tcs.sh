#!/bin/bash

templates=( /home/drcc/ADHD200/templates/aal_mask_pad.nii.gz \
            /home/drcc/ADHD200/templates/ez_mask_pad.nii.gz \
            /home/drcc/ADHD200/templates/ho_mask_pad.nii.gz \
            /home/drcc/ADHD200/templates/tt_mask_pad.nii.gz \
            /home/drcc/ADHD200/templates/ADHD200_parcellate_200.nii.gz \
            /home/drcc/ADHD200/templates/ADHD200_parcellate_400.nii.gz )
ids=( aal ez ho tt cc200 cc400 )

for restfile in s*rest*.nii.gz
do
    for i in $( seq 0 $((${#templates[@]} - 1 )))
    do
        echo "$restfile ${templates[$i]} --> ${ids[$i]}";
        3dROIstats -mask ${templates[$i]} ${restfile} > ${restfile%%.nii.gz}_${ids[$i]}_TCs.1D
    done

    if [[ ${restfile} == sn* ]]
    then
        echo "fsl multiregression $restfile"
        fsl_glm -i ${restfile} -d ${restfile%%.nii.gz}_TCs.1D -o fc2_${restfile}
    else
        echo "skipping fsl multiregression $restfile"
    fi
done 
