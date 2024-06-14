#!/bin/bash

list=${1}
for f in $(cat ${list}) 
  do
  sample_id=$( echo ${f} | sed -e 's/^.*A/A/g' |sed -e 's/.tar$//g' |uniq)	
  id=$( echo ${f}|sed -e 's/.tar$//g' |uniq)
	
  if [[ $(echo ${f} | tr "_" "\t" | cut -f 2 |sed -e 's/.tar$//g') != "${sample_id}" ]]; then
      batch=$(ls ${f} | tr "_" "\t" | cut -f 2)
      sample=$(echo ${sample_id}"_"${batch})
      echo ${batch}
      echo "sample =" ${sample}
    else 
		  echo "sample =" ${sample}
  fi
	
	sshpass -p $(cat ${dir}) \
	rsync -avP ./${f} -e \
	ssh yh@172.:/mnt/nas2/GBM/ONT/${sample}/
	wait
	echo $sample_id "upload done" >> ${id}.txt
	sshpass -p $(cat ${dir}) \
	rsync -avP ./${id}.txt -e \
	ssh yh@172:/mnt/nas2/GBM/ONT/${sample}/
done &
