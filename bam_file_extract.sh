#!/bin/bash

bam_file=${1}

samtools view -@ 12 -F 16 $bam_file |awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$10,$12}' > test_for.bam
samtools view -@ 12 -f 16 $bam_file |awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$10,$12}' > test_rev.bam
echo `date`
for bam in test_for.bam test_rev.bam
        do
            echo create pipe
            [ ! -p tmp ] && mkfifo tmp
            exec 9<>tmp
        
        # fill in the pipe
        for ((i=0;i<12;i++))
            do
                echo >&9
                done

        split -n l/100 ${bam} /tmp/ont_${bam%.bam}_split$$

    for file in /tmp/ont_${bam%.bam}_split$$*
        do
            #echo $file
            {
                read -u 9
                {
                    awk -v name=${file} \
                        '{for (i = 1; i <= NF; i++) { \
                            if (i == 2 || i == 6 || i == NF ) { \
                                if (i == 2) { \
                                    if ( name ~ "rev" ) {strand = "-"} \
                                    else {strand = "+"}; \
                                    {printf strand"\t"}; \
                                    strand = 0 } \
                                else if (i == 6) { \
                                    igar=$i; \
                                    gsub(/[0-9]*H|[0-9]*S|[0-9]*I/,"",cigar); \
                                    gsub("[A-Z]","\n",cigar); \
                                    len = split(cigar, arr, "\n"); \
                                    for (j = 1; j <= len; j++) {sum += arr[j]}; \
                                        mis_cigar=split($i,arr,"D|I"); \
                                        if (mis_cigar == "1" ) {mis_sum = 0} \
                                        else {for (n = 1; n <= mis_cigar ; n++) { \
                                            gsub(/([0-9]*[A-Z])/,"", arr[n]); \
                                            mis_sum += arr[n] }}; \
                                            {printf sum"\t"mis_sum"\t"$i"\t"}; \
                                            sum = 0 ; \
                                            mis_sum = 0 ;\
                                            } \
                                else { \
                                    mis=$i; \
                                    if (mis ~ "NM") {gsub("NM:i:","",mis)} \
                                    else {mis=0} \
                                    {printf mis"\n"} \
                                    } \
                                } \
                             else {printf $i"\t"} \
                                } \
                                }' \
                        $file \
                    | awk 'OFS="\t" { print $1,$3,($4-1),($4+$6-1),$2,($10-$7),$5,$8,$9}' \
                        > ${file}.list


                echo >&9
                        } &
                }
        done
        wait
    done
    wait
echo  close pipe
exec 9>&-
rm tmp
echo `date`
cat /tmp/ont_*_split$$*.list > temp.list
rm /tmp/ont_*_split*
