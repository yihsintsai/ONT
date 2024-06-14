!/bin/bash

bam_file=${1}
out_dir=${2}
[[ ! -f $bam_file ]] && { echo "file not exists"; exit 1; }


samtools view -@ 12 -F 16 $bam_file |awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$10,$12}' > ${bam_file%.bam}_for.bam
samtools view -@ 12 -f 16 $bam_file |awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$10,$12}' > ${bam_file%.bam}_rev.bam
echo `date` "transfer format start"
for bam in  ${bam_file%.bam}_for.bam ${bam_file%.bam}_rev.bam
        do
            echo $bam
            echo create pipe
            [ ! -p /mnt/nas2/yh/tmp ] && mkdir -p /mnt/nas2/yh/tmp
            exec 9<>tmp
        
        # fill in the pipe
        for ((i=0;i<12;i++))
            do
                echo >&9
                done

        split -n l/100 ${bam} /mnt/nas2/yh/tmp/ont_${bam%.bam}_split$$

    for file in /mnt/nas2/yh/tmp/ont_${bam%.bam}_split$$*
        do
            #echo $file
            {
                read -u 9
                {
                    awk -v name=${file} \
                        '{for ( i = 1 ; i <= NF ; i++) { \
                            if ( i == 2 || i == 6 || i == NF ) { \
                                if (i == 2) { \
                                    if ( name ~ "rev" ) {strand = "-"} \
                                    else {strand = "+"}; \
                                    {printf strand"\t"}; \
                                    strand = 0 
                                    } \
                            else if (i == 6) { \
                                cigar=$i; \
                                leng=$i; \
                                gsub(/[0-9]*H|[0-9]*S|[0-9]*I/,"",cigar); \
                                gsub(/[0-9]*H|[0-9]*S|[0-9]*D/,"",leng); \
                                gsub("[A-Z]","\n",leng); \
                                gsub("[A-Z]","\n",cigar); \
                                len = split(cigar, arr, "\n"); \
                                leng_len = split(leng, arr_1, "\n"); \
                                for (j = 1 ; j <= len; j++) {sum += arr[j]}; \
                                for (g = 1 ; g <= leng_len ; g++) {leng_sum += arr_1[g]}; \
                                    mis_cigar=split($i,arr,"D|I"); \
                                    if (mis_cigar == "1" ) {mis_sum = 0} \
                                    else {for (n = 1; n <= mis_cigar ; n++) { \
                                        gsub(/([0-9]*[A-Z])/,"", arr[n]); \
                                        mis_sum += arr[n] }}; \
                                    {printf sum"\t"leng_sum"\t"$i"\t"}; \
                                    sum = 0 ; \
                                    leng_sum =0 ; \
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
                    | awk -v bpsite=${bam%.bam}  'OFS="\t" { print $1,$3,($4-1),($4+$6-1),$2,$7,$5,$8}' \
                        > ${file}.list


                echo >&9
                        } &
                }
        done
        wait
    done
    wait
echo  `date` "transfer format end"
exec 9>&-
rm tmp
echo `date`
cat /mnt/nas2/yh/tmp/ont_${bam_file%.bam}*_split$$*.list > ${out_dir}/${bam_file%.bam}.bed
rm /mnt/nas2/yh/tmp/ont_${bam_file%.bam}*_split$$*.list
rm /mnt/nas2/yh/tmp/ont_${bam_file%.bam}*_split$$*
rm ${bam_file%.bam}_for.bam
rm ${bam_file%.bam}_rev.bam
