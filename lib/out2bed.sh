# bash out2bed.sh xxx.repeatmasker.out > xxx.repeatmasker.bed
cat $1 | awk 'NR>=4 {print $5"\t"$6-1"\t"$7"\t"$11"\t"$10}'
