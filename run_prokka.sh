rm -r annotations
p='/x/assemblies_downstream'
for file in $p/*.fasta
do
  output=${file/$p/''}
  output=${output/'.fasta'/''}
  sbatch prokka.sh $file $output
done
