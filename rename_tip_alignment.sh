alignment="/x/x.fasta"
renamed="${alignment}_renametip.fasta'"

sed 's/#/_/g' $alignment > $renamed
