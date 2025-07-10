
FOLDER=/home/sympoesis/annotationWD/genomes

echo "Busco and Quast"
for FILE in ${FOLDER}/*.fa
do
BASE=$(basename -s .fa $FILE)
cd /home/sympoesis/annotationWD/
echo ">${BASE} Begin"
echo "${BASE} BUSCO AND QUAST"
mkdir ${BASE}_results
cd /home/sympoesis/annotationWD/${BASE}_results/
mkdir bandq_results
cd bandq_results/
source /home/julianna/miniconda3/bin/activate busco
busco -i $FILE -l ascomycota_odb10 -o busco_$BASE -m genome -c 20 -f
conda deactivate

source /home/julianna/miniconda3/bin/activate quast
quast $FILE -o quast_$BASE
conda deactivate
echo "${BASE} BUSCO AND QUAST END"
done
conda deactivate

echo "Funannotate"
for FILE in ${FOLDER}/*.fa
do
BASE=$(basename -s .fa $FILE)
cd /home/sympoesis/annotationWD/${BASE}_results/
echo "${BASE} Funannotate"
funannotate clean -i $FILE -o ${BASE}_clean.fa
funannotate mask -i ${BASE}_clean.fa --cpus 20 -o ${BASE}_masked.fa
funannotate predict -i ${BASE}_masked.fa -o fun_predict_$BASE --species "$BASE" --protein_evidence ${FOLDER}/Lichen.prot.evidence.pep --busco_seed_species aspergillus_fumigatus --cpus 20
echo "${BASE} Funannotate END"
done
conda deactivate

source /home/julianna/miniconda3/bin/activate antismash
echo "antismash"
for FILE in ${FOLDER}/*.fa
do
BASE=$(basename -s .fa $FILE)
echo "${BASE} ANTISMASH"
cd /home/sympoesis/annotationWD/${BASE}_results/fun_predict_${BASE}/predict_results/
antismash --taxon fungi --cb-general --cb-knownclusters --cb-subclusters --asf \
    --pfam2go --smcog-trees --genefinding-tool none --output-dir ../../antismash_out_${BASE} --output-basename ${BASE} ${BASE}.gbk
echo "${BASE} ANTISMASH END"
done
conda deactivate

#Annotate commands contained in seperate script
bash /home/sympoesis/annotationWD/annotate.sh