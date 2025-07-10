FOLDER=/home/sympoesis/annotationWD/genomes

source /home/julianna/miniconda3/bin/activate funa2
echo "final annotation with funannotate"
for FILE in ${FOLDER}/*.fa
do
BASE=$(basename -s .fa $FILE)
echo "${BASE} FINAL ANNOTATION"
cd /home/sympoesis/annotationWD/${BASE}_results/
#eggnog should run by default...
funannotate annotate -i fun_predict_${BASE}/ --antismash antismash_out_${BASE}/${BASE}.gbk --cpus 20 # --sbt
echo "${BASE} FINAL ANNOTATION END"
done
echo "all completed steps"
conda deactivate