
screen -S datasets #Para crear
screen -r prokka #Para regresar
screen -S 2146139.prokka -X quit
screen -ls


##############
###DATASETS###
##############
#Datasets (version: 16.31.0)


setwd ~/Documentos/SP_segundo_intento
screen -S datasets

pwd
nano accessions.txt #Creas el documento donde vas a pegar todas las accessiones que bajaste manualmente de PATRIC




#Ahora ten listo el bucle ejectubale de datasets, lo puedes descargar de esta pagina: https://github.com/geraldmoreno993/genomica_evolutiva/blob/main/command_ncbidatasets.sh


sh command_ncbidatasets.sh accessions_722.txt



#Listar screen y volver a entrar
screen -ls
screen -r 2018479.proyecto_spseudintermedius

#generando archiv JSON (metadata)
xargs -a accessions_722.txt -I {} datasets summary genome accession {} --as-json-lines >> metadata.jsonl


#scp ubigem@100.84.255.113:~/Documentos/proyectoSP /media/gerald/DC604FBB604F9B62/staphylococus_pseidintermedius_patric

#Dataformat
#Una vez tengas tu archivo .jsonl puedes ejecutar dataformat
dataformat tsv genome --help

dataformat tsv genome --inputfile

dataformat tsv genome \
    --fields accession,organism-name, geo-loc-name, \
    --inputfile genomes.jsonl > genome_metadata_with_sequences.tsv

#usar toda la metadata
jq -r 'to_entries | map(.key) | @tsv' genomes.jsonl | head -n 1



##############
###OrthoAnI###
##############
#orthoani 0.6.0 (latest)

conda create --name orthoani
conda activate orthoani
pip install orthoani


orthoani --version


cd /home/ubigem/Documentos/SP_segundo_intento

conda activate mlst
mlst *.fasta

mlst 2.23.0
#Se usa una combinacion de MLST. Orthoani y peso de archivo para filtrar y eliminar aquellas secuencias que no son Staphylococcus pseudintermedius
rm GCA_002627685.1.fasta GCA_002627745.1.fasta GCA_002615325.1.fasta GCA_002615305.1.fasta GCA_007998355.1.fasta GCA_002627705.1.fasta GCA_004769005.1.fasta GCA_002627725.1.fasta GCA_004768965.1.fasta GCA_002615345.1.fasta GCA_002615285.1.fasta GCA_004768985.1.fasta GCA_004790815.1.fasta GCA_004790795.1.fasta GCA_004790875.1.fasta GCA_004790775.1.fasta GCA_004790835.1.fasta GCA_900183575.1.fasta


#Eliminar todas las sesiones de screen
screen -ls | grep Detached | awk '{print $1}' | xargs -I{} screen -S {} -X quit
ps aux | grep screen
screen -ls


#Anotación con Prokka (version 1.14.6)
conda env list
conda create --name prokka
conda activate prokka
conda install -c conda-forge -c bioconda -c defaults prokka

mkdir prokka
cd prokka


for a1 in ../*.fasta
do
    prefix=$(basename $a1 .fasta)
    prokka --cpus 12 $a1 -o ${prefix} --prefix ${prefix} --kingdom Bacteria --genus Staphylococcus --species pseudintermedius;
    echo "Procesando archivo: $a1"
done ;






#copiar todos los .gff al directorio panaroo
mkdir panaroo
find . -type f -name "*.gff" -exec cp {} panaroo/ \;
cd panaroo
ls

#Panaroo (version v1.5.1)
cd /home/ubigem/Documentos/proyectoSP

conda create -n panaroo python=3.9
conda activate panaroo

conda install -c conda-forge -c bioconda -c defaults 'panaroo>=1.3'



mkdir -p results_panaroo

panaroo -i *.gff -o results_panaroo --clean-mode strict -t 11 -a core --aligner mafft


ls

cd ~/Documentos/proyectoSP/panaroo




cd ~/Documentos/SP_segundo_intento/prokka/panaroo/results_panaroo
mkdir phylogeny

cp core_gene_alignment_filtered.aln phylogeny

cd phylogeny



#Gubbins Version: 2.4.1
conda  create --name gubbins
conda activate gubbins
conda install bioconda::gubbins
gubbins -v
cd ~/Documentos/proyectoSP/panaroo/results_panaroo
mkdir gubbins
cp core_gene_alignment.aln gubbins/
cd ~/Documentos/proyectoSP/panaroo/results_panaroo/gubbins
ls
conda activate gubbins


screen -S gubbins
screen -r gubbins
screen -S gubbins -X quit #Para eliminar

run_gubbins.py --prefix gubbins --threads 15 core_gene_alignment_filtered.aln



#Reconstructing ancestral sequences with raxmlHPC-PTHREADS-SSE3...
#raxmlHPC-PTHREADS-SSE3 -f A -p 1 -m GTRGAMMA -s /home/ubigem/Documentos/SP_segundo_intento/prokka/panaroo/results_panaroo/phylogeny/core_gene_alignment_filtered.aln.phylip -n core_gene_alignment_filtered.iteration_5.internal -T 15 -t /home/ubigem/Documentos/SP_segundo_intento/prokka/panaroo/results_panaroo/phylogeny/tmp44ez4g7l/core_gene_alignment_filtered.iteration_5.tre.rooted > /dev/null 2>&1
#...done. Run time: 108348.77 s

#Reinserting gaps into the alignment...
#...done. Run time: 108729.22 s

#Running Gubbins to detect recombinations...
#gubbins -r -v core_gene_alignment_filtered.aln.gaps.vcf -a 100 -b 10000 -f /home/ubigem/Documentos/SP_segundo_intento/prokka/panaroo/results_panaroo/phylogeny/tmp44ez4g7l/core_gene_alignment_filtered.aln -t core_gene_alignment_filtered.iteration_5.tre -m 3 core_gene_alignment_filtered.aln.gaps.snp_sites.aln
#...done. Run time: 108941.25 s

#Checking for convergence...
#...done. Run time: 108942.06 s
#Maximum number of iterations (5) reached.

#Exiting the main loop.

#Creating the final output...
#...finished. Total run time: 108945.49 s


#MLST T.Seeman 2.23.0

conda create --name mlst
screen -S mlst #Para crear
screen -r mlst #Para regresar
conda activate mlst

conda install -c conda-forge -c bioconda -c defaults mlst


mlst --list

cd ~/Documentos/proyectoSP
mkdir mlst

cd ~/Documentos/SP_segundo_intento
mkdir mlst


for fasta in *.fasta; do
        mlst --scheme spseudintermedius "$fasta" | cut -f 1,3 >> mlst/ST.tsv
done

mlst --scheme spseudintermedius GCA_036409575.1.fasta

mlst --scheme spseudintermedius GCA_016804965.1.fasta | cut -f 1,3 >> mlst/ST.tsv

cd ~/Documentos/SP_segundo_intento/mlst









cd ~/Documentos/proyectoSP/panaroo/results_panaroo/gubbins
mkdir tree
cp core_gene_alignment.aln.snp_sites.aln tree
cd tree
ls #tiene que aparecer solamente este archivo core_gene_alignment.aln.gaps.snp_sites.aln y arrancamos con el tree



mkdir amrfinderplus
mkdir vfdb
#AMRfinder plus y vfdb

conda create --name abricate

#instalar abricate a desde https://github.com/tseemann/abricate

conda install -c conda-forge -c bioconda -c defaults abricate
abricate --check
abricate --list

#abricate 1.0.1
#vfdb    2597    nucl    2023-Nov-4 
screen -S vfdb #*
screen -S 3212769.vfdb -X quit #Para eliminar

mkdir vfdb

conda activate abricate

##########Ejecutar todo este bloque#########################
for f in *.fasta; do 
        abricate --db vfdb --mincov 80 --minid 80 "$f" > "vfdb/${f/%.fasta}."tab;  
done;


######################## Aqui me quede

screen -r vfdb #Para regresar 
# Combina todos los resultados generados en la carpeta "vfdb" en un archivo único
abricate --summary vfdb/*.tab > vfdb/summary_vfdb.tab



############AMRFinderPlus v4.0.19#################################################


conda create -y -c conda-forge -c bioconda -n amrfinder --strict-channel-priority ncbi-amrfinderplus

screen -S amrfinder #*
screen -r amrfinder #Para regresar
mkdir amrfinderplus

conda activate amrfinder
amrfinder -u


# Iterar sobre los archivos .fasta
for file in *.fasta; do
    # Extraer el nombre base del archivo (sin la ruta y la extensión .fasta)
    base_name=$(basename "$file" .fasta)
    
    # Crear el archivo de salida con el nombre extraído
    output="${base_name}_output.tsv"
    
    # Ejecutar AMRFinderPlus
    amrfinder -n "$file" -o "$output" \
              --organism Staphylococcus_pseudintermedius \
              --plus -i 0.7 -c 0.6 \
              --mutation_all "${base_name}_mutations.tsv"
    
    # Mover los archivos TSV a la carpeta "resultados_tsv"
    mv "${base_name}_output.tsv" "${base_name}_mutations.tsv" amrfinderplus/
done


#Ejecutar esto cuando te encuentres arriba de la carpeta de los tsv
head -n 1 amrfinderplus/*output.tsv | head -n 1 > amrfinderplus/todos_los_genes.tsv
tail -n +2 -q amrfinderplus/*output.tsv >> amrfinderplus/todos_los_genes.tsv

head -n 1 amrfinderplus/*mutations.tsv | head -n 1 > amrfinderplus/todas_las_mutaciones.tsv
tail -n +2 -q amrfinderplus/*mutations.tsv >> amrfinderplus/todas_las_mutaciones.tsv



##########Ejecutar todo este bloque#########################
for f in *.fasta; do 
        abricate --db card --mincov 80 --minid 80 "$f" > "card/${f/%.fasta}."tab;  
done;

# Combina todos los resultados generados en la carpeta "vfdb" en un archivo único
abricate --summary card/*.tab > card/summary_card.tab




# snippy 4.6.0###########################################################

conda create --name snippy

conda activate snippy
conda install -c conda-forge -c bioconda -c defaults snippy

mkdir snippy


#Descargar genoma de referencia S.pseudintermedius ST71 mecA+

datasets download genome accession GCF_000478385.1 --include genome,gbff #al final necesitaba el gbk por eso lo anote con prokka

prokka --cpus 12 GCF_000478385.1_MRSP_E140_v1.0_genomic.fna --prefix GCF_000478385.1 --kingdom Bacteria --genus Staphylococcus --species pseudintermedius

#Ejecutar snippy

screen -S snippy

#snippy-multi input.tab --ref ref.gbk --cpus 16



#snippy --cpus 8 --outdir mysnps_prueba --ref ref.gbk --ctgs mutant.fasta

for f in *.fasta; do 
    snippy --outdir "${f%.fasta}" --ref GCF_000478385.1.gbk --ctgs "$f"
done

#snippy-core CFSAN033867 CFSAN033868 CFSAN033888 CFSAN033889 CFSAN033890 CFSAN033891 CFSAN033892 CFSAN033893 CFSAN033894 CFSAN033895 CFSAN033896 CFSAN033897 CFSAN033899 CFSAN033901 CFSAN033903 CFSAN033904 CFSAN033905 CFSAN033911 CFSAN033912 CFSAN033913 CFSAN033916 CFSAN033917 CFSAN033923 CFSAN033924 CFSAN033925 CFSAN033926 CFSAN033927 CFSAN033928 CFSAN033937 CFSAN033939

mkdir vcf
find . -maxdepth 2 -type f -name "snps.vcf" -exec sh -c 'cp "$1" "vcf/$(basename $(dirname "$1"))_snps.vcf"' _ {} \;
cd vcf 
ls
for f in *snps.vcf; do mv "$f" "${f/_snps.vcf/.vcf}"; done

#Ejecutar el snippy-core\
cd ~/Documentos/SP_segundo_intento/snippy
#El ejecutable runme_2.sh tiene este comando: snippy-core --ref GCF_000478385.1.fasta GCA_027749375.1 y todas las carpetas que se han creado
chmod +x runme_2.sh
sh runme_2.sh


#vcftools para filtrar variantes con frecuencia menor a 0.05

conda create --name vcftools
conda activate vcftools
conda install bioconda::vcftools

vcftools --vcf core.vcf --recode --maf 0.05 --out core_fil 
#Eliminar manualmente l;as 3 primeras lineas del archivo core_fil_recode.vcf
#El archivo modificado se guarda como snp_vcf.txt 


###Isescanv1.7.2.3############################################################################

screen -S iscan
conda create -y --name isescan isescan && conda activate isescan
mkdir isescan
isescan.py --version


# Recorrer todos los archivos .fna en la carpeta de genomas

#crear seqfile
ls *.fasta > seqfile.txt
sed -i 's|^|./|' seqfile.txt



#Isescan para multiples archivos
mkdir -p isescan  

for file in *.fasta; do
    base_name=$(basename "$file" .fasta)
    isescan.py --seqfile "$file" --output "isescan/$base_name" --nthread 6
done



screen -r iscan


#IQ-tree
conda create -n iqtree
conda activate iqtree

conda install bioconda::iqtree

mkdir gtr
cp core_gene_alignment.aln.snp_sites.aln gtr
screen -S tree

#Model test
iqtree -s core_gene_alignment.aln.snp_sites.aln -m TEST -nt AUTO

#Resultado: 88 modelos evaluados donde ...
#Perform fast likelihood tree search using GTR+ASC+G model...
#Estimate model parameters (epsilon = 5.000)
#Perform nearest neighbor interchange...
#Estimate model parameters (epsilon = 1.000)
#1. Initial log-likelihood: -8048143.587
#Optimal log-likelihood: -8048142.653
#Rate parameters:  A-C: 1.07479  A-G: 5.49495  A-T: 1.41623  C-G: 0.65110  C-T: 5.61111  G-T: 1.00000
#Base frequencies:  A: 0.262  C: 0.245  G: 0.237  T: 0.256
#Gamma shape alpha: 0.213
#Parameters optimization took 1 rounds (112.806 sec)
#Time for fast ML tree search: 1590.871 seconds

#NOTE: ModelFinder requires 20932 MB RAM!
#ModelFinder will test up to 88 DNA models (sample size: 390360) ...
# No. Model         -LnL         df  AIC          AICc         BIC
#  1  GTR+F         9853009.087  1387 19708792.175 19708802.073 19723875.556
#  2  GTR+F+ASC     9847677.971  1387 19698129.941 19698139.840 19713213.323
#  3  GTR+F+G4      8183336.758  1388 16369449.517 16369459.430 16384543.774
#  4  GTR+F+ASC+G4  8048143.792  1388 16099063.584 16099073.497 16114157.840
#  7  SYM+G4        8152172.679  1385 16307115.359 16307125.229 16322176.991
#  8  SYM+ASC+G4    8049488.018  1385 16101746.035 16101755.905 16116807.668
# 11  TVM+F+G4      8150762.880  1387 16304299.760 16304309.659 16319383.142
# 12  TVM+F+ASC+G4  8048175.764  1387 16099125.528 16099135.427 16114208.910
# 15  TVMe+G4       8152210.423  1384 16307188.846 16307198.702 16322239.603
# 16  TVMe+ASC+G4   8049557.510  1384 16101883.020 16101892.875 16116933.777
# 19  TIM3+F+G4     8157815.216  1386 16318402.431 16318412.316 16333474.938
# 20  TIM3+F+ASC+G4 8055317.067  1386 16113406.134 16113416.018 16128478.641
# 23  TIM3e+G4      8162326.401  1383 16327418.801 16327428.643 16342458.684
# 24  TIM3e+ASC+G4  8059725.875  1383 16122217.750 16122227.592 16137257.633
# 27  TIM2+F+G4     8155550.520  1386 16313873.040 16313882.925 16328945.547
# 28  TIM2+F+ASC+G4 8052894.977  1386 16108561.954 16108571.838 16123634.461
# 31  TIM2e+G4      8158048.214  1383 16318862.428 16318872.270 16333902.310
# 32  TIM2e+ASC+G4  8055286.016  1383 16113338.032 16113347.874 16128377.915
# 35  TIM+F+G4      8162412.784  1386 16327597.569 16327607.453 16342670.076


#Conclusión

#Best Model:

#GTR+F+ASC+G4 (row 4):
#Log-Likelihood: -8048143.792
#AIC: 16099063.584
#AICc: 16099073.497
#BIC: 16114157.840

#El modelo GTR+F+ASC+G4 es superior porque incorpora una corrección por sesgo de sitios invariantes (ASC), 
#además de la heterogeneidad en las tasas de sustitución (Gamma). Esto le permite ajustarse mejor a los datos 
#en comparación con el modelo GTR+F+G4 o cualquier otro modelo probado.

#Arbol (comenzar a las 6:48 pm 15 dic 2024)
#IQ-TREE multicore version 2.1.4-beta COVID-edition for Linux 64-bit built Jun 24 2021

cd ~/Documentos/proyectoSP/panaroo/results_panaroo/gubbins/tree/gtr
iqtree -s core_gene_alignment.aln.snp_sites.aln -m GTR+F+ASC+G4 -bb 1000 -nt 15

screen -r tree

#la sesión screen se llama 4544.tree con abreviación tree



#SNPsites snp-sites 2.5.1
conda create --name snpsites
conda activate snpsites
conda install bioconda::snp-sites
snp-sites -V
cd ~/Documentos/proyectoSP/panaroo/results_panaroo
mkdir ./../../snpsites




##Nullarbor

conda create --name nullarbor
screen -S nullarbor #Para crear
conda activate nullarbor
conda install -c conda-forge -c bioconda -c defaults nullarbor

nullarbor.pl --check

#KRAKEN
cd ~/Documentos/instaladores
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
tar -C $HOME -zxvf minikraken_20171019_8GB.tgz
#KRAKEN2
screen -S nullarbor2 #Para crear
screen -R nullarbor
#Ctr + D ... screen terminating
conda activate nullarbor

wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz
tar -C $HOME -zxvf minikraken2_v2_8GB_201904.tgz
cd kraken2

mkdir centrifugue
cd /home/ubigem/Documentos/instaladores/centrifugue
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz
mkdir $HOME/centrifuge-db
tar -C $HOME/centrifuge-db -zxvf p_compressed+h+v.tar.gz






