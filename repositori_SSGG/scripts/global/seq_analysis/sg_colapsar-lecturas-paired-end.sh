#!/bin/bash

EXPECTED_ARGS=3 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments



if [ ${#} -ne $EXPECTED_ARGS ]; then
        echo 'Nombre: sg_colapsar-lecturas-paired-end.sh'
        echo 'Descripcion: Colapsa las lecturas redundantes a partir de un fichero con paired-ends'
        echo 'Parametros:'
        echo '       $1 : Grupo de lecturas 1'
        echo '       $2 : Grupo de lecturas 2'
	echo '       $3 : Directorio de salida'
        echo 'Output: ficheros fasta, uno para cada grupo de lecturas, que contienen el conteo de las lecturas en el nombre de secuencia (numero generado por el script) y una sola copia de cada pareja de secuencias redundantes'
        exit $E_BADARGS
fi

echo $(date) "sg_colapsar-lecturas-paired-end.sh se esta ejecutando"
mkdir collapse_tmp
cd collapse_tmp
paste $1 $2 > tmp1
echo $(date) "paso 1/3"
sg_collapse-reads.pl tmp1 > tmp2 
echo $(date) "paso 2/3"
sort -T . tmp2 > tmp3
valor=$(wc -l tmp3)
echo $(date) "paso 3/3"
sg_comparar-reads.pl tmp3 $valor
mv p1-read1.fastq $3
mv p1-read2.fastq $3
echo $(date) "borrando temporal"
cd ..
rm -rf collapse_tmp
echo $(date) "sg_colapsar-lecturas-paired-end.sh ha terminado"
