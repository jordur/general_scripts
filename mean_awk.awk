### Script per calcular la mitjana de la columna que vulguis: output= mean, min i max ###############
awk '{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} END {print total/count, min, max}'
