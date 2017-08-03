#!/usr/bin/perl -w

sub usage
{
        print "\nEste script aplica la función sg_filtering_hits_by_mismatches_Q.pl sobre las muestras presentes en el argumento <ruta>.\n";
	print "Los valores por defecto que se pasan a la función sg_filtering_hits_by_mismatches son 5 y 10 para F3, y 3 y 10 para F5 (además del archivo csfasta.ma).\n";
	print "Si se requieren otros valores, debe modificarse el presente código.\n";
	print "Se realiza una copia de los ficheros csfasta.ma de partida, poniéndoles la terminación _original al final.\n";
	print "Si en uno de los directorios se detecta que ya existe un fichero con la terminación .csfasta.ma_original, no se sobreescribe, si no que se emplea como fichero de partida sobre el que se aplica el filtro, guardándolo en el correspondiente .csfasta.ma \n";
        print "\nCOMO SE USA: sg_filtrar_por_mismatch_Qvalues_muestra.pl_ <ruta>\n";
        print "ejemplo: sg_filtrar_por_mismatch_Qvalues_muestra.pl /path_a_directorio_con_muestras \n\n";
        print "INPUT: El path que contiene las diferentes muestras (por ejemplo, /data/results/Solid0065/BF11_Cardio-cancer/cancer). A ser posible, no incluir el carácter / del final del path, pues puede causar problemas.";
	print "Es necesario que dentro del directorio se tengan las siguientes estructuras de carpetas: \n"; 
	print "- path/muestra/secondary/bioscope_XXXXXX/output/F3/s_mapping/\n";
	print "- path/muestra/secondary/bioscope_XXXXXX/output/F5/s_mapping/\n";
	print "Se buscan los ficheros .csfasta.ma en el interior de dichas carpetas, y se les aplica el filtro. Si se detecta que ya se había ejecutado el filtro, no se sobreescribe el fichero .csfasta.ma_original \n";
        print "OUTPUT: Ficheros filtrados sustituyendo al existente .csfasta.ma y ficheros originales con extensión .csfasta.ma_original. Además se imprimenlas incidencias encontradas, como por ejemplo cuando ya existe el fichero .csfasta.ma_original, el nombre de los ficheros .csfasta.ma filtrados. \n";
        exit(1);
}

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
        usage();
}

# Check arguments and extensions of the files
if (substr($ARGV[0],length($ARGV[0]),1) eq '/')
{
	$samples_path=scalar($ARGV[0]);
}
else
{
	$samples_path=scalar(($ARGV[0]) . "/");
}

# main body of the programm

# First of all we list the samples in the path directory
opendir(DIR, "$samples_path");

@samples= readdir(DIR);
closedir(DIR);

# Now each folder in the path directory is analyzed
foreach $sample (@samples) {
	#print (-d $sample) . "\n";

	if ((-d ($samples_path . $sample)) && $sample ne '.' && $sample ne '..'){
		#print "$samples_path" . "$sample" . "\n";	
		
		# If the proper file structure is found, then the filtering is carried out on the .csfasta files
		#Needed file structure: $samples_path/secondary/bioscope_110216/output/F3/s_mapping
		#			$samples_path/secondary/bioscope_110216/output/F5/s_mapping
		
		$full_name_sec="$samples_path" . "$sample" . "/secondary/";

		opendir(DIR,"$full_name_sec");
		
		@bioscopes= readdir(DIR);
		closedir(DIR);
		
		# Now the bioscope folders are analyzed
		
		@strands= ("F3", "F5");

		foreach $bioscope (@bioscopes) {
			#print "$samples_path" . "$sample" . "/secondary/" . "$bioscope" . "\n";

			

			if ((-d "$full_name_sec" . "$bioscope") && $bioscope ne '.' && $bioscope ne '..'){

				foreach $strand (@strands) {

					if ($strand eq "F3")
					{
						$MismThreshold = 5;
						$QThreshold = 10;
					}
					else
					{
						$MismThreshold = 2;
                                                $QThreshold = 10;
					}

					$full_name= "$full_name_sec" . "$bioscope" . "/output/" . "$strand" . "/s_mapping/";
					opendir(DIR,"$full_name");

					@map_files= readdir(DIR);
					closedir(DIR);

					# Finally the csfasta files will be filtered
					foreach $map_file (@map_files){
						if (substr($map_file,length($map_file)-10,10) eq "csfasta.ma")
						{
							$full_name_file= "$full_name" . "$map_file";
							$full_name_renamed_file= "$full_name_file" . "_original";
							
							# It will be checked if the file was already filtered:
							if (-e $full_name_renamed_file)
							{
								print "The file $full_name_file was already filtered! \n";
							}
							else
							{
							
								# The file to process is copied with "_original" appended at its end
								@args=("mv", "$full_name_file", "$full_name_renamed_file");
								system (@args)==0 or die "Problemas al mover $full_name_file \n";
 							}
								# The task to filter the file is released
								print "Filtering file $full_name_file with a Mismatch Threshold of $MismThreshold and a Quality Threshold of $QThreshold  \n";
								system ("sg_filtering_hits_by_mismatches_Q.pl $full_name_renamed_file $MismThreshold $QThreshold > $full_name_file");
								
								#system (nohup sg_filtering_hits_by_mismatches_Q.pl solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F3_07S1148.csfasta.ma_original 5 10 > solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F3_07S1148.csfasta.ma &);
								#print "$samples_path" . "$sample" . "/secondary/" . "$bioscope" . "/output/F3/s_mapping/" . "$map_file" . "\n";
							#}
						}
					}

				}

				#print $bioscope . "\n";
				
				#nohup sg_filtering_hits_by_mismatches_Q.pl solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F3_07S1148.csfasta.ma_original 5 10 > solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F3_07S1148.csfasta.ma &

			}
		}
	}
}

