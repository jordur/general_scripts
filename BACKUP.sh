#! /bin/bash
## Copia de seguretat a GitHub
read -p "Vols fer backup a GitHub (y/n)?" dev
	case $dev in
		[Nn]* ) continue;;
		[Yy]* ) cd /media/bec2-jcalvete/Disc_Dur/Feina_Jordi/scripts
			read -p  "Escriu comentari per commit: " com
			git pull origin master
			git add --ignore-removal .
    			git commit -m "$com"
    			git push origin master; break;;
		* ) "Digues si (y) o no (n)"
		
	esac

## Copia de seguretat a /Elements
cd
if [ -f /media/bec2-jcalvete/Elements/autorun.inf ];then
	read -p "Vols fer BACKUP a /media/Elements (y/n)?" dev
		case $dev in
			[Nn]* ) exit;;
			[Yy]* ) cd /media/bec2-jcalvete/Disc_Dur/Feina_Jordi/glandulas_mejicanas; rsync -avz --exclude '*.bam' --exclude '*.sam' --exclude '*.fastq' --exclude '*.fasta' * /media/bec2-jcalvete/Elements; break;;
			* ) "Digues si (y) o no (n)"
		esac
else
	echo "WARNING : No existeix /media/Elements"
fi
exit
