#!/usr/bin/env python
import sys
import time

def obtener_chr(fichero_bed):

	"""Obtiene el numero del chromosoma que estamos procesando
	Esto sirve posteriormente para escribir el nombre del fichero resultante de output
	"""

	fh = open(fichero_bed, 'r')
	linea = (list(fh)[0]).strip().split()
        num_chr = linea[0]
	fh.close()
	return num_chr

def obtener_primera_linea(fichero_bed):

	"""Obtiene los valores del primer intervalo, obtenido del archivo de origen .*bed
	Esto sirve posteriormente para tratar de forma especial a la primera linea
	"""

	fh = open(fichero_bed, 'r')
	linea = (list(fh)[0]).strip().split()
	ant = int(linea[1])
	ant_final = int(linea[2])
	fh.close()
	return ant, ant_final

def obtener_ultimo_pileup(fichero_pileup):

	"""Obtiene el ultimo valor del archivo *.pileup
	Esto sirve posteriormente para detectar si hemos llegado a la ultima linea
	del fichero pileup y poder tratarla de forma especial
	"""

	fh = open(fichero_pileup, 'r')
	linea = (list(fh)[-1]).strip().split()
	ult_pileup = int(linea[1])
	fh.close()
	return ult_pileup

def obtener_ultimo_bed(fichero_bed):

	"""Obtiene el ultimo valor de la ultima linea del ficero *.bed
	Esto sirve posteriormente para tratar de forma especial la ultima linea
	"""

        fh = open(fichero_bed, 'r')
        linea = (list(fh)[-1]).strip().split()	
        ult_intervalo_inicio = int(linea[1])
	ult_intervalo_final = int(linea[2])
	fh.close()
        return ult_intervalo_inicio, ult_intervalo_final 

def obtener_num_lineas_totales_bed (fichero_bed):

	fichero_bed_lineas = open(fichero_bed)
	total_lineas_bed = len(fichero_bed_lineas.readlines())
	return total_lineas_bed - 1

def obtener_array_intervalos(fichero_bed):

	fichero_bed_lineas = open(fichero_bed)
	total_lineas_bed = fichero_bed_lineas.readlines()
	array_lineas_intervalos = []
	for linea in total_lineas_bed:
        	linea_corregida = linea.strip().split()
	        array_lineas_intervalos.append(linea_corregida)
	return array_lineas_intervalos              #Devolvemos una lista. Cada elemento de la lista es una lista con 3 elementos del tipo [chr10, 0, 20]

def obtener_intervalo_mid_ant_s(mid_ant_s, fichero_bed):
	fh = open(fichero_bed, 'r')
        linea = (list(fh)[mid_ant_s]).strip().split()
        mid_ant_s_intervalo_inicio = int(linea[1])
        mid_ant_s_intervalo_final = int(linea[2])
        fh.close()
        return mid_ant_s_intervalo_inicio, mid_ant_s_intervalo_final

def comprobar_salto(fichero_pileup):

	"""Recorremos el fichero *.pileup y comprobamos  si hay algun salto superior a 1
	En el caso de que esto ocurra se llama a la funcion actual_en_intervalo() para comprobar
	dentro de que rango se encuentra.
	
	Se tratan como casos especiales la primera y ultima linea.
	"""

	num_chr = obtener_chr(sys.argv[2])

	output = "%s_coordenadas_no_cubiertas" %(num_chr)
	str(output)

	fh = open(output, 'w') 		#abrimos el fichero para borrarlo en caso de que ya exista
	fh.close()

	ant, ant_final = obtener_primera_linea(sys.argv[2])
	ult_intervalo_inicio, ult_intervalo_final = obtener_ultimo_bed(sys.argv[2])
	ult_pileup = obtener_ultimo_pileup(sys.argv[1])

	total_lineas = obtener_num_lineas_totales_bed(sys.argv[2])
	array_intervalos = obtener_array_intervalos(sys.argv[2])
	
	mid = 0
	mid_ant_s_old = 0
	
	primera_linea = True
	ultima_linea = False
	primera_iteracion = True
	fo = open(output, 'a')

	for linea in open (fichero_pileup):
		linea_completa = linea.strip().split()
		actual = int(linea_completa[1])

		if (primera_linea and actual >= ant):
			mid_ant_s = 0
                        res_busqueda, ant_intervalo_inicio, ant_intervalo_final, intervalo_inicio, intervalo_final, mid = busqueda_en_fichero_bed(array_intervalos, actual, 0, total_lineas)

			if (mid - mid_ant_s >= 1):
				intervalos_intermedios = array_intervalos[mid_ant_s:mid]

				for intervalo in intervalos_intermedios:
					fo.write("%s_%s_%s_+\t%s\t%s\t+\n" %(num_chr, intervalo[1], intervalo[2], intervalo[1], intervalo[2]))

			if (actual > intervalo_inicio):
				fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, intervalo_inicio, actual - 1, intervalo_inicio, actual - 1))

			mid_ant_s = mid
			mid_ant = mid - 1  
			primera_linea = False
		
                if (actual - ant > 1 and primera_iteracion == False or actual == ult_pileup):     #caso normal, para cualquier linea != primera or ultima
                        res_busqueda, ant_intervalo_inicio, ant_intervalo_final, intervalo_inicio, intervalo_final, mid = busqueda_en_fichero_bed(array_intervalos, actual, 0, total_lineas)
                        mid_ant_s_intervalo_inicio, mid_ant_s_intervalo_final = obtener_intervalo_mid_ant_s(mid_ant_s, sys.argv[2])
                        if (ant	 < mid_ant_s_intervalo_final and mid - mid_ant_s > 1):	#en lineas colindantes
                            fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, ant + 1, mid_ant_s_intervalo_final, ant + 1, mid_ant_s_intervalo_final))
                        if (ant < ant_intervalo_final and mid - mid_ant_s == 1):	#en lineas colindantes
                            fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, ant + 1, ant_intervalo_final, ant + 1, ant_intervalo_final))
                        if (mid_ant == mid_ant_s and mid - mid_ant_s == 1 and intervalo_inicio < actual < intervalo_final):		#caso que hayan dos saltos consecutivos
                            fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, intervalo_inicio, actual - 1, intervalo_inicio, actual - 1))
                        if (mid == mid_ant_s and actual != ult_pileup):                           #en la misma linea
                            fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, ant + 1, actual - 1, ant + 1, actual - 1))
                        if (mid == mid_ant_s and actual == ult_pileup and actual - ant > 1):                           #varios saltos ult linea
							fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, ant + 1, actual - 1, ant + 1, actual - 1))

                        if (ant < ant_intervalo_final and mid - mid_ant_s > 1):	#caso en el que hay intervalos de por medio sin cubrir entre los 2 saltos
				if (actual != ult_pileup):
					intervalos_intermedios = array_intervalos[mid_ant_s+1:mid]
				else:
					intervalos_intermedios = array_intervalos[mid_ant_s+1:mid]
					mid_ant_s_intervalo_inicio, mid_ant_s_intervalo_final = obtener_intervalo_mid_ant_s(mid_ant_s, sys.argv[2])
					fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, ant + 1, mid_ant_s_intervalo_final, ant + 1, mid_ant_s_intervalo_final))
                                for intervalo in intervalos_intermedios: #escribimos los intervalos enteros no utilizados que se encuentran entre el ant y actual
                                        fo.write("%s_%s_%s_+\t%s\t%s\t+\n" %(num_chr, intervalo[1], intervalo[2], intervalo[1], intervalo[2]))
                        	if (actual > intervalo_inicio and actual < intervalo_final):
                                    fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, intervalo_inicio, actual - 1, intervalo_inicio, actual - 1))

			if (actual == ult_pileup and intervalo_inicio <= actual < intervalo_final):
				fo.write("%s_%d_%d_+\t%d\t%d\t+\n" %(num_chr, actual + 1, intervalo_final, actual + 1, intervalo_final))
				if (ult_intervalo_inicio <= actual <= ult_intervalo_final):
					fo.close()
					return output
				else:
					intervalos_finales = array_intervalos[mid+1:]
        	                        for intervalo in intervalos_finales: #escribimos los intervalos enteros no utilizados que se encuentran entre el ant y actual
        	        	                fo.write("%s_%s_%s_+\t%s\t%s\t+\n" %(num_chr, intervalo[1], intervalo[2], intervalo[1], intervalo[2]))
					fo.close()
					return output
			mid_ant_s_old = mid_ant_s
			mid_ant_s = mid
		mid_ant = mid
		ant = actual	#actualizamos ant en cada iteracion		
		primera_iteracion = False
	fo.close()
	return output

def busqueda_en_fichero_bed(A, value, low, high):	#Se pasa la lista, valor a buscar, y los limites inferior y superior de la lista
	while low <= high:				#loop hasta que el valor superior de la lsita es menor que el valor inferior (ejem. busqueda se finaliza)
		mid = int((low + high) / 2)             #determina el punto central de la lista
		if (int(A[mid][1]) > value):              #es este valor inferior al punto central? redefine los parametros de la busqued
			high = mid - 1

        	elif (int(A[mid][1]) < value):            #es el valor central inferiror al valor biuscado? Si es asi, se redefinen los parametros
			low = mid + 1
        	else:
			ant_intervalo_inicio = int(A[mid-1][1])
			ant_intervalo_final = int(A[mid-1][2])
		        intervalo_inicio = int(A[mid][1])
       			intervalo_final = int(A[mid][2])
			return True, ant_intervalo_inicio, ant_intervalo_final, intervalo_inicio, intervalo_final, mid

        intervalo_inicio = int(A[mid][1])
        intervalo_final = int(A[mid][2])

	if not (int(A[mid][1]) <= value <= int(A[mid][2])):   #hay que restarle si ya estamos en el medio
		mid = mid - 1

        ant_intervalo_inicio = int(A[mid-1][1])
        ant_intervalo_final = int(A[mid-1][2])
        intervalo_inicio = int(A[mid][1])
        intervalo_final = int(A[mid][2])

	return True, ant_intervalo_inicio, ant_intervalo_final, intervalo_inicio, intervalo_final, mid


if __name__ == '__main__':

	if len(sys.argv) == 1:	
		print "\nEste script obtiene la lista de las coordenadas no cubiertas a partir de los ficheros *.pileup y *.bed de un mismo cromosoma\n"
		print "Ej archivo *.pileup:"
		print "chr10   61993   T       T       30      0       37      1       .     B"
		print "chr10   61994   T       T       30      0       37      1       .     E"
		print "chr10   61995   A       A       30      0       37      1       .     S\n"
		print "Ej archivo *bed:"
		print "chr15   61900           61950\nchr10   61951           61972\nchr10   61973           62001\n"
		print "El uso correcto de este script es el siguiente:\nEjemplo:\tcoordenadas_no_cubiertas.py chr15.pileup chr15.bed\n"
		print "\t\t### A T E N C I O N ###\n"
		print "Asegurate que el archivo *.pileup y el archivo *.bed son del MISMO cromosoma\n"
		print "El output para el cromosoma15 con los archivos chr15.pileup y chr15.bed seria el siguiente:"
		print "(el archivo se genera en el directorio en el cual te encuentras)\n"
		print "chr15_coordenadas_no_cubiertas\n"
		sys.exit()

	if len(sys.argv) < 3:
		print "Necesito 2 argumentos, el uso correcto de este script es el siguiente:\ncoordenadas_no_cubiertas.py chrX.pileup chrX.bed"

	else:	#entramos al programa si el numero de argumentos proporcionados es correcto
		inicio = time.time()
		print "\n########## TRABAJANDO #########\n"
                print "Los ficheros input proporcionados son:\t",sys.argv[1],"\t",sys.argv[2],"\n"
		output = comprobar_salto(sys.argv[1])
		final = time.time()
		duracion = final - inicio
                print "El fichero output que se ha generado es: ---->",output,"<----\n"
                print "Este resultado se ha generado en",duracion,"segundos\n"
                print "########## TERMINADO  #########\n"
	



