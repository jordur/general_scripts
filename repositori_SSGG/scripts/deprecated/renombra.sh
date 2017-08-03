#!/bin/bash
#
# renombra.sh 0.4.1 - GPL
# (c) Guimi
# Script para el renombrado masivo de ficheros
# http://guimi.net
#
# Ultima modificacion: Guimi 2006-01
#

####################################
#### ERRORES
E_NOPARAM=64        # Error en el numero de parametros
E_PARAMINVAL=65     # Error en los parametros
E_SALINVAL=66       # Los parametros generan una salida invalida


####################################
#### AYUDA

# Funciones que muestran la ayuda
muestra_ayuda ()
{
  echo "Use `basename $0` --help  para obtener ayuda" >&2
}

muestra_ayuda_larga ()
{
  more << EOF

renombra.sh
(c) 2005 Guimi (http://guimi.net)
Este software se publica bajo la licencia GNU Public License (GPL) v.2
This software is under the GNU Public License (GPL) v.2

Uso:
----
renombra.sh 'filtro' [prefijo] [opciones]
  Renombra los ficheros del directorio incluyendo opcionalmente prefijos,
  sufijos y contadores.
  Se puede filtrar los ficheros que deseamos renombrar. Ver ejemplos al final.

mrename 'filtro' prefijo -c | -m
  mrename esta en desuso, utilice renombra.sh
  El uso de mrename implica los parametros -cn 1 -d 1 -v -i


Opciones:
---------
  -c
  --copy
     Copia los ficheros en vez de renombrarlos
  -cn [n]
  --contador [n]
     Utiliza un contador, iniciado en 'n' o en '0'
  -cs
  --contador_sufijo
     Situa el contador como sufijo
     Implica --contador
  -d n
  --digitos n
     Establece el numero minimo de digitos del contador. Por omision 3
     Implica --contador
  -e
  --extensiones
     No modifica las extensiones -parte posterior al primer punto-
  -f
  --force
     No pide confirmacion en caso de que el fichero ya exista
  -h
  --help
     Muestra esta ayuda
  -i
  --interactive
     Pide confirmacion en caso de que el fichero ya exista
  -m
  --move
     Mueve (renombra) los ficheros. Es la opcion por omision
  -n
  --no_original
     Elimina el nombre original
  -p cadena
  --prefijo cadena
     Incluye 'cadena' como prefijo del nombre actual
  -s cadena
  --sufijo cadena
     Incluye 'cadena' como sufijo del nombre actual
  -sed 'expresion1' 'expresion2'
     El script sustituye 'expresion1' por 'expresion2' utilizando sed
  -t
  --test
     El script no realiza ninguna accion
     Implica -v
  -v
  --verbose
     Muestra informacion de las acciones del script


Error de Salida invalida:
-------------------------
  Determinadas combinaciones de parametros pueden hacer que el nombre
  obtenido para todos los ficheros sea ""
  El script no permite dichas combinaciones
  Por ejemplo
    renombra.sh '*' -n

Ejemplos:
---------

renombra.sh '*' -cn  -t
  Muestra como modificaria todos los nombres de los ficheros anteponiendo un
  contador, pero sin llegar a hacerlo

renombra.sh '*' cambiado_
  Modifica todos los nombres de los ficheros anteponiendo
  'cambiado_' en el nombre

renombra.sh '*.jpg' -s _cambiado -e
  Modifica todos los nombres de los ficheros .jpg posponiendo '_cambiado'
  en el nombre, respetando las extensiones

renombra.sh '*.mp3' -p Serrat_ -i -e -cn 108 -cs -d 4 -n
  Modifica todos los nombres de los ficheros '*.mp3' por 'Serrat_xxxx.mp3'
  El contador xxxx tiene un minimo de 4 digitos y comienza en el 108
  Pide confirmacion antes de sobreescribir ficheros en caso de que ya existan

renombra.sh '*' -sed ' ' ''
renombra.sh '*' -sed ' '
  (Equivalentes) Elimina todos los espacios de todos los nombres de archivo

renombra.sh '*.mpg' -sed 'pg$' 'peg'
  Modifica todos los nombres de los ficheros '*.mpg', sustituyendo
  la pg por peg solo a final de nombre (consulte man sed)

Guimi
http://guimi.net

EOF
}

muestra_ayuda_mrename ()
{
  echo "mrename is DEPRECATED, please use renombra.sh --help" >&2
  echo "  to know the new features" >&2
  echo "mrename esta anticuado, por favor use renombra.sh --help" >&2
  echo "  para conocer las nuevas opciones" >&2
}

####################################
####################################
#### COMPROBACION DE PARAMETROS

# Parametros iniciales
OPCIONES=""
COMANDO="mv"

# Mas parametros iniciales para mrename
if [ `basename $0` = "mrename" ]
then
  contador="1"
  DIGITOS="1"
  OPCIONES="-i"
  VERBOSE="1"
  muestra_ayuda_mrename
fi

# Comprobamos los parametros
while [ $# -gt 0 ]; do    # Mientras tengamos parametros
  case "$1" in
    -c|--copy) # Copia los ficheros en vez de moverlos
      COMANDO="cp"
      ;;
    -cn|--contador) # Iniciamos el contador
      # Comprobamos si hay otro parametro despues
      if [ -n "$2" ]
      then
        case $2 in
          -*)  # Si el siguiente parametro es una opcion ponemos el contador a 0
            contador='0'
            ;;
          *[!0-9]*) # Si algun caracter no es un digito generamos error
            echo "Error: Contador invalido" >&2
            muestra_ayuda
            exit $E_PARAMINVAL
            ;;
          *)  # Tomamos el contador
            contador="$2"
            shift
            ;;
        esac
      else # Si no hay otro parametro iniciamos el contador a 0
        contador='0'
      fi
      ;;
    -cs|--contador_sufijo) # Situa el contador como sufijo
      CONTADOR_SUFIJO="1"
      # Si no esta definido el contador lo ponemos a 0
      if [ -z "$contador" ]; then
        contador="0"
      fi
      ;;
    -d|--digitos) # Establece el numero de digitos del contador
      # Comprobamos si hay otro parametro despues
      if [ -z "$2" ]; then
        echo "Error: Numero de digitos invalido" >&2
        muestra_ayuda
        exit $E_NOPARAM
      else
        case $2 in
          -*)  # Si el siguiente parametro es una opcion generamos error
            echo "Error: Numero de digitos invalido" >&2
            muestra_ayuda
            exit $E_PARAMINVAL
            ;;
          *[!0-9]*) # Si algun caracter no es un digito generamos error
            echo "Error: Numero de digitos invalido" >&2
            muestra_ayuda
            exit $E_PARAMINVAL
            ;;
          *)  # Tomamos el numero de digitos
            DIGITOS="$2"
            # Si no esta definido el contador lo ponemos a 0
            if [ -z "$contador" ]; then
              contador="0"
            fi
            ;;
        esac
        shift
      fi
      ;;
    -e|--extensiones) # No modificamos las extensiones
      EXTENSIONES="1"
      ;;
    -f|--force) # Parametro para el comando
      OPCIONES="-f"
      ;;
    -h|--help) # Mostramos la ayuda completa
      muestra_ayuda_larga
      exit 0
      ;;
    -i|--interactive) # Parametro para el comando
      OPCIONES="-i"
      ;;
    -m|--move) # Mueve (renombra) los ficheros
      COMANDO="mv"
      ;;
    -n|--no_original) # Eliminamos el nombre original
      NO_NOMBRE="1"
      ;;
    -p|--prefijo) # Tomamos un prefijo
      if [ -z "$2" ]
      then
        echo "Error: Falta el prefijo" >&2
        muestra_ayuda
        exit $E_NOPARAM
      else
        cadena_prefijo="$2"
        shift
      fi
      ;;
    -s|--sufijo) # Tomamos un sufijo
      if [ -z "$2" ]
      then
        echo "Error: Falta el sufijo" >&2
        muestra_ayuda
        exit $E_NOPARAM
      else
        cadena_sufijo="$2"
        shift
      fi
      ;;
    -sed) # Traduce cadena1 en cadena2
      # Comprobamos si hay otro parametro despues
      #  cadena_tr_2 puede ser nula (cambiar cadena_tr_1 por nada)
      if [ -z "$2" ]; then
        echo "Error: Falta(n) cadena(s) del parametro tr" >&2
        muestra_ayuda
        exit $E_NOPARAM
      else
        cadena_tr_1=$2
        cadena_tr_2=$3
        shift
        shift
      fi
      ;;
    -t|--test) # No realizamos el cambio
      TEST="1"
      VERBOSE="1"
      ;;
    -v|--verbose) # Hacemos que el programa de informacion
      VERBOSE="1"
      ;;
    -*) # Parametro incorrecto
      echo "Error: Parametro invalido" >&2
      muestra_ayuda
      exit $E_PARAMINVAL
      ;;
    *) # Obtenemos la cadena de filtro
      if [ -z "$filtro" ]; then
        filtro=$1
      elif [ -z "$cadena_prefijo" ]; then
        cadena_prefijo="$1"
      else
        # Si ya tenemos una cadena y un prefijo mostramos error
        echo "Error: Numero de parametros invalido" >&2
        muestra_ayuda
        exit $E_NOPARAM
      fi
      ;;
  esac
  shift       # Pasamos al siguiente parametro
done

# Si usamos el script como mrename necesitamos un prefijo
if [ `basename $0` = "mrename" ] && [ -z "$cadena_prefijo" ]
then
  muestra_ayuda
  exit 0
fi

# Si no hemos indicado cadena_a_buscar
# generamos error
if [ -z "$filtro" ]
then
  ##filtro=""
  echo "Error: Falta por definir el filtro" >&2
  muestra_ayuda
  exit $E_NOPARAM
fi

# Iniciamos el numero de digitos si no lo esta ya
if [ -z "$DIGITOS" ]; then
  DIGITOS=3
fi
FORMATO='%0'$DIGITOS'd\n'

####################################
####################################
#### RENOMBRAMOS FICHEROS

# Comprobamos si hay ficheros a renombrar
ls $filtro >/dev/null 2> /dev/null
if [ "$?" -gt 0 ]
then
  if [ "$VERBOSE" = "1" ]; then
    echo "No hay ficheros que renombrar"
  fi
  exit 0
fi

# BUCLE DE RENOMBRADO de los ficheros
for fichero_actual in $filtro
do

# Las siguientes alternativas no funcionan correctamente
# con ficheros que contienen espacios
#
#> lista_ficheros=`ls $filtro`
#> for fichero_actual in $lista_ficheros
#
#> for fichero_actual in `ls $filtro`

  # Establecemos el contador con el formato pedido
  if [ -n "$contador" ]
  then
    nuevo_contador=`printf ${FORMATO} $contador`
    contador=$(($contador+1))
  else
    nuevo_contador=""
  fi

  # Comprobamos si hay que separar la extension
  if [ "$EXTENSIONES" = "1" ]
  then
    # Tomamos la posicion del punto
    posicion_punto=`expr index "$fichero_actual" \.`
    if [ $posicion_punto -gt 0 ] # Si hay punto
    then
      # Tomamos el nombre y la extension
      longitud_fichero=${#fichero_actual}
      cadena_extension=${fichero_actual:posicion_punto}
      longitud_extension=${#cadena_extension}
      longitud_nombre=$(($longitud_fichero-$longitud_extension-1))
      nombre_fichero=${fichero_actual:0:$longitud_nombre}
    else # Si NO hay punto
      # Tomamos el nombre y anulamos la extension anterior que hubiese
      nombre_fichero=$fichero_actual
      unset cadena_extension
    fi
  else # Si NO tratamos las extensiones
    # No separamos la extension
    nombre_fichero=$fichero_actual
  fi

  # Eliminamos el nombre original si procede
  if [ -n "$NO_NOMBRE" ]
  then
    nombre_fichero=""
  elif [ -n "$cadena_tr_1" ]
  # Realizamos la sustitucion sobre el nombre
  then
    nombre_fichero=$(echo $nombre_fichero | sed s/"${cadena_tr_1}"/"${cadena_tr_2}"/g)
  fi

  # Generamos el nuevo nombre con el contador donde corresponda
  if [ -n "$CONTADOR_SUFIJO" ]; then
    cadena_final=$cadena_prefijo$nombre_fichero$nuevo_contador$cadena_sufijo
  else
    cadena_final=$cadena_prefijo$nuevo_contador$nombre_fichero$cadena_sufijo
  fi

  # Anyadimos la extension al nuevo nombre
  if [ -n "$cadena_extension" ]; then
    cadena_final=$cadena_final'.'$cadena_extension
  fi

  # Determinadas combinaciones de parametros pueden hacer que cadena_final sea ""
  # El resultado final seria un solo fichero sin nombre (el ultimo)
  # Prevenimos dicha situacion
  if [ "$cadena_final" = "" ]
  then
    echo "Error: el comando solicitado no genera una cadena valida"
    muestra_ayuda
    exit $E_SALINVAL
  fi

  # Si procede, realizamos el comando
  if [ "$fichero_actual" != "$cadena_final" ]; then
    # Mostramos lo que vamos a hacer
    if [ "$VERBOSE" = "1" ]; then
      echo "$COMANDO $OPCIONES '$fichero_actual' '$cadena_final'"
    fi

    if [ -z "$TEST" ]; then
      $COMANDO $OPCIONES "$fichero_actual" "$cadena_final"
      SALIDA=$?

      if [ "$SALIDA" -gt 0 ]
      then
        echo "$? - Error en $COMANDO $OPCIONES '$fichero_actual' '$cadena_final'" >&2
      fi
    fi
  fi

done


####################################
# Terminamos
exit 0
