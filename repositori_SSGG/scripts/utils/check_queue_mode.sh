#!/bin/bash
echo "For changing queue mode, visit http://jenkins.sistemasgenomicos.com:8080"
echo
if [ `qconf -sq vaishia.q | grep slots | grep compute- | wc -l` == "1" ]; then 
	echo "INFO: currently on Production queues mode"; 
else 
	echo "INFO: currently on BioScope queues mode"; 
fi
