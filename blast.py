#!/usr/bin/env python

from Bio import Seq
from Bio import SeqIO
from subprocess import Popen, PIPE
import sys

def GenBank2Fasta(genbank_file,output="genbank.fasta"):
	#Esta función convierte el archivo genbank en un multifasta.
	my_output = open(output,"w")
	with open(genbank_file, "r") as input_handle:
		for record in SeqIO.parse(input_handle, "genbank"):
			for feature in record.features:
				if feature.type == 'CDS':
					#Algunas proteías no tiene translation y dan error.
					try: 
						locus = feature.qualifiers["locus_tag"][0]
						prot = feature.qualifiers["translation"][0]
						my_output.write(">"+locus+"\n"+prot+"\n\n")
					except:
						pass
	my_output.close()
	#Obtenemos el número de proteínas que tiene el archivo Genbank para sacarlo en pantalla.
	#El número de proteínas será el número de líneas entre tres, ya que la secuencia del archivo
	#es ID, secuencia, línea en blanco. Pongo división exacta // porque el split me da un elemento "" al final de más.
	my_input = open(output,"r")
	my_input = my_input.read()
	my_input = my_input.split("\n")
	return len(my_input)//3

def QueryCounter(query):
	#Devuelve el número de querys distintas.
	my_input = open(query,"r")
	querynames = []
	lines = my_input.read().split("\n")
	for line in lines:
		if line != "":
			if line[0] == ">":
				querynames.append(line[1:])
	return querynames

def Blaster(query,subject,output_file="Blaster_output_raw.tsv"):
	#Realiza el blast llamando al comando de linux.
	proceso = Popen(['blastp','-query',query,'-subject',subject,'-outfmt','6 qseqid sseqid pident qcovs evalue qstart qend'], stdout=PIPE, stderr=PIPE)
	error_encontrado = proceso.stderr.read().decode("utf-8")
	listado = proceso.stdout.read().decode("utf-8")

	proceso.stderr.close()
	proceso.stdout.close()

	my_output = open(output_file,"w")
	my_output.write(listado)
	my_output.close()

	P1 = open(output_file,"r")
	p1 = Popen(['wc','-l'], stdin=P1, stdout=PIPE)
	raw_hits = int(p1.stdout.read().decode("utf-8"))
	#Si no hay error devuelve un booleano True para continuar el análisis y el número de hits totales antes de filtrarlos.
	#Si falla se informa y se acaba el proceso.
	if error_encontrado:
		print("Se produjo el siguiente error:\n%s" % error_encontrado)
		return False,0
	return True,raw_hits
	
def BlastFilter(input_file="Blaster_output_raw.tsv",qcovs=0,pident=0,evalue=0.000001,output="Blaster_output.tsv"):
	#Esta función filtra los resultados de BLAST en función de los parámetros del usuario.
	my_input = open(input_file,"r")
	my_output = open(output,"w")
	#Ponemos la cabecera empezando en '#'
	my_output.write("#qseqid\tsseqid\tqcovs\tpident\tevalue\tqstart\tqend\n")
	#Los hits saldrán luego en pantalla y sirven para que:
	# - Si hits de una query = 1, informe de que no se puede hacer un árbol
	# - si hits de una query = 0, informe de que no se pueden analizar los dominios ni graficar.
	hits = 0
	for hit in my_input:
		split_hit = hit.split("\t")
		if float(split_hit[2]) >= qcovs and float(split_hit[3]) >= pident and float(split_hit[4]) < evalue:
			my_output.write(hit)
			hits +=1

	my_input.close()
	my_output.close()
	return hits

def QuerySpliter(query,input_file="Blaster_output.tsv"):
	#Esta función genera un archivo por query distinta y almacena en ellos los ID de las proteínas 
	#hits que han superado el filtro.
	my_input = open(input_file,"r")
	split = my_input.read().split("\n")
	double_split = []
	for hit in split:
		double_split.append(hit.split("\t"))
	my_input.close()
	#Llamo a los archivos con '€€€' para eliminarlos después fácilmente cuando tengamos los multifasta
	querys = QueryCounter(query)
	IDquerys = []
	for query in querys:
		IDquerys.append("€€€"+query)
		my_output = open("€€€"+query,"w")
		my_output.write("#Proteínas para %s\n" %query)	
		for hit in double_split:
			if hit[0] == query:
				my_output.write(hit[1]+"\n")
	#Además nos devuelve el nombre de las querys, lo que nos servirá a lo largo de todo el main.py para
	#llamar y crear archivos.
	return IDquerys

def GenbankFilter(input_file1,input_file2="genbank.fasta"):
	#Esta función usa los archivos de cada query (input_file1) con los ID de las proteínas hits para buscar dichos ID
	#en el archivo multifasta que hemos generado con Genbank2Fasta() y generar multifastas pequeños con las secuencias
	#asociadas a los ID.
	my_input = open(input_file1,"r")
	distinct = []
	for hit in my_input:
		if hit[0]!="#":
			if hit not in distinct:
				distinct.append(hit[:-1])
	my_input.close()
	
	my_input = open(input_file2,"r")
	my_output = open(input_file1[3:],"w")

	gbsplit = my_input.read().split("\n")
	for ID in distinct:
		for i in range(len(gbsplit)//3):
			if gbsplit[3*i] == ">"+ID:
				my_output.write(gbsplit[3*i]+"\n"+gbsplit[3*i+1]+"\n\n")
	
	my_output.close()
	#Se devuelven los hits de cad query para sacarlos en pantalla.
	return len(distinct)