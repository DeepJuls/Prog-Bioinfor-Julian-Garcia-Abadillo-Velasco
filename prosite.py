#!/usr/bin/env python

import re
from Bio.ExPASy import Prosite,Prodoc

def PatternParser():
	#El output lo saca en .tsv para luego convertirlo en una matriz usando dos split, uno de retorno de carro y otro de tabulador.
	handle = open("PrositeDB/prosite.dat","r")
	records = Prosite.parse(handle)
	my_output = open("PrositeDB/prosite_patterns","w")
	my_output.write("#Name\tAccesion\tpattern\tDescription\n")

	for record in records:
		my_output.write(record.name+"\t"+record.pdoc+"\t"+record.pattern+"\t"+record.description+"\n")
	my_output.close()
	return

def DescriptionParser():
	#El output se saca en texto plano y las consideraciones para su manejo se tienen en cuenta en DescriptionInfo().
	handle = open("PrositeDB/prosite.txt","r")
	records = Prodoc.parse(handle)
	my_output = open("PrositeDB/prosite_description","w")
	for record in records:
		my_output.write(record.accession+"\n")
		my_output.write(record.text)
	my_output.close()
	return

def DescriptionInfo(PDOC):
	my_input = open("PrositeDB/prosite_description","r")
	archivo = my_input.read()
	write = False
	texto = ""
	
	for line in archivo.split("\n"):
		#Para no tener en cuenta líneas en blanco
		if line != "":
			if line == PDOC:
				#Si la línea coincide con el argumento de entrada, entonces querremos escribir las líneas posteriores...
				write = True
			elif line[0:4] == "PDOC":
				#... Hasta que aparezca otra PDOC distinta.
				write = False
			if write and line != PDOC:
				#Write controla que se estén escribiendo líneas correspondientes al PDOC del argumento...
				#... y que sean distintas a la propia línea de PDOC.
				texto += "\n"+line
	return texto

def Prosite2Python(PrositePattern):
	#Para convertir el modelo de patrón de prosite al de python.
	if PrositePattern != "":
		#Esta condición porque hay tuplas que tienen vacío el campo de patrón.
		#Quitamos el punto.
		pattern = PrositePattern[:-1]
		#Quitamos el separador '-'.
		pattern = pattern.replace("-","")
		#Cualquier aa
		pattern = pattern.replace("x",".")
		#Exclusiones
		pattern = pattern.replace("{","[^")
		pattern = pattern.replace("}","]")
		#Repeticiones
		pattern = pattern.replace("(","{")
		pattern = pattern.replace(")","}")
		#Restricciones de extremos Nt y Ct.
		pattern = pattern.replace("<","^")
		pattern = pattern.replace(">","$")
		return pattern
	#Si no hay un patrón asociado a una determinada tupla, se devuelve una string que no va a aparecer nunca,
	return "ÑÑÑ"

def OBL_MegaParser(nombre_querys,ID):
	#Esta función es prácticamente un script entero, pero es para no repetir el cuádruple bucle (querys,hits de cada query,patrones para cada hit, matchs de cada patrón)
	#Generamos una matriz con la información de patrones. Se encarga de recorrer los archivos mulitfasta y anotar los dominios de las proteínas.
	my_input = open("PrositeDB/prosite_patterns","r")
	lines = my_input.read().split("\n") 
	patterns = []
	for line in lines:
		if line!="":
			if line[0]!="#":
				patterns.append(line.split("\t"))
	my_input.close()

	#Recorremos los archivos.
	uniquepatterns = []
	for query in nombre_querys:
		#Este contador tendrá en cuenta los matchs de todos los hits de cada query.
		contador_global = 0
		my_input = open(query,"r")
		fasta = my_input.read().split("\n")
		my_input.close()
		#Se guardará la información en archivos con el nombre de la query.
		my_output = open("Results/"+ID+"/Prosite/"+query+".txt","w")
		#Recorremos los hits de cada archivo.
		#Fasta será una lista con la siguiente estructura:
		#Índices 3k serán cabeceras de secuencias en formato fasta.
		#índices 3k+1 contendrán las secuencias asociadas a dichas cabeceras.
		#índices 3k+2 contienen líneas en blanco para separar.
		for i in range(len(fasta)//3):
			my_output.write(("Se han encontrado los siguientes patrones para la proteína %s:\n" %fasta[3*i][1:]))
			my_output.write("%s\n\n" %fasta[3*i+1])
			#Recorremos los patrones para cada hit.
			for pattern in patterns:
				#Tengo que usar re.search porque el booleano de re.finditer siempre da True, aún no habiendo matchs.
				haymatch = re.search(Prosite2Python(pattern[2]),fasta[3*i+1])
				if haymatch:
					if not pattern[0] in uniquepatterns:
						uniquepatterns.append(pattern[0])
					#Si existe algún match con este patrón...
					#Se reinicia el contador del patrón.
					contador_parcial = 0
					#Se escriben las caraterísticas del patrón.
					my_output.write("Nombre: {}\nAccesión: {}\nDescripción: {}\nPatrón: {}\nPosiciones del match en la proteína:\n".format(pattern[0],pattern[1],pattern[3],pattern[2]))
					matchs = re.finditer(Prosite2Python(pattern[2]),fasta[3*i+1])
					#Se recorren todos los matchs.
					for match in matchs:
						contador_global += 1
						contador_parcial += 1
						#Indexamos las secuencias a partir del uno para que quede coherente con su longitud en los archivos de prosite.
						my_output.write("{}. [{}-{}]\t".format(contador_parcial,match.start()+1,match.end()))
						#Estos saltos de carro son para que queden bien los outputs de prosite.
						if contador_parcial % 5 == 0:
							my_output.write("\n")
					my_output.write("\n")
					if contador_parcial % 5 != 0:
						my_output.write("\n")
		print("\tSe han encontrado {} coincidencias.\t[{}]".format(contador_global,query))
	my_output.close()
	return patterns,uniquepatterns

def OPT_MegaParser(nombre_querys,ID):
	#Esta función saca las proteínas únicas y los dominios en el formato para el módulo graph.
	my_input = open("PrositeDB/prosite_patterns","r")
	lines = my_input.read().split("\n") 
	patterns = []
	for line in lines:
		if line!="":
			if line[0]!="#":
				patterns.append(line.split("\t"))
	my_input.close()

	#Recorremos los archivos.
	#Estas listas se devuelven en el Return para temas de graficación.
	uniqueprots = []
	lenprots = []
	for query in nombre_querys:
		my_input = open(query,"r")
		fasta = my_input.read().split("\n")
		my_input.close()
		#Recorremos los hits de cada archivo.
		#Fasta será una lista con la siguiente estructura:
		#Índices 3k serán cabeceras de secuencias en formato fasta.
		#índices 3k+1 contendrán las secuencias asociadas a dichas cabeceras.
		#índices 3k+2 contienen líneas en blanco para separar.
		for i in range(len(fasta)//3):
			if not fasta[3*i][1:] in uniqueprots:
				uniqueprots.append(fasta[3*i][1:])
				lenprots.append(len(fasta[3*i+1]))
				pd_output = open("Graph/Dominios/"+fasta[3*i][1:],"w")
				#La primera línea tiene la longitud de la proteína para hacer la escala cuando se grafique.
				pd_output.write("#"+str(len(fasta[3*i+1]))+"\n")
				#Recorremos los patrones para cada hit.
				for pattern in patterns:
					#Tengo que usar re.search porque el booleano de re.finditer siempre da True, aún no habiendo matchs.
					haymatch = re.search(Prosite2Python(pattern[2]),fasta[3*i+1])
					if haymatch:
						pd_output.write(pattern[0]+",")
						#Si existe algún match con este patrón...
						#Se escribe el nombre del dominio.
						matchs = re.finditer(Prosite2Python(pattern[2]),fasta[3*i+1])
						#Se recorren todos los matchs y se escriben en la misma linea del nombre del dominio, separado todo por comas.
						for match in matchs:
							#Se indexa desde uno y se guarda el punto central del dominio.
							pd_output.write(str((match.start()+1+match.end())/2)+",")
						pd_output.write("\n")
		pd_output.close()
	#Aunque la longitud de las proteínas las tengo en los archivos las saco en una lista para graficar BLAST.
	return uniqueprots,lenprots

def PatternInformer(uniquepatterns,patternsDB,ID):
	#Esta función escribe el archivo que contiene la información extendida (prosite.txt) de los dominios detectados al menos una vez. 
	my_output = open("Results/"+ID+"/Prosite/_DomainInformation.txt","w")
	my_output.write("##### INFORMACIÓN ADICIONAL DE TODOS LOS PATRONES DETECTADOS #####\n")
	for pattern in uniquepatterns:
		for entry in patternsDB:
			if pattern == entry[0]:
				my_output.write("\n{} [{}]{}\n".format(entry[0],entry[1],DescriptionInfo(entry[1])))
	return len(uniquepatterns)

