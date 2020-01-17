#!/usr/bin/env python

import sys
import os

from subprocess import call, Popen, PIPE
import matplotlib.pyplot as plt

from time import time
from datetime import date,timedelta

import intro
import blast
import muscle
import prosite
import graph

### PRESENTANCIÓN ###
#Este try abarca todo el script por si el usuario interrumpe manualmente que quede bonito.
try:
	call("clear")
	print("Bienvenido a GVelProt.\nSiga las intrucciones que se le indican a continuación.\n")
	#Llamamos a la ayuda.
	intro.Help()

### PEDIR DATOS ###
	#Recogemos los parámetros en el diccionario data.
	data = intro.Data()
	#Guardamos en ID el nombre del folder donde se debe almacenar todo.
	ID = intro.IDproceso()
	#Sacamos la fecha para el log.
	date = date.today().strftime("%d-%b-%Y")
	#Abro el log para ir guardando los avances.
	log = open("log.txt","w")
	log.write("[ID Proceso]\n\t{}\n\n[Fecha]\n\t{}\n\n[Parámetros]\n".format(ID[8:],date))
	for param in list(data.keys()):
		log.write("\t - {}: {}\n".format(param,data[param][0]))

	call("clear")
	print("ID del proceso: %s\n\n[PROGRESO DEL ANÁLISIS] - Blast & Muscle\n" %ID[8:])
	print("Parseando archivo genbank...")

### BLAST & MUSCLE ###
	log.write("\n[PROGRESO]\n")
	Ti_blast = time()
	try:
		#Intentamos hacer el parseo de Genbank, si no es un GB fallará y se informará.
		number_prots = blast.GenBank2Fasta(data["subject"][0],output="genbank.fasta")

		if number_prots == 0:
			print("El archivo genbank no ha podido convertirse a fasta, revise dicho input.")
			sys.exit()
		else:
			#Este bloque hace el BLAST, filtra, y muestra iterativamente el número de hits filtrados que tendrá cada query.
			log.write("\t- Multifasta obtenido con éxito desde archivo genbank: %d proteínas obtenidas.\n" %number_prots)
			print("\tArchivo genbank convertido a fasta con éxito:\n\t%d proteínas contenidas en total.\n" %number_prots)
			querys = blast.QueryCounter(data["query"][0])
			print("Parseando archivo query...\n\tSe han detectado %d querys distintas\n" %len(querys))
			print("Ejecutando Blast...")
			try_blast,raw_hits = blast.Blaster(data["query"][0],"genbank.fasta")
			print("\tSe han obtenido %d hits antes de aplicar el filtrado. \n" % raw_hits)
		if try_blast:
			#Si el booleano que se obtiene de BLAST es True hace el filtrado según los parámetros del usuario.
			print("Aplicando el filtrado...")
			blast.BlastFilter(qcovs=data["qcovs"][0],pident=data["pident"][0],evalue=data["evalue"][0])
			IDquerys = blast.QuerySpliter(data["query"][0])
			#El diccionario hits_totales nos informará de qué análisis podemos hacer con cada query en función del
			#número de hits que nos aporte.
			hits_totales = []
			nombre_querys = []
			for query in IDquerys:
				hits = blast.GenbankFilter(query)
				hits_totales.append(hits)
				nombre_querys.append(query[3:])
				print("\tSe han conservado {} hits tras el filtrado. [{}]".format(hits,query[3:]))
			print("")
			hits_totales = dict(zip(nombre_querys,hits_totales))
			#Se crean las carpetas específicas de este análisis: GVelProt00034 por ejemplo.
			call(["mkdir","Data/"+ID])
			call(["mkdir","Results/"+ID])
			#Ya podemos guardar el log en esta carpeta para no tenerla en medio.
			call(["mv","log.txt","Results/"+ID])
			#Creamos la carpeta específica para los árboles antes de correr MUSCLE.
			call(["mkdir","Results/"+ID+"/Trees"])

			print("Ejecutando Muscle...")
			#Contamos el número de árboles para ver que mensaje sacamos en pantalla.
			num_arboles = 0
			for query in nombre_querys:
				call(["rm","€€€"+query])
				if hits_totales[query] < 2:
					print("\tNo hay suficientes secuencias para generar el árbol. [%s]" %query)		
				else:
					num_arboles += 1
					print("\tAlineando las secuencias...\t[%s]" %query)
					muscle.Aligmenter(query)
					print("\tGenerando el árbol...\t\t[%s]" %query)
					muscle.TreeMaker(query)
					call(["mv",query+".tre","Results/"+ID+"/Trees"])
			
			Tt_blast = time() - Ti_blast
			#Si num_arboles es distinto de cero se evalúa si hay árboles para todas las querys o no y se escribe una cosa u otra en el log.
			if bool(num_arboles):
				if num_arboles == len(nombre_querys):
					log.write("\t- Se ha generado un árbol filogenético para todas las querys.\n")
				elif num_arboles < len(nombre_querys):
					log.write("\t- Se ha generado un árbol filogenético para {} de las {} querys.\n".format(num_arboles,len(nombre_querys)))
				print("\nPROCESO TERMINADO [{}]. Compruebe el directorio Results/{}/Trees".format(str(timedelta(seconds = int(Tt_blast))),ID))
			else:
				log.write("\t- No se ha podido generar ningún árbol filogenético: todas las querys tenían menos de dos hits.\n")
			#Movemos los archivos intermedios a su carpeta correspondiente en Data, a la espera de ser eliminados o no...
			call(["mv","Blaster_output_raw.tsv","genbank.fasta","Data/"+ID])
			#Menos el arhcivo con los resultados de BLAST filtrados que nos necesito para graficar.
			call(["cp","Blaster_output.tsv","Data/"+ID]) 

		else:
			print("Se ha producido un error, revise los archivos y parámetros de entrada.")
			log.write("\tSe ha producido un error.\n")
			log.close()
			sys.exit()

	except FileNotFoundError:
		#Excepción relativa a la no existencia de los archivos de entrada.
		print("Se ha producido un error, revise los archivos y parámetros de entrada.")
		log.write("\tSe ha producido un error.\n")
		log.close()
		sys.exit()

### PROSITE ###
	#En primer lugar guardo las querys cuyos hits son mayores de cero, ya que son con las que podemos trabajar.
	#Si no hay ninguno se acaba el análisis.
	nombre_querys_prosite = []
	for query in nombre_querys:
		if hits_totales[query] != 0:
			nombre_querys_prosite.append(query)

	if len(nombre_querys_prosite) == len(nombre_querys):
		log.write("\t- Se ha realizado un análisis de dominios en los hits de todas las querys.\n")
	elif bool(len(nombre_querys_prosite)):
		log.write("\t- Se ha realizado un análisis de dominios en los hits de {} de las {} querys.\n".format(len(nombre_querys_prosite),len(nombre_querys)))
	else:
		log.write("\t- No se ha podido realizar un análisis de dominios debido a la ausencia de hits.\n")
	if len(nombre_querys_prosite) == 0:
		print("\nNo hay ninguna proteína. No se puede continuar al servicio PROSITER.")
		sys.exit()

	input("\nPulse INTRO para acceder a la aplicación PROSITE para la búsqueda de dominios. >>  ")
	Ti_prosite = time()
	print("\n[PROGRESO DEL ANÁLISIS] - Prosite\n")
	#Creamos la carpeta para almacenar los resultados de PROSITE.
	call(["mkdir","Results/"+ID+"/Prosite"])

	#Estas condiciones evaluan si existe o no los archivos parseados que contienen la información de .txt y .dat
	#ya que es un gasto de recursos parsearlo cada vez siendo la misma base de datos.
	if not os.path.exists("PrositeDB/prosite_patterns"):
		prosite.PatternParser()
		log.write("\t- Se ha generado un archivo con la información general de los dominios.\n")
	if not os.path.exists("PrositeDB/prosite_description"):
		prosite.DescriptionParser()
		log.write("\t- Se ha generado un archivo con la información específica de los dominios.\n")
	#Borramos una hipotética carpeta /Graph que no debe existir porque la voy a borrar al final, pero por si acaso
	#hay un error y se queda el proceso a medias.
	call(["rm","-r","-f","Graph"])
	#Creamos la carpeta temporal /Graph que va a almacenar los archivos intermedios para la graficación.
	call(["mkdir","Graph"])
	call(["mv","Blaster_output.tsv","Graph"])
	call(["mkdir","Graph/Dominios"])
	print("Buscando dominios de la base de datos...")
	#Parseamos las querys y prosite.dat para la parte obligatoria y sacamos esas dos variables para escribir
	#la información de prosite.dat una única vez por dominio encontrado.
	patternsDB,uniquepatterns = prosite.OBL_MegaParser(nombre_querys_prosite,ID)
	print("\nRecopilando información sobre los dominios encontrados...")
	#PatternInformer() escribe la información de prosite.dat una vez por dominio. Dominios es el número de dominios distintos.
	dominios = prosite.PatternInformer(uniquepatterns,patternsDB,ID)
	#Este parser vuelve a recorrer las querys y los patrones, pero nos saca los archivos necesarios para graficar las proteías
	#con sus correspondientes dominios a lo largo de la secuencia.
	uniqueprots,lenprots = prosite.OPT_MegaParser(nombre_querys_prosite,ID)
	print("\tSe han encontrado {} dominios distintos en {} proteínas.\n".format(dominios,len(uniqueprots)))
	Tt_prosite = time() - Ti_prosite 
	print("PROCESO TERMINADO [{}]. Compruebe el directorio Results/{}/Prosite".format(str(timedelta(seconds = int(Tt_prosite))),ID))

	call(["mkdir","Data/"+ID+"/FastaFiles"])
	for file in nombre_querys:
		call(["mv",file,"Data/"+ID+"/FastaFiles"])

	if not data["keep_data"][0]:
		call(["rm","-r","Data/"+ID]) 

### GRAFICACIÓN ###
	input("\nPulse INTRO para continuar con la graficación de los resultados >>\n  ")

	Ti_graf = time()
	#Creamos la carpeta donde se almacenarán los archivos temporales de BLAST antes de graficar.
	call(["mkdir","Graph/Blast"])
	#Esta función genera los csv necesarios para la grafiación de BLAST.
	graph.Blast2Plotter(uniqueprots,lenprots)
	#Creamos las carpetas de resultados dodnde se van a almacenar los gráficos de BLAST.
	call(["mkdir","Results/"+ID+"/Graphics"])
	call(["mkdir","Results/"+ID+"/Graphics/BLAST"])
	print("Graficando resultados de BLAST...")

	lista_blast_global = []
	#Recorremos las proteínas que aparecen al menos una vez y graficamos sus resultados de BLAST asociados.
	for prot in uniqueprots:
		path = "Graph/Blast/"+prot
		lista_blast_global += graph.CSV2plt_BLAST(path,data["pident"][0])
		graph.BlastPlotter(path,data["pident"][0],ID)
		print("\tImagen PNG generada con éxito [%s]" %prot)
		plt.close()
	#Ahora graficamos todos los resultados en conjunto para verlos con perspectiva.
	graph.BlastPlotter(None,data["pident"][0],ID,Global=True,array=graph.GlobalBlaster(lista_blast_global,uniqueprots))
	log.write("\t- Se han graficado los resultados de BLAST.\n")
	print("\tImagen PNG generada con éxito [Resultado global]")
	#Ahora creamos la carpeta para guardar los resultados de PPROSITE.
	call(["mkdir","Results/"+ID+"/Graphics/PROSITE"])
	print("\nGraficando dominios...")
	#Graficación de los dominios de las proteínas que aparecen al menos una vez.
	for prot in uniqueprots:
		path = "Graph/Dominios/"+prot
		graph.DomainPlotter(path,ID)
		plt.close()
		print("\tImagen PNG generada con éxito [%s]" %prot)
	log.write("\t- Se han graficado los resultados de PROSITE.\n")
	Tt_graf = time() - Ti_graf

	print("\nPROCESO TERMINADO [{}]. Compruebe el directorio Results/{}/Graphics".format(str(timedelta(seconds = int(Tt_graf))),ID))
	#Sumamos todos los intervalos de tiempo para obtener el total.
	Tt = str(timedelta(seconds = int(Tt_blast + Tt_prosite + Tt_graf)))
	#Borramos la carpeta /Graph de archivos csv intermedios.
	call(["rm","-r","Graph"])
	log.write("\t- Se han eliminado los ficheros intermedios de graficación y se han cerrado todos los procesos iniciados.")
	log.write("")
	log.close()
	print("\nANÁLISIS COMPLETADO [%s]." %Tt)
except KeyboardInterrupt:
	print("\nSe ha interrumpido el proceso.")
	#Borramos forzadamente (puede estar o no) los archivos intermedios que hayan podido quedar en la carpeta principl.
	call(["rm","-f","-r","Graph","Blaster_output_raw.tsv","genbank.fasta","aligment","Blaster_output.tsv"])
	try: #Puede que no hayamos llegado a crear esto.
		for query in nombre_querys:
			call(["rm","-f","€€€"+query,query])
	except:
		pass
print("GVelProt 2020 - Desarrollado por Julián García-Abadillo Velasco")

