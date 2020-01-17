#!/usr/bin/env python

import os
from subprocess import call, Popen, PIPE

def Data():
	#Esta función parsea todos los parámetros que se necesitan, un argv.parser casero pero ya lo tenía hecho.
	print("A continuación se le pedirá que seleccione algunos parámetros.\nSi desea mantener los valores por defecto, cuando sea posible, pulse la tecla intro.\n")
	print("Estos son los archivos en el directorio: ")
	#Muestra en pantalla todos los archivos que hay para que el usuario no se tenga que acordar de cuales son los nombres exactos
	#de query y genbank. Se excluyen los propios módulos (*.py) y los directorios.
	p1 = Popen(["ls","-p"], stdout=PIPE)
	p2 = Popen(["grep","-v","/"], stdin=p1.stdout,stdout=PIPE)
	no_files = p2.stdout.read().decode("utf-8")
	
	for ls in no_files.split("\n"):
		if not ".py" in ls:
			print(ls)
	
	#Este diccionario almacenará todos los parámetros. En la lista que tienen por value, el primer elemento es el parámetro por defecto
	#o, si es cambiado, el que nos da el usuario. El segundo elemento es el tipo de dato que tiene que ser (está controlado ese error)
	#El tercero es la información que debe mostrar el input() para pedir el parámetro.
	questions = ["query","subject","qcovs","pident","evalue","keep_data"]
	answers = [
	      [None,str,"Introduzca el nombre del archivo .fasta que contiene la(s) query(s):"],
	      [None,str,"Introduzca el archivo genbank:"],
	      [0,float,"Introduzca el valor mínimo de coverage (%):"],
	      [0,float,"Introduzca el valor mínimo de identidad (%):"],
	      [0.000001,float,"Introduzca el valor máximo de e-value:"],
	      [False,bool,"¿Desea conservar los archivos intermedios? [En caso afirmativo escriba 'True']:"]
	      ]
	data = dict(zip(questions,answers))
	
	def subData(question):
		#Esta subfunción es el parser del diccionario, para recopilar los datos.
		while True:
			#Bucle infinito del que solo se sale si el usuario mete algo que sea válido.
			print ("\n"+data[question][2])
					
			if data[question][0] != None:
				print("El valor por defecto es %s" %data[question][0])
			answer = input("Nuevo valor para %s >> " %question)
					
			if answer == "":
				if data[question][0] != None:
					return
				else:
					print ("Es necesario introducir un valor.\n")

			else:
				try: 
					#Intentamos convertir el dato al tipo de dato que tiene asignado en el diccionario.
					data[question][0] = data[question][1](answer)
					return
				except:
					print ("Error en el tipo de dato, inténtelo de nuevo.\n")

	for question in questions:
	 	subData(question)
	#Devolvemos el diccionario para utilizar durante todo el main.
	return data

def IDproceso():
	#Esta función controla por una parte la creación de las carpetas principales Data y Results
	#y por otra me da el número de proceso que corresponde en función del número máximo que exista.
	if not os.path.exists("Data"):
		call(["mkdir","Data"])

	if not os.path.exists("Results"):
		call(["mkdir","Results"])
		return "GVelProt00001"
	else:
		proceso = Popen(["ls","Results"],stdout=PIPE)
		ls = proceso.stdout.read().decode("utf-8")
		ls = ls.split("\n")[:-1]
		IDs = []
		for ID in ls:
			IDs.append(int(ID[8:]))
		nuevoID = str(max(IDs)+1)
		
		while len(nuevoID)<5:
			nuevoID = "0" + nuevoID
		return "GVelProt" + nuevoID

def Help():
	#Solo se muestra si se especifica en el input "y" o "Y"
	ayuda = input("¿Desea acceder a la ayuda? Y/[N]>> ")
	if ayuda.upper() == "Y":
		call("clear")
		print("\nGVelProt - CÓDIGO LIBRE - Desarrollado por Julián García-Abadillo Velasco. 2020")
		print("\n[INTRODUCCIÓN DE DATOS]\n")
		print("\t* El programa solicitará al usuario seis parámetros:")
		print("\t\t - Query: archivo de tipo fasta con una o más secuencias que será lanzada contra el subject.")
		print("\t\t - Subject: archivo extraído de genbank que será convertido a un archivo multifasta y se usará como subject.")
		print("\t\t - Qcovs: valor de corte para el porcentaje de cobertura de BLAST, su valor por defecto es 0%. [Filtro desactivado]")
		print("\t\t - Pident: valor de corte para el porcentaje de identidad de BLAST, su valor por defecto es 0%. [Filtro desactivado]")
		print("\t\t - Evalue: valor de corte para e-value en los resultados de BLAST, su valor por defecto será 1e-6.")
		print("\t\t - keep_data: booleano que determina si borrar los archivos intermedios, almacenados en /Data.")
		print("\t* Si el nombre de los archivos está mal escrito o no existe, el programa terminará.")
		print("\n[BLAST & MUSCLE]\n")
		print("\t* En primer lugar, se parseará el archivo de genbank y se informará de las proteínas que contiene.")
		print("\t* El programa realizará un BLAST y filtrará los resultados para cada una de las querys en archivos multifasta con las secuencias hit.")
		print("\t* En el terminal se irá informando del avance del proceso.\n")
		print("\t* A continuación se alinea y genera el árbol de cada una de las querys. Si no existen al menos dos secuencias para alinear, no se continua con esa query. ")
		print("\t* Los resultados se almacenan en /Trees")
		input("\n[1/4]")
		call("clear")
		print("\nGVelProt - CÓDIGO LIBRE - Desarrollado por Julián García-Abadillo Velasco. 2020")
		print("\n[PROSITE]\n")
		print("\t* Se recorrerá la base de datos PROSITE para obtener los dominios que se encuentran en cada una de las secuencias hit de cada query.")
		print("\t* Para cada proteína, se dará información sobre el nombre, el código, el patrón aminoacídico y una breve descripción, así como la posición de todos los cortes.")
		print("\t* Encontrará esta información en la carpeta /Prosite\n")
		print("\t* Adicionalmente, la carpeta prosite contendrá un archivo llamado '_DomainInformation'\n\t  con más información sobre cada uno de los dominios encontrados al menos una vez en todo el análisis.\n")
		print("\t* La información por terminal corresponde al número total de coincidencias para cada hit,\n\t  por lo que un mismo patrón puede aparecer más de una y en más de un hit.")
		print("\t* Finalmente se mostrará el número de dominios diferentes que se han encontrado para el número de proteínas diferentes en todo el análisis.")
		print("\n[GRAFICACIÓN]\n")
		print("\t* Se graficará cada uno de los hits obtenidos en BLAST. Las imágenes generadas se almacenarán en /Graphics/BLAST .")
		print("\t* Se graficará cada una de las proteínas que ha aparecido al menos una vez como hit en el análisis,\n\t  con sus correspondientes dominios.")
		print("\t* Las imágenes generadas se almacenarán en /Graphics/PROSITE")
		input("\n[2/4]")
		call("clear")	
		print("\nGVelProt - CÓDIGO LIBRE - Desarrollado por Julián García-Abadillo Velasco. 2020")
		print("\n[DATA]\n")
		print("\t* Este directorio contendrá los archivos intermedios que se han necesitado para realizar los análisis, se incluye el archivo BLAST sin filtrar,")
		print("\t  El multifasta original, los multifasta secundarios de cada query...")
		print("\t* Para mantener esta carpeta debe escribir como parámetro para 'keep_data' la opción 'True', en cualquier otra cosa se eliminará el directorio.")
		print("\t* Tenga en cuenta que las carpetas /Results y /Data se crean autómáticamente,\n\t  la información de cada análisis se almacena en subcarpetas con su ID correspondiente.")
		print("\t* Por favor, no elimine los archivos con extensión .txt y .dat de la carpeta /PrositeDB, son necesarios para el análisis de Prosite.")
		print("\n[ALGUNAS CONSIDERACIONES]\n")
		print("\t*Para el correcto funcionamiento del programa se pide que se siga el siguiente formato:")
		print("\t\t - Los nombres de las querys no deben contener el símbolo '|' ni espacios en medio, ya que pueden confundir a python\n\t\t   al abrir archivos y escribir sobre ellos.")
		print("\t\t - El parser de Genbank solo reconoce archivos cuya extensión es .gbff, no acepta .gb.")
		print("\t\t - No interrumpir manualmente durante el proceso de graficación.")
		input("\n\n\n\n[3/4]")
		call("clear")
		print("\nGVelProt - CÓDIGO LIBRE - Desarrollado por Julián García-Abadillo Velasco. 2020")
		print("\n[NOTA DEL DESARROLLADOR]\n")
		print("\t* En la graficación de BLAST se grafica cada una de las proteínas y las querys que han matcheado con ella, el color indica\n\t   el porcentaje de identidad y la posición determina la cobertura.")
		print("\t* No he conseguido ponerle etiqueta al primer elemento de cada gráfico de BLAST, pero he visto que es un error recurrente en StackOverflow.")
		print("\t* Si ha llegado hasta aquí solo me queda darle las gracias y pedirle disculpas por el spanglish utilizado en el código\n\t  y la nomenclatura de los directorios y ficheros.")
		input("\n\n\n\n\n\n\n\n\n\n\n\n\n[4/4]")
	call("clear")	
	return
