#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np 
from random import randint

def CSV2plt_BLAST(file,min_pident):
	#Esta función convierte los archivos csv que salen de Blast2Plotter() que a su vez vienen del output de blast
	#en listas para el graficador BlastPlotter().
	my_input = open(file,"r")
	lines = my_input.read()
	lines = lines.split("\n")[:-1]
	nbars = len(lines)
	vector = []
	for i in range(nbars):
		if lines[i][0] == "#":
			#El primer elemento determina la posición Y en el gráfico, el segundo la longitud de la proteína,
			#el tercero donde empieza(superposición con hit), el cuarto es el nombre de la query y el último es el color.
			vector.append([0,int(lines[i][1:]),1,file[12:],100])
		else:
			tupla = lines[i].split(",")
			#Al nombre le cambio "_" por saltos de carro para que quede mejor en el plot.
			vector.append([-i,int(tupla[3])-int(tupla[2])+1,int(tupla[2]),tupla[0].replace("_","\n"),int(tupla[1])])
	#Una vez que tenemos la matriz de información hay que transponerla para que la reconozca plt.barh, que es el que tiene que graficar.
	vector = np.array(vector)
	vector = np.transpose(vector)
	vector = vector.tolist()
	my_input.close()
	return vector

def GlobalBlaster(recopilacion,uniqueprots):
	#A esta función le entra la suma de todas las salidas del convertidor CSV2plt_BLAST() y lo convierte en un array para BlastPlotter().
	array,array1,array2,array3,array4 = [],[],[],[],[]
	for i in range(len(recopilacion)//5):
		array1 += list(map(int,recopilacion[5*i+1]))
		array2 += list(map(int,recopilacion[5*i+2]))
		array3 += recopilacion[5*i+3]
		array4 += list(map(int,recopilacion[5*i+4]))
		array0 = list(range(len(array1),0,-1))
	for i in range(len(array3)):
		if array3[i] not in uniqueprots:
			array3[i] = ""
	array=[array0,array1,array2,array3,array4]
	return array

def Blast2Plotter(uniqueprots,lenprots,ID=0):
	#Convierte el archivo TSV (blast_results.tsv) de BLAST en un csv que tras pasar por CSV2plt_BLAST() se convierte en listas para BlastPlotter()
	my_input = open("Graph/Blaster_output.tsv","r")
	blastresults = my_input.read()
	blastresults = blastresults.split("\n")
	
	for i in range(len(uniqueprots)):
		my_output = open("Graph/Blast/"+uniqueprots[i],"w")
		my_output.write("#"+str(lenprots[i])+"\n")
		for hit in blastresults:
			if hit != "":
				if hit[0] != "#":
					tupla = hit.split("\t")
					if tupla[1] == uniqueprots[i]:
						my_output.write(tupla[0]+","+tupla[3]+","+tupla[5]+","+tupla[6]+"\n")
	my_input.close()
	my_output.close()
	return

def BlastPlotter(file,minpident,ID,Global=False,array=None):
	#Esta función grafica todos los resultados de BLAST: los de las proteínas y el global con todos los resultados. 
	if Global:
		array=array
	else:
		array = CSV2plt_BLAST(file,minpident)

	#Tras trasponer la matriz con numpy todos los elementos se convierten en strings, hay que convertirlas de nuevo a enteros.
	Ycor = list(map(int,array[0]))
	Witdh = list(map(int,(array[1])))
	Xcor = list(map(int,(array[2])))
	#Creamos la figura.
	fig,ax = plt.subplots(figsize=(30,18))	
	#Normalizamos el rango de colores (Spectral) para valores entre el mínimo de pident y el máximo (100%).
	normalize = mcolors.Normalize(vmin=minpident, vmax=100)
	colormap = cm.Spectral

	Spectral=[]
	#Con el siguiente bloque convertimos el valor de pident en una tupla de cuatro elementos con el color en formato RGBA.
	for pident in array[4]:
		Spectral.append(list(colormap(normalize(int(pident)))))	

	fig,ax = plt.subplots(figsize=(30,18))
	ax.barh(Ycor,Witdh,left=Xcor,color=Spectral)

	#Hacemos un poco más grandes los índices que indican la posición de las proteínas y eliminamos los ejes.
	plt.rcParams.update({'font.size':18})
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)

	#Con esto se crea la barra de color, que nos servirá de leyenda.
	scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
	cbar = plt.colorbar(scalarmappaple)
	cbar.set_label("Porcentaje de identidad",fontsize=15)

	#NOTA: El label correspondiente a la proteína hit no sale por algún fallo de matplotlib y es recurrente en StackOverFlow, 
	#pero no le encuentro solución. El nombre que corresponde a esa proteína está en el título.
	if Global:
		ax.set_title("Resultados globales", fontsize=30)
		#En el global no pongo labels porque no se ven.
		plt.savefig("Results/"+ID+"/Graphics/BLAST/Global")
	else:
		ax.set_title("Matchs de %s" %file[12:], fontsize=30)
		ax.set_yticklabels(array[3],fontsize=18,zorder=100)
		plt.savefig("Results/"+ID+"/Graphics/BLAST/"+file[12:])
	plt.close()
	return

def CSV2plt_PROSITE(file):
	#Esta función convierte los archivos csv que salen de prosite en listas para el graficador DomainPlotter().
	vector = []
	my_input = open(file,"r")
	lines = my_input.read()
	lines = lines.split("\n")
	for line in lines:
		if line != "":
			if line[0] == "#":
				lenprot = int(line[1:])
			else:
				tupla = line.split(",")
				for i in range(1,len(tupla)-1):
					vector.append([float(tupla[i]),tupla[0]])
	vector.sort()
	posiciones,dominios = [],[]
	for dupla in vector:
		posiciones.append(dupla[0])
		dominios.append(dupla[1])
	#Devuelvo las listas para DomainPlotter()
	return posiciones,dominios,lenprot

def DomainPlotter(file,ID):
	posiciones,dominios,lenprot = CSV2plt_PROSITE(file)
	#Altura de las etiquetas, para solaparlas lo menos posible.
	levels = np.array([-3,3,-6,6,-9,9,-12,12,-15,15])
	fig,ax = plt.subplots(figsize=(30,18))
	ax.plot((1,lenprot),(0,0),"k")
	color = {}
	#La primera aparición de cada dominio determinará que color tendrá en ese gráfico, pero puede variar de uno a otro.
	coloresDB=["darkviolet","darkblue","goldenrod","darkgreen","red","maroon","blue","indigo",
	           "aqua","aquamarine","chocolate","chartreuse","cyan","fuchsia","green","coral",
	           "lime","lightblue","lightgreen","navy","olive","crimson","pink","magenta"]
	
	for i in range(len(np.unique(dominios))):
		color[np.unique(dominios)[i]] = coloresDB[i % len(coloresDB)]
	for ii, (iname, ipos) in enumerate(zip(dominios,posiciones)):
		level = levels[ii % len(levels)]+randint(-3,3)
		#El level es el número que le toque y un entero aleatorio entre -3,3; para cubrir todo el rango entre un level y otro.
		if level > 0:
			vert = 'top'
		else:
			vert = 'bottom'
		#Funciones para graficar.
		ax.scatter(ipos, 0, s=100, facecolor=color[iname], edgecolor=color[iname], zorder=9999)
		ax.plot((ipos, ipos), (0, level), c=color[iname], alpha=.7)
		ax.text(ipos, level, iname,
	            horizontalalignment='center', verticalalignment=vert, fontsize=14,
	            backgroundcolor=color[iname],color="w",fontstyle="italic")
	plt.title(label="DOMINIOS DE "+file[15:],fontsize=50,style="italic",backgroundcolor="darkviolet",color="white")
	plt.setp((ax.get_yticklabels() + ax.get_yticklines() +
	list(ax.spines.values())), visible=False)
	plt.savefig("Results/"+ID+"/Graphics/PROSITE/"+file[15:])
	plt.close()
	return
