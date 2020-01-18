#!/usr/bin/env python

from subprocess import call, Popen, PIPE

def Aligmenter(fasta_file,output="aligment"):
	#Esta función alinea los multifasta de cada query.
	proceso = Popen(["muscle","-in",fasta_file ], stdout=PIPE, stderr=PIPE)
	error_encontrado = proceso.stderr.read().decode("utf-8")
	listado = proceso.stdout.read().decode("utf-8")

	proceso.stderr.close()
	proceso.stdout.close()

	my_output = open(output,"w")
	my_output.write(str(listado))
	my_output.close()
	return 


def TreeMaker(output, align_file="aligment"):
	#Esta función hace los árboles, que tendrán el mismo nombre que la query pero con formato .tre
	output += ".tre"
	proceso = Popen(["muscle","-maketree","-in",align_file,"-out",output,"-cluster", "neighborjoining"], stderr=PIPE)
	error_encontrado = proceso.stderr.read().decode("utf-8")
	#Borro el archivo de alineamiento intermedio.
	call (["rm",align_file])
	return
