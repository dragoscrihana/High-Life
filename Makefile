all: homework

homework: homework.c
	mpicc -o homework homework.c