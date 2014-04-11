default: 
	(cd src && make)

clean:
	@echo " Making clean ... "
	rm -fr src/*.o
	rm -fr src/*.c
	rm -fr src/*.pyf