default: 
	make dirs
	(cd src && make)
	mv src/time_combine ./bin
	mv src/ts_combine ./bin
	mv src/cgns_scale ./bin
	mv src/cgns_split ./bin
	mv src/cgns_create.so ./bin
	mv src/cgns_coarsen ./bin
	mv src/cgns_divide ./bin
	mv src/cgns_refine ./bin

dirs:
	mkdir -p bin

clean:
	@echo " Making clean ... "

	(cd src && make clean)
	rm bin/cgns_create.so
	rm bin/cgns_scale
	rm bin/cgns_split
	rm bin/time_combine
	rm bin/ts_combine
	rm bin/cgns_coarsen
	rm bin/cgns_refine
	rm bin/cgns_mirror
