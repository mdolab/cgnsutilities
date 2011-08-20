default: 
	make dirs
	(cd src && make)
	cp src/time_combine ./bin
	cp src/ts_combine ./bin
	cp src/cgns_scale ./bin
	cp src/cgns_split ./bin
	cp src/cgns_create.so ./bin
	cp src/cgns_coarsen ./bin
	cp src/cgns_divide ./bin
	mv src/time_combine ./bin
	mv src/ts_combine ./bin
	mv src/cgns_scale ./bin
	mv src/cgns_split ./bin
	mv src/cgns_create.so ./bin
	mv src/cgns_coarsen ./bin
	mv src/cgns_refine ./bin
	mv src/cgns_mirror ./bin
	mv src/time_combine ./bin
	mv src/ts_combine ./bin
	mv src/cgns_scale ./bin
	mv src/cgns_split ./bin
	mv src/cgns_create.so ./bin
	mv src/cgns_coarsen ./bin
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
