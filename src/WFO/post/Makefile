xfer_laps : xfer_laps.o fill_laps.o /usr/local/netcdf/lib/libnetcdf.a
	cc -o xfer_laps -g xfer_laps.o fill_laps.o /usr/local/netcdf/lib/libnetcdf.a

xfer_laps.o : xfer_laps.c
	cc -c -g -Aa -D_INCLUDE_POSIX_SOURCE xfer_laps.c

fill_laps.o : fill_laps.c
	cc -c -g -Aa -D_INCLUDE_POSIX_SOURCE fill_laps.c
