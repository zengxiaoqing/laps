#  Make file for Nexrad RISC RPG Project
#
#  Version 1.00 - Preliminary program development
#  Version 1.01 - First Registered version in SCCS
#  Version 1.02 - Added error_log to file list and renamed the output to rpg
#  Version 1.03 - Added tapesys, tape thread simulator
#
rpg: rpg_tape monitor widebd_server widebd_client serial_HSI.o minsys.o bx25.o sysmon.o  wideband.o mess.o error_log.o riscrpg.h mess.h x25.h makefile
	cc -g -o rpg minsys.o serial_HSI.o bx25.o sysmon.o  wideband.o mess.o error_log.o -I/usr/local/include -L/nsslsun/nssl/bens/lib -lnssl

widebd_server: bd_server.c wideband.o  makefile 
	cc -g -o widebd_server bd_server.c wideband.o  -lnssl

widebd_client: bd_client.c wideband.o  makefile
	cc -g -o widebd_client bd_client.c wideband.o  -lnssl

rpg_tape: rpg_tape.c get_tape_mess.o sysmon.o  wideband.o mess.o riscrpg.h mess.h makefile
	cc -g -o rpg_tape rpg_tape.c get_tape_mess.o sysmon.o  wideband.o mess.o -I/usr/local/include -L/nsslsun/nssl/bens/lib -lnssl

monitor:  monitor.c sysmon.o riscrpg.h makefile
	/usr/5bin/cc -g -o monitor monitor.c -I/usr/local/include -lnssl -lcurses -ltermcap

minsys.o: minsys.c riscrpg.h x25.h makefile
	cc -g -c minsys.c -I/usr/local/include

mess.o: mess.c mess.h makefile
	cc -g -c mess.c -I/usr/local/include

serial_HSI.o: serial_HSI.c serial_HSI.h makefile
	cc -g -c serial_HSI.c -I/usr/local/include

wideband.o: wideband.c riscrpg.h makefile
	cc -g -c wideband.c -I/usr/local/include

bx25.o: bx25.c x25.h makefile
	cc -g -c bx25.c -I/usr/local/include

sysmon.o: sysmon.c riscrpg.h makefile
	cc -g -c sysmon.c -I/usr/local/include

error_log.o: error_log.c riscrpg.h makefile
	cc -g -c error_log.c -I/usr/local/include

get_tape_mess.o: get_tape_mess.c mess.h riscrpg.h makefile
	cc -g -c get_tape_mess.c -I/usr/local/include

riscrpg.h: mess.h x25.h makefile

# shm_attach.o: shm_attach.c makefile
#	cc -g -c shm_attach.c -I/usr/local/include -lnssl

testshm: testshm.o makefile 
	cc -g -o testshm testshm.o -I/usr/local/include
testcircw: testcircw.o  makefile 
	cc -g -o testcircw testcircw.o  -I/usr/local/include
testcircr: testcircr.o  makefile 
	cc -g -o testcircr testcircr.o  -I/usr/local/include


sample_wide_reader: sample_wide_reader.o  wideband.o makefile 
	cc -g -o sample_wide_reader sample_wide_reader.o  wideband.o -I/usr/local/include

sample_wide_writer: sample_wide_writer.o  wideband.o mess.h makefile 
	cc -g -o sample_wide_writer sample_wide_writer.o  wideband.o -I/usr/local/include

wideband_echo: wideband_echo.o  wideband.o mess.h makefile 
	cc -g -o wideband_echo wideband_echo.o  wideband.o -I/usr/local/include

wideband_scan: wideband_scan.o  wideband.o mess.h makefile 
	cc -g -o wideband_scan wideband_scan.o  wideband.o -I/usr/local/include

wideband_rda: wideband_rda.o  wideband.o mess.h makefile 
	cc -g -o wideband_rda wideband_rda.o  wideband.o -I/usr/local/include




