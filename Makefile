
default: app

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

app: buildAnd.o
	$(LINK.C) -o $@ $^ ${SLEPC_EPS_LIB}
	${RM} app.o

main.o: buildAnd.c
	$(LINK.C) -c buildAnd.c ${SLEPC_EPS_LIB}


