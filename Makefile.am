ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = eca-shared

bin_PROGRAMS = tpm beta_sim mc3demo inbred_p newhybs pritch



# CFLAGS for each target
tpm_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) -Ieca-shared/include -Wno-deprecated
beta_sim_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) -Ieca-shared/include -Wno-deprecated
mc3demo_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) -Ieca-shared/include -Wno-deprecated
inbred_p_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) -Ieca-shared/include -Wno-deprecated
newhybs_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) -Ieca-shared/include -Wno-deprecated
pritch_CFLAGS =  $(GL_CFLAGS) $(GLU_CFLAGS) $(GLUT_CFLAGS) \
	-Isrc -Inewhybrids/src -Ieca-shared/include -Wno-deprecated



# LDFLAGS for each target
tpm_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm
beta_sim_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm
mc3demo_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm
inbred_p_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm
newhybs_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm
pritch_LDFLAGS = $(GL_LIBS) $(GLU_LIBS) $(GLUT_LIBS) -lm


# which convenience libraries to link to each target
tpm_LDADD = eca-shared/libecashared.a  eca-shared/libgfmcmc.a
beta_sim_LDADD = eca-shared/libecashared.a  eca-shared/libgfmcmc.a
mc3demo_LDADD = eca-shared/libecashared.a  eca-shared/libgfmcmc.a
inbred_p_LDADD = eca-shared/libecashared.a  eca-shared/libgfmcmc.a eca-shared/libecaopt3.a
newhybs_LDADD = eca-shared/libecashared.a  eca-shared/libgfmcmc.a eca-shared/libecaopt3.a
pritch_LDADD = eca-shared/libecashared.a  eca-shared/libgfmcmc.a eca-shared/libecaopt3.a



# sources for each target
tpm_SOURCES = src/tpm.c
beta_sim_SOURCES = src/beta_sim.c 
mc3demo_SOURCES = src/mc3demo.c
inbred_p_SOURCES = src/inbred_p.c
newhybs_SOURCES = newhybrids/src/GLUT_for_NewHybrids.c  \
	newhybrids/src/DataInput.c  \
	newhybrids/src/DataOutput.c  \
	newhybrids/src/NewHybrids.c  \
	newhybrids/src/RunWOG_func.c
pritch_SOURCES = src/GLUT_for_PritchEtc.c \
	src/MAZ4NewHybs.h \
	src/MargAllocZ.c \
	src/MargAllocDataProb.c \
	src/PritchFuncts.c  \
	newhybrids/src/NewHybrids.c \
	newhybrids/src/DataInput.c \
	newhybrids/src/DataOutput.c