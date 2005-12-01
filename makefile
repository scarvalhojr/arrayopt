# compilers and tools

# Solaris (CeBiTec network)
GCC = /vol/gcc-3.4.3/bin/gcc
JDK_PATH = /vol/java-1.5.0

# Windows (Cygwin)
#GCC = gcc
#JDK_PATH = /usr

# Java tools
JAVAC = $(JDK_PATH)/bin/javac
JAVADOC = $(JDK_PATH)/bin/javadoc
JAVAH = $(JDK_PATH)/bin/javah

# directories
SRC_DIR = ./src
BIN_DIR = ./bin
DOC_DIR = ./doc
INC_DIR = ./include
QAP_DIR = $(SRC_DIR)/qap

# windows libraries
W32_LIB = ./lib/win32
W32_INC = -I$(INC_DIR) -I$(JDK_PATH)/include -I$(JDK_PATH)/include/win32

# solaris libraries
SOL_LIB = ./lib/solaris
SOL_INC = -I$(INC_DIR) -I$(JDK_PATH)/include -I$(JDK_PATH)/include/solaris

# source files and packages
JAVA_SRC = $(SRC_DIR)/arrayopt/layout/*.java $(SRC_DIR)/arrayopt/textui/*.java
JAVA_BIN = $(SRC_DIR)/arrayopt/layout/*.class $(SRC_DIR)/arrayopt/textui/*.class
JAVA_PCK = arrayopt.layout arrayopt.textui
JNI_CLS = arrayopt.layout.QAPGraspDense

# phony targets (always execute)
.PHONY: classes apidoc jni cleanbin cleandoc cleanlib

# ==============================================================================

# build Java classes
classes: $(BIN_DIR)
	$(JAVAC) -d $(BIN_DIR) $(JAVA_SRC)

# build API documentation
apidoc: $(DOC_DIR)
	$(JAVADOC) -d $(DOC_DIR) -sourcepath $(SRC_DIR) @apidoc_options $(JAVA_PCK)

# generate JNI header files
jni: arrayopt
	$(JAVAH) -jni -d $(INC_DIR) -classpath $(BIN_DIR) $(JNI_CLS)

# build everything (libraries too?)
all: classes apidoc

# ==============================================================================

# under development...

# build Windows libraries of external implementations
# libwin: $(W32_LIB) $(W32_LIB)/qap_graspd.dll

# external library (Windows): GRASP for dense QAP
# $(W32_LIB)/qap_graspd.dll: $(QAP_DIR)/graspd/qap_graspd.c $(QAP_DIR)/graspd/qap_graspd.f
#	$(GCC) -Df2cFortran -Xlinker -G -O3 -o $@ $(SOL_INC) $^ -lfrtbegin -lg2c

# ==============================================================================

# build Solaris libraries of external implementations
libsol: $(SOL_LIB) $(SOL_LIB)/libqap_graspd.so

# external library (Solaris): GRASP for dense QAP
$(SOL_LIB)/libqap_graspd.so: $(QAP_DIR)/graspd/qap_graspd.c $(QAP_DIR)/graspd/qap_graspd.f
	$(GCC) -G -O3 -o $@ $(SOL_INC) $^

# To do:
# Use '-Xlinker -G' or '-shared' (or '-symbolic'?) instead of just '-G'
# for building shared libraries

# ==============================================================================

# make sure directories are created
$(BIN_DIR) $(DOC_DIR) $(SOL_LIB) $(W32_LIB):
	mkdir -p $@

# remove files
clean: cleanbin cleandoc cleanlib

cleanbin:
	rm -r -f $(BIN_DIR)/*

cleandoc:
	rm -r -f $(DOC_DIR)/*

cleanlib:
	rm -r -f $(W32_LIB)/*
	rm -r -f $(SOL_LIB)/*
