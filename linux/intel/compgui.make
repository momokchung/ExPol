#
#
#  ################################################################
#  ##                                                            ##
#  ##  compgui.make  --  compile Tinker routines needed for FFE  ##
#  ##              (Intel Fortran for Linux Version)             ##
#  ##                                                            ##
#  ################################################################
#
#
#  This assumes $JAVA_HOME is set to the JAVA install directory,
#  which is /usr/lib/jvm/default-java on many Linux systems
#
#
icc -c -O3 -no-prec-div -static -w server.c -I $JAVA_HOME/include -I $JAVA_HOME/include/linux
