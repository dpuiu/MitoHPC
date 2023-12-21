FROM ubuntu:23.04

ENV DEBIAN_FRONTEND=noninteractive

ENV HP_HDIR=/MitoHPC/
ENV HP_SDIR=/MitoHPC/scripts/
ENV HP_BDIR=/MitoHPC/bin/
ENV HP_ADIR="bams" 
ENV HP_ODIR="out"
ENV HP_IN="in.txt" 

ENV PATH="$HP_SDIR:$HP_BDIR:$PATH"

###########################################

RUN apt-get -y update
RUN apt-get install -y wget tar nano curl git

###########################################
RUN \
  git clone https://github.com/dpuiu/MitoHPC && \
  chmod a+x $HP_SDIR/*.* && \
  . $HP_SDIR/init.sh && \
  $HP_SDIR/install_sysprerequisites.sh && \
                                                                     $HP_SDIR/install_prerequisites.sh && \
  export HP_MT=rCRS && export HP_MTC=rCRSC && export HP_MTR=rCRSR && $HP_SDIR/install_prerequisites.sh && \
  export HP_MT=RSRS && export HP_MTC=RSRSC && export HP_MTR=RSRSR && $HP_SDIR/install_prerequisites.sh && \
  HP_MT=RSRS $HP_SDIR/install_prerequisites.sh && \
  $HP_SDIR/checkInstall.sh && \
  rm -fr /MitoHPC/prerequisites/ /MitoHPC/examples*
