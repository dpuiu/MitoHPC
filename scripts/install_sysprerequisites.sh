#!/usr/bin/env bash
set -x

##############################################################################################################

# Program that checks/installs packages needed by bwa, samtools, bedtools ...

# run using sudo; sudo privileges required !!!

##############################################################################################################

command -v apt-get
if [ "$?" == 0 ] ; then
  apt-get -y update # && apt-get upgrade
  apt-get install -y git wget default-jdk default-jre zlib1g libz-dev libncurses5-dev libbz2-dev pkg-config liblzma-dev build-essential unzip  parallel # make gcc
  apt-get install -y python   # pip
  apt-get install -y python-is-python3
fi

command -v dnf
if [ "$?" == 0 ] ; then
  dnf -y update
  dnf install -y which nano git wget java-1.8.0-openjdk bzip2 gcc gcc-c++ zlib-devel ncurses-devel bzip2-devel xz-devel unzip perl perl-Data-Dumper perl-ExtUtils-MakeMaker perl-Test-Simple python3 python3-pip make
  alternatives --install /usr/bin/python python /usr/bin/python3 60
fi

command -v yum
if [ "$?" == 0 ] ; then
  yum -y update
  #yum install -y which nano git wget java-1.8.0-openjdk bzip2 gcc gcc-c++ zlib-devel ncurses-devel bzip2-devel xz-devel  unzip perl perl-Data-Dumper  perl-ExtUtils-MakeMaker perl-Test-Simple python parallel
  yum install -y which nano git wget java-1.8.0-openjdk bzip2 gcc gcc-c++ zlib-devel ncurses-devel bzip2-devel xz-devel  unzip perl perl-Data-Dumper  perl-ExtUtils-MakeMaker perl-Test-Simple python3 python3-pip make # removed parallel python
  alternatives --install /usr/bin/python python /usr/bin/python3 60
fi
