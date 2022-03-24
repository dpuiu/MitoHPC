#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program that checks/installs packages needed by bwa, samtools, bedtools ...

# sudo privileges required !!!

##############################################################################################################

sudo -v
if [ "$?" == 1 ] ; then
  echo "ERROR: sudo permission required"
  exit 1
fi

which apt-get
if [ "$?" == 0 ] ; then
  sudo apt-get -y update && apt-get upgrade
  sudo apt-get install -y git wget default-jdk default-jre zlib1g libz-dev libncurses5-dev libbz2-dev pkg-config liblzma-dev python  build-essential unzip parallel # make gcc
fi

which yum
if [ "$?" == 0 ] ; then
  sudo yum update
  sudo yum install -y which nano git wget java-1.8.0-openjdk bzip2 gcc gcc-c++ zlib-devel ncurses-devel bzip2-devel xz-devel  unzip perl perl-Data-Dumper  perl-ExtUtils-MakeMaker perl-Test-Simple python parallel
fi
