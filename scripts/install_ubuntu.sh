#!/bin/bash -eux

sudo apt-get -y update && apt-get upgrade
sudo apt-get install -y git wget default-jdk make zlib1g libz-dev libncurses5-dev libbz2-dev pkg-config liblzma-dev
