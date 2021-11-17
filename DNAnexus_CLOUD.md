# DNAnexus CLOUD SETUP #

### Open Web Browser ###

    Login into your DNAneux account: https://platform.dnanexus.com/login

    Create a DNAnexus project (Ex: MyProject)

### Install DNAnexus Cloud Workstation App(dxpy) on Local Machine ### 

    # install/load python3
    module avail
    module load python/3.8

    # update pip
    pip install --upgrade pip --user   

    # install dxpy
    pip3 install dxpy

    # configure ssh
    PATH=$PATH:~/.local/bin/
    dx ssh_config

    # update ~/.bashrc
    nano ~/.bashrc
      PATH=$PATH:~/.local/bin/
      module load python/3.8

### Run dxpy on Local Machine => Temporary DNAnexus Client ###

    dx select MyProject

    dx run cloud_workstation --ssh -imax_session_length=8h 	# Ex: 8 hrs (default 1)
    # job-....
    # you will be logged in the DNAnexus client

### Install HP on Temporary DNAnexus Client ###

    # download precompiled version of HP; previously compiled on a DNAnexus ubuntu client
    wget ftp://ftp.ccb.jhu.edu/pub/dpuiu/HP.tgz			

    # untar
    tar -xzvf HP.tgz 

    # install sysprerequisites : java, perl, python ...
    cd HP/scripts/
    export HP_SDIR=`pwd`
    ./install_sysprerequisites.sh 

    # download hs38DH.fa from public DNAnexus project
    cd $HP_SDIR/RefSeq/
    dx download "project-BQpp3Y804Y0xbyG4GJPQ01xv:/H. Sapiens - GRCh38 with alt contigs - hs38DH/hs38DH.fa*"
    gunzip hs38DH.fa.gz

    # check install; should return Success!!!
    cd $HP_SDIR/scripts/
    ./checkInstall.sh	                                                

### Download Alignment Files on Temporary DNAnexus Client ####

     cd ~
     mkdir -p MyProject/bams
     cd MyProject/bams

     # transfer alignment files from a different DNAnexus project to DNAnexus Machine
     dx download "project-???/path-???/files-???.bam"
     cd -

### Run Pipeline on Temporary DNAnexus Client ####

     # init HP variables; create $HP_IN file
     cp $HP_SDIR/init.sh .
     . ./init.sh
    
     # create & run command file
     $HP_SDIR/run.sh > filter.all.sh                            
     bash ./filter.all.sh &

     # wait till completes ... => out/{mutect2,count,cvg}.*

### Upload Results from Temporary DNAnexus Client to DNAnexus Project ####

     dx mkdir -p out
     dx cd out
     dx upload out/mutect2.* out/count.* out/cvg.tab  

     
     

