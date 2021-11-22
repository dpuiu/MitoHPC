# DNAnexus CLOUD SETUP #

### Open Web Browser ###

    Login into your DNAneux account: https://platform.dnanexus.com/login

    Create a DNAnexus project (Ex: MyProject)

### Install DNAnexus Cloud Workstation App(dxpy) on Local Machine ### 

    # open terminal app

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
  
    # select project
    dx select MyProject

    # create client : Ex: 8 hrs (default 1hr)
    dx run cloud_workstation --ssh -imax_session_length=8h 	
    # job-...
    # you will be logged in the DNAnexus client

    # install dxfuse for easy file access
    wget https://github.com/dnanexus/dxfuse/releases/download/v1.0.0/dxfuse-linux

    # create mount points & mount remote project directories
    mkdir ~/Ref/ ~/UKbiobank
    dxfuse-linux  ~/Ref/ project-BQpp3Y804Y0xbyG4GJPQ01xv  # reference assemblies
    dxfuse-linux ~/UKbiobank ...                           # alignment files

### Install HP pipeline ###

    # download precompiled version of HP; previously compiled on a DNAnexus ubuntu client
    wget ftp://ftp.ccb.jhu.edu/pub/dpuiu/HP.tgz			

    # untar
    tar -xzvf HP.tgz 

    # install sysprerequisites : java, perl, python ...
    cd HP/scripts/
    export HP_SDIR=`pwd`
    ./install_sysprerequisites.sh 
    cd $HP_SDIR/RefSeq/

    # download or link the hs38DH.fa reference
    dx download "project-BQpp3Y804Y0xbyG4GJPQ01xv:/H. Sapiens - GRCh38 with alt contigs - hs38DH/hs38DH.fa*"       
    #or
    ln -s Ref/'Reference Genome Files: AWS US (East)'/'H. Sapiens - GRCh38 with alt contigs - hs38DH'/hs38DH.fa.gz

    zcat hs38DH.fa.gz > hs38DH.fa

    # check install; should return Success!!!
    cd $HP_SDIR/scripts/
    ./checkInstall.sh	                                                


### Download/Link Alignment Files  ####

     cd ~
     mkdir -p MyProject/bams
     cd MyProject/bams

     # transfer alignment files from a different DNAnexus project to DNAnexus Machine
     dx download "project-???/path-???/files-???.bam"

     # or just link then (if dxfuse-linux previously setup)
     ln -s ~/UKbiobank/.../*bam .

     cd -

### Run Pipeline ####

     # init HP variables; create $HP_IN file
     cp $HP_SDIR/init.sh .
     . ./init.sh
    
     # create & run command file
     $HP_SDIR/run.sh > run.all.sh                            
     bash ./run.all.sh &

     # wait till completes ... => out/{mutect2,count,cvg}.*

### Save Results ###

     # upload Results from Temporary DNAnexus Client to DNAnexus Project
     dx mkdir -p out
     dx cd out
     dx upload out/mutect2.* out/count.* out/cvg.tab  

     # or scp copy to remote machine
     scp out/mutect2.* out/count.* out/cvg.tab user@remotehost:...
