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
  
    # login using DNAnexus username and password
    dx login
  
    # list & select project
    dx select
    
    # create client : Ex: 4 hrs (default 1hr)
    dx run cloud_workstation --ssh -y -imax_session_length=4h # -instance-type mem1_ssd1_v2_x8 --name "default:4hr"
    # job-...
    # you will be logged in the DNAnexus client

    # reconnect terminal if the connection drops
    # dx ssh job-...

    # terminate early
    # dx terminate job-...
  
### Install Software on DNANexus Client ###

    # install dxfuse for easy file access
    wget https://github.com/dnanexus/dxfuse/releases/download/v1.0.0/dxfuse-linux
    chmod u+x dxfuse-linux

    sudo apt-get update 
    sudo apt-get install parallel samtools bwa bcftools samblaster bedtools fastp tabix

 ### Mount Project(s) Data ###

    mkdir ~/Ref ~/UKB
    ./dxfuse-linux ~/Ref/ project-BQpp3Y804Y0xbyG4GJPQ01xv        # reference assemblies
    ./dxfuse-linux ~/UKB/ project-G7KB5zQJz55qG0158ZXb2J5p        # UKB project which includes the 200K samples

### Install HP pipeline  ###

#### From Github ####

    # git clone
    git clone https://github.com/dpuiu/HP.git

    # install
    cd HP/scripts/
    export HP_SDIR=$PWD
    sudo ./install_sysprerequisites.sh
    ./install_prerequisites.sh  

#### Precompiled ####
    
    # download precompiled version of HP; previously compiled on a DNAnexus ubuntu client
    wget ftp://ftp.ccb.jhu.edu/pub/dpuiu/HP.tgz			
    tar -xzvf HP.tgz 

    # install sysprerequisites : java, perl, python ...
    cd HP/scripts/
    export HP_SDIR=`pwd`
    sudo ./install_sysprerequisites.sh 
    cd $HP_SDIR/RefSeq/

    # download or link the hs38DH.fa reference
    dx download "project-BQpp3Y804Y0xbyG4GJPQ01xv:/H. Sapiens - GRCh38 with alt contigs - hs38DH/hs38DH.fa*"       
    #or
    ln -s Ref/'Reference Genome Files: AWS US (East)'/'H. Sapiens - GRCh38 with alt contigs - hs38DH'/hs38DH.fa.gz

    zcat hs38DH.fa.gz > hs38DH.fa

 #### Check Pipeline Install ####

    cd $HP_SDIR/scripts/

    ./checkInstall.sh	          
    # should returm Success!                                      

    cd
    tar -czvf HP.tgz HP/
    dx upload HP.tgz

### Download/Link Alignment Files  ####

     cd ~
     mkdir crams
     cd crams

     # transfer alignment files from a different DNAnexus project to DNAnexus Machine
     dx download "project-???/path-???/files-???.cram*"

     # or just link then (if dxfuse-linux previously setup)
     ln -s UKB/UKB\ -mtDNA_2022/Bulk/Whole\ genome\ sequences/Whole\ genome\ CRAM\ files/ crams
     wc -l crams/newlistofallfiles.txt 
     200029 


### Run Pipeline ####

     # init HP variables; create $HP_IN file
     cp $HP_SDIR/init.sh .

     # edit init.sh
     nano init.sh
       HP_CN=
       HP_L=111000
       ...

     # generate in.txt
     head ... crams/newlistofallfiles.txt | sed 's|/home/dnanexus/wgs|crams|' | $HP_SDIR/ls2in.pl -out out > in.tx

     # re-init
     . ./init.sh
    
     # create command file
     $HP_SDIR/run.sh > run.all.sh                            

     # run in the background
     bash ./run.all.sh &

     #or run in paralled
     grep filter.sh ./run.all.sh | parallel
     getSummary.sh

     # wait till completes ... => out/{mutect2,count,cvg}.*

### Save Results ###

     # upload Results from Temporary DNAnexus Client to DNAnexus Project
     dx mkdir -p out
     dx cd out
     dx upload out/mutect2.* out/count.* out/cvg.tab  

     # or scp copy to remote machine
     scp out/mutect2.* out/count.* out/cvg.tab user@remotehost:...
