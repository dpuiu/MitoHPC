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

### View DNAnexus Files ###
  
    # login using DNAnexus username and password
    dx login
  
    # list & select project
    dx select

    # view UKB results
    dx ls
    dx ls dpuiu/
    dx ls dpuiu/bams/
    dx ls dpuiu/out/
    dx ls dpuiu/out/TAR/

    dx cat dpuiu/out/mutect2.03.vcf
    
### Create Temporary DNAnexus/AWS Client ###

    # create client : Ex: 4 hrs (default 1hr)
    dx run cloud_workstation --ssh -y -imax_session_length=4h # -instance-type mem1_ssd1_v2_x8 --name "default:4hr"
    # job-...
    # you will be logged in the DNAnexus client

    # reconnect terminal if the connection drops
    dx ssh job-...

    # terminate early
    dx terminate job-...
  
### Install Software on DNANexus Client ###

    # install dxfuse for easy file access
    wget https://github.com/dnanexus/dxfuse/releases/download/v1.0.0/dxfuse-linux
    chmod u+x dxfuse-linux

    sudo apt-get update 
    sudo apt-get install parallel samtools

 ### Mount Project(s) Data ###

    mkdir ~/Ref ~/UKB
    ./dxfuse-linux ~/Ref/ project-BQpp3Y804Y0xbyG4GJPQ01xv        # reference assemblies
    ./dxfuse-linux ~/UKB/ project-G7KB5zQJz55qG0158ZXb2J5p        # UKB project which includes the 200K samples

    # create symlink to cram, bam file locations
    ln -s UKB/UKB\ -mtDNA_2022/Bulk/Whole\ genome\ sequences/Whole\ genome\ CRAM\ files/ crams
    ln -s UKB/UKB\ -mtDNA_2022/dpuiu/bams  

### Install HP pipeline  ###

#### From Github ####

    # git clone
    git clone https://github.com/dpuiu/HP.git

    # install
    cd HP/scripts/
    export HP_SDIR=$PWD
    . ./init.sh
    sudo ./install_sysprerequisites.sh
    ./install_prerequisites.sh  
    ./check_install.sh
    # should returm Success!

#### Precompiled ####
    
    # login into DNANeux account
    dx select 
    dx select                               # run twice to generate ~/.dnanexus_config/unsetenv
    source ~/.dnanexus_config/unsetenv
    dx clearenv
    dx logout
    dx login

    # download precompiled version of HP; previously compiled on a DNAnexus ubuntu client
    dx download dpuiu/HP.tgz			
    tar -xzvf HP.tgz 

    # install sysprerequisites : java, perl, python ...
    cd HP/scripts/
    export HP_SDIR=`pwd`
    . ./init.sh
    sudo ./install_sysprerequisites.sh 
    ./checkInstall.sh
    # should returm Success!

### Run Pipeline ####

     # init HP variables; create $HP_IN file
     cp $HP_SDIR/init.sh .

     # edit init.sh
     nano init.sh
       HP_CN=
       HP_L=111000
       ...

     # generate in.txt
     cat crams/newlistofallfiles.txt | sed 's|/home/dnanexus/wgs|crams|' | $HP_SDIR/ls2in.pl -out out > in.txt

     # re-init
     . ./init.sh
    
     # create command file
     $HP_SDIR/run.sh > run.all.sh                            

     # run in the background
     bash ./run.all.sh &

     # or run in paralled
     grep filter.sh ./run.all.sh | parallel
     
     # or run only a subset of the jobs 
     grep filter.sh ./run.all.sh | grep -P "/10/' | parallel

     # get summaries
     getSummary.sh

     # wait till completes ... => out/{mutect2,count,cvg}.*

### Save Results ###

     # upload Results from Temporary DNAnexus Client to DNAnexus Project
     dx cd ...
     dx upload ...

