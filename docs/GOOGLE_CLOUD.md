# GOOGLE CLOUD SETUP #

### Open Web Browser ###
    Login into google account: user@gmail.com

### Go to "Google Console" ### 
    https://console.cloud.google.com/home/dashboard?project=topmed&authuser=0

### Go to "Compute Engine" ###

### Create/Select instance ###
    Machine family: General purpose
    Series: E2

### Connect using ssh ###

### Setup the machine (once) ### 

    $ ls
 
    # create software directories; set PATH
    mkdir ~/bin
    export PATH=~/bin/:$PATH
    echo "export PATH=~/bin/:\$PATH" >> ~/.bashrc
 
    # install packages : wget & samtools
    sudo apt-get update -y
    sudo apt-get install wget samtools 
 
    # download fusera
    cd bin/
    wget https://github.com/mitre/fusera/releases/download/v2.0.0/fusera 
    wget https://github.com/mitre/fusera/releases/download/v2.0.0/sracp
    chmod a+x fusera sracp 
    cd -
 
    # google config
    gsutil config 
    gcloud auth login
 
    # install gcsfuse (for bucket mounting)
    export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s`
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -
    sudo apt-get update
    sudo apt-get install gcsfuse
    sudo groupadd fuse
    sudo usermod -a -G fuse user_gmail_com

    # mount genomics-public-data
    gcsfuse --implicit-dirs genomics-public-data ~/genomics-public-data/ 
    ls ~/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta*

    # create project directory
    mkdir ~/TOPMED
    cd ~/TOPMED
    scp user@localhost:localpath/prj_17293_D27121.ngc .    

    # create mounting and output directories
    mkdir -p runs/ runs.filter/

### Filter chrM/NUMT alignments ###

    cd ~/TOPMED
    
    # create input file: 1 column containing the TOPMed Run names; should contain up to ~100 lines 
    nano runs.txt 
    head runs.txt 
    SRR...
    ...
    wc -l runs.txt 

    # mount runs
    fusera  mount --verbose --ngc  ~/TOPMED/prj_17293_D27121.ngc --accession runs.txt runs/ & 

    # find alignment files (.cram & .cram.crai)  & run "samtools view" 
    find runs/ -name "*crai" | perl -ane '/(runs\/(\w+).*).crai/; print "samtools view $1 -T ~/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta chrM -b > runs.filter/$2.bam\n";'  | tee | sh
    find runs/ -name "*crai" | perl -ane '/(runs\/(\w+).*).crai/; print "samtools view $1 -T ~/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta chr1:629084-634422 chr17:22521366-22521502 chrM -C > runs.filter/$2.cram\n";'  | tee | sh
    
    # copy filtered files to local machine
    scp -r runs.filter/ user@localhost:localpath
    
    # unmount mounting directory
    fusera unmount runs/
