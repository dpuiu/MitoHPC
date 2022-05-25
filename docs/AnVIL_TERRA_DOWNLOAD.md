# DOWNLOADING GTEx DATA #

## WEB BROWSER ##

* Go to https://anvil.terra.bio/ and login 
* Go to Workspaces , and select AnVIL_GTEx_V8_hg38,GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files ... then Select a Billing Project
* Select multiple files to download, copy the gsutil command 

## LOCAL LINUX MACHINE ##

    # Download and install gcloud CLI
    wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-371.0.0-linux-x86_64.tar.gz
    tar -xzvf google-cloud-sdk-371.0.0-linux-x86_64.tar.gz 
    cd google-cloud-sdk/
    ./install.sh 
    cd ..

    #update gcloud CLI
    gcloud components update

    # Login and Configure
    PROJECT=topmed                          # gtexproject
    gcloud auth login --no-launch-browser   # a broser window will pop up; select the account you are planning to use
      
    gcloud auth list

    gcloud config set project $PROJECT
    gcloud config list
  
    # Run the gsutil commands , adding the -u option for billing project
    gsutil -u $PROJECT ls
    gsutil -u $PROJECT ls gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/

    # WGS and WES files   
    gsutil -u $PROJECT ls gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/ 
      gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files/
      gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WES_BAM_files/
    gsutil -u $PROJECT -m cp  ""gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files/GTEX-1117F-0003-SM-6WBT7.cram" "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files/GTEX-1117F-0003-SM-6WBT7.crai"   .
   
    # RNAseq files
    gsutil -u $PROJECT ls gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/
      
    gsutil -u $PROJECT -m cp "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/GTEX-1128S-0005-SM-5P9HI.Aligned.sortedByCoord.out.patched.md.bam" "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/GTEX-1128S-0005-SM-5P9HI.Aligned.sortedByCoord.out.patched.md.bam.bai" .
