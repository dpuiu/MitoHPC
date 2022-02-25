# DOWNLOADING GTEx DATA #

## WEB BROWSER ##

* Go to https://anvil.terra.bio/ and login 
* Go to Workspaces , and select AnVIL_GTEx_V8_hg38,GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files ... then Select a Billing Project
* Select multiple files to download, copy the gsutil command 

## LOCAL LINUX MACHINE ##

* Download and install gcloud CLI

    wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-371.0.0-linux-x86_64.tar.gz
    tar -xzvf google-cloud-sdk-371.0.0-linux-x86_64.tar.gz 
    cd google-cloud-sdk/
    ./install.sh 
    cd ..

* Login and Configure

    gcloud auth login
    gcloud auth list
    gcloud config set project ...
    gcloud config list

* Run the gsutil command , adding the -u option for billing project

    gsutil -u $PROJECT -m cp  "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/file1.bam" "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/file2.bam.bai"   .
