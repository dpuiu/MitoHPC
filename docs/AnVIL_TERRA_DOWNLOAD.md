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

    # Login
    PROJECT=topmed                          # gtexproject
    gcloud auth login --no-launch-browser   # a broser window will pop up; select the account you are planning to use
    gcloud auth list

    # Configure
    gcloud config set project $PROJECT
    gcloud config list

### GTEX ###

    # Run the gsutil commands , adding the -u option for billing project
    gsutil -u $PROJECT ls gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/

    # WGS and WES files: list and transfer
    gsutil -u $PROJECT ls gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/
      gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files/
      gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WES_BAM_files/
    gsutil -u $PROJECT -m cp  ""gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files/GTEX-1117F-0003-SM-6WBT7.cram" \
      "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files/GTEX-1117F-0003-SM-6WBT7.crai"   .
   
    # RNAseq files : list and transfer
    gsutil -u $PROJECT ls gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/
    gsutil -u $PROJECT -m cp "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/GTEX-1128S-0005-SM-5P9HI.Aligned.sortedByCoord.out.patched.md.bam" \
      "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/GTEX-1128S-0005-SM-5P9HI.Aligned.sortedByCoord.out.patched.md.bam.bai" .

### CCDG_ARIC ###

    # List files
    gsutil -u $PROJECT ls gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7

    gsutil -u $PROJECT ls gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7 | head
      gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/aric.txt
      gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/A00010/
      gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/A00016/
      gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/A00068/
      ...

    gsutil -u $PROJECT cat gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/aric.txt | grep cram$ | wc -l
      8976

    gsutil -u $PROJECT cat gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/aric.txt | head
      gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/CRAMs/A00003/A00003.hgv.cram
      gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/CRAMs/A00004/A00004.hgv.cram
      gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/CRAMs/A00014/A00014.hgv.cram
      ...

     # View CRAM file
     gsutil -u $PROJECT cat gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/CRAMs/A00003/A00003.hgv.cram | samtools view

     # Mount bucket (not working yet)
     mkdir fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7
     gcsfuse --billing-project $PROJECT  gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7 fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7
     ls fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/ | wc -l
        0
     # empty ; check billing enabled https://cloud.google.com/billing/docs/how-to/verify-billing-enabled ???


     # Download CRAM and CRAI files
     gsutil -u $PROJECT cp  gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/CRAMs/A00003/A00003.hgv.cram .
     gsutil -u $PROJECT cp  gs://fc-secure-f5d884c0-a24c-46e6-8c29-cad7f5b158c7/CRAMs/A00003/A00003.hgv.cram.crai .
     ...

     # Filter chrM reads
     samtools view A00003.hgv.cram chrM -b > A00003.hgv.chrM.bam

