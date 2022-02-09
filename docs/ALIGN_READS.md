# ALIGN READS USING BWA # 

## INPUT ##
 
    Illumina paired-end reads .fq.gz
    Reference genome

### SETUP ENVIRONMENT ###

    # init environment (check README.md) 
    test -s init.sh  
    . ./init.sh

    # index reference if necessary    
    bwa_index.sh

    mkdir -p $HP_ADIR

    # generate alignment script
    find $HP_FDIR -name "*_1.fq*" | perl -lane '/(.+\/(.+))_1.fq/; print "$2\t$1\t$ENV{HP_ADIR}/$2";' | sort | cut -f2,3 | sed 's|^|bwa_mem.sh |' > align.all.sh
    bash align.all.sh
