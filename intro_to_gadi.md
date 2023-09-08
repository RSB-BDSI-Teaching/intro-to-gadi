# Introduction to Gadi

## Workshop Contents 

* Logging in to Gadi
* Run FastQC analysis on DNA sequencing data 
* Align reads to reference genome 
* 


## Logging in to Gadi

There are two ways to log in to Gadi, SSH and are@NCI's web terminal. 

__1. SSH__ 

When you created an account with NCI, you will be assigned a __username__. To logging into Gadi, we can run the ssh command in our local terminal.

```sh
ssh username@gadi.nci.org.au 
```

The password would be your NCI account password. Once the password is accepted, the ssh connection will be established on one of the ten Gadi login nodes.

__2. are@NCI's web terminal__

If you don't have a local terminal and SSH installed on your laptop. You can use are@NCI's web terminal tool to log into Gadi. 

Go to this [page](https://are.nci.org.au/) and use your NCI username and password to log in. 

![are-nci-login](figures/are-nci-login.png)

After logging in, please click Gadi Terminal to start using Gadi. 

![gadi-terminal](figures/gadi-terminal.png)

## Welcome Message 

Once the connection is established, you will see a welcome message on your screen:

```
###############################################################################
#                  Welcome to the NCI National Facility!                      #
#      This service is for authorised clients only. It is a criminal          #
#      offence to:                                                            #
#                - Obtain access to data without permission                   #
#                - Damage, delete, alter or insert data without permission    #
#      Use of this system requires acceptance of the Conditions of Use        #
#      published at http://nci.org.au/users/nci-terms-and-conditions-access   #
###############################################################################
|         gadi.nci.org.au - 260,760 processor InfiniBand x86_64 cluster       |
===============================================================================

Oct 4 2022 Account Status Information on Login
   Whenever logging in from outside of Gadi, all users will now be presented
   with a brief summary about file expiry and projects near compute and/or
   storage quotas. This can be configured using the "login-info-conf" utility.
   For more information, please see https://opus.nci.org.au/x/HoABCw.

Mar 21 Sapphire Rapids Nodes Available
   An expansion to the CPU capacity of Gadi is now available with the addition
   of 720 Intel Xeon Sapphire Rapids nodes to Gadi. These nodes are available
   in the "normalsr" and "expresssr" queues.

   For more information, please see https://opus.nci.org.au/x/gIDAD
===============================================================================
Project a00 is at 97.79% of gdata inode capacity (580.71 K)
[username@gadi-login-05 ~]$ 
```

## Create a directory under /scratch/vp91

```sh
mkdir /scratch/vp91/username
cd /scratch/vp91/username
mkdir fastqc-results
nano run_fastqc.sh 
```

```sh
#!/bin/bash

#PBS -l ncpus=4
#PBS -l mem=10GB
#PBS -l jobfs=10GB
#PBS -q normal
#PBS -P vp91
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/vp91
#PBS -l wd

module load bwa/0.7.17
module load samtools/1.9
module load bcftools/1.9

fastq_dir=/scratch/vp91/ANU-Bioinfomatics-2023/data
genome=/scratch/vp91/ANU-Bioinfomatics-2023/ref-genome/ecoli_rel606.fa
out_dir=/scratch/vp91/jl6157/vc-results

for i in ${fastq_dir}/*_1.trim.fastq.gz
do
        base=$(basename $i _1.trim.fastq.gz)
        fq1=${fastq_dir}/${base}_1.trim.fastq.gz
        fq2=${fastq_dir}/${base}_2.trim.fastq.gz

        echo "Aligning sample $base"

        bam=${out_dir}/${base}.aligned.bam
        sorted_bam=${out_dir}/${base}.aligned.sorted.bam

        bwa mem -t 4 $genome $fq1 $fq2 | samtools view -S -b > $bam
        samtools sort -o $sorted_bam $bam

        echo "Calling variants for sample $base"

        raw_bcf=${out_dir}/${base}_raw.bcf
        variants=${out_dir}/${base}_variants.vcf
        final_variants=${out_dir}/${base}_final_variants.vcf

        bcftools mpileup --threads 4 -O b -o $raw_bcf -f $genome $sorted_bam
        bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
        vcfutils.pl varFilter $variants > $final_variants
done
```



## Gadi Resources 

Gadi has a lot of resources we can use to perform various tasks, and here I will introduce a few that are most relavant to us biologists. 

#### Compute Hours - granted to projects 

To run jobs on Gadi, users need to have sufficient compute hours available. The compute allocations are granted to projects instead of directly to users. Only members of a project can look up and use its compute allocation. The amount of compute hours is set by scheme manager and is belong to a project. Allocation is valid for a quarter and will be reset in the next quarter. 

To look up how much compute allocation in your project, you can run:

```sh
nci_account 
```

You should see message like this:

```
Usage Report: Project=a00 Period=2023.q3
=============================================================
    Grant:   150.00 KSU
     Used:    76.86 KSU
 Reserved:     0.00 SU
    Avail:    73.14 KSU


Storage Usage Report: Project=a00
=============================================================
Filesystem        Used     iUsed    Allocation    iAllocation
gdata1b        2.31 TiB  586.71 K     4.00 TiB       600.00 K
scratch2     520.17 GiB    4.53 K     1.00 TiB       369.63 K
=============================================================
```

__SU (Service Unit)__: is the unit that measures Gadi compute hours. Jobs run on the Gadi normal queue are charged 2 SUs to run for an hour on a single core with a memory allocation of up to 4 GiB. Jobs on Gadi are charged for the proportion of a node's total memory that is used, see more about job cost [here](https://opus.nci.org.au/display/Help/2.2+Job+Cost+Examples). In addition, different Gadi quenes have different charge rates, see the breakdown of charge rates [here](https://opus.nci.org.au/display/Help/Queue+Structure). 

#### Home Directory 



















# References 

NCI User Guide - [0. Welcome to Gadi](https://opus.nci.org.au/display/Help/0.+Welcome+to+Gadi) 



