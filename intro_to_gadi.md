# Introduction to Gadi

## Learning Objectives 

* Logging in to Gadi via SSH or are@NCI's web terminal 
* Explore scratch and gdata folder 
* Understand how to use software on Gadi 

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

__Welcome Message:__

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
[us1234@gadi-login-05 ~]$ 
```

## Explore your `home`, `scratch`, and `gdata` directories

At login, your landing point is your home directory. You will see the tilde sign `~` before the command prompt indicating you are in your home directory. 

You can also use `pwd` to check what the full path of your home directory is. 

```sh
[us1234@gadi-login-05 ~]$ pwd
/home/001/us1234
```

There are 3 places where you can store your data on Gadi:

* Your home directory, which has 10 GiB storage.
* The `/scratch/project` directory, which has 72 GiB storage by default. Files not accessed for more than 100 days are automatically moved from project directories. 
* The `/g/data/project` directory, the amount of storage is set by the scheme manager. You can log in to [MyNCI](http://my.nci.org.au/) to check your project's allocation. 

__Let's go to the `/scratch` directory.__

```sh
cd /scratch 
ls
```

You should see results look like this:

```
public  project1  project2
```

The different directories belong to different projects. 

__Let's go to the `/g/data` directory.__

```sh
cd /g/data
ls
```

The result should look like this:

```
a01  a02  a03  a04  a05  a06  a07  a08  a09  a10  
a11  a12  a13  a14  a15  a16  a17  a18  a19  a20
```

The directories are belong to different projects, and if you have enrolled in a project, you should find your project directory here. The `training project` we enrolled in does not have a gdata folder so you won't find a directory here. 

__Check how much storage you have access to using `lquota`:__

```
[us1234@gadi-login-06 ~]$ lquota
--------------------------------------------------------------------------
           fs       Usage      Quota      Limit   iUsage   iQuota   iLimit
--------------------------------------------------------------------------
   a00 scratch 520.17 GiB   1.00 TiB   1.05 TiB     4531   369634   388115
   a01 scratch   2.11 TiB   5.00 TiB   5.25 TiB   948753 10485760 11010048
   a00   gdata   2.26 TiB   4.00 TiB   4.20 TiB   588867   600000   630000
--------------------------------------------------------------------------
```

Apart from storage limitation, there is also a quota called `iQuota` applied to the storage allocation on `/scratch` and `/g/data`. It sets the maximum number of files allowed in the project which `iUsage` shows the existing number of files and `iQuota` shows the file-number limit. 

__Checking personal usage on projects using `nci-files-report`:__

```
[us1234@gadi-login-03 ~]$ nci-files-report
------------------BREAKDOWN BY PROJECT, GROUP, AND USER-------------------
FILESYSTEM  SCAN DATE   PROJECT  GROUP  USER    SPACE USED  TOTAL SIZE  COUNT
scratch     2023-09-14  a00      a00    us1234      882.5M      882.5M  15
scratch     2023-09-14  a01      a01    us1234        8.0K        8.0K  2
--------------------------------------------------------------------------

------------------SUMMARIES-------------------
FILESYSTEM  USER    SPACE USED  TOTAL SIZE  COUNT
scratch     us1234      882.5M      882.5M  17
----------------------------------------------
```

It returns a summary of your usage on each project and on each folder. If you don't have any files created, it will return nothing. 

You can also use options to limit the results shown, try `nci-files-report -h` to learn more about the command. 






## Job Folder $PBS_JOBFS 




## Run the variant-calling workflow 

```sh
mkdir /scratch/vp91/username
cd /scratch/vp91/username
mkdir fastqc-results
nano run_fastqc.sh 
```

```sh
#!/bin/bash

#PBS -l ncpus=8
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

fastq_dir=/scratch/vp91/ANU-Bioinformatics-2023/data
genome=/scratch/vp91/ANU-Bioinformatics-2023/ref-genome/ecoli_rel606.fa
out_dir=~/variant-calling/results

mkdir -p $out_dir 

for i in ${fastq_dir}/*_1.trim.fastq.gz
do
        base=$(basename $i _1.trim.fastq.gz)
        fq1=${fastq_dir}/${base}_1.trim.fastq.gz
        fq2=${fastq_dir}/${base}_2.trim.fastq.gz

        echo "Aligning sample $base"

        bam=${out_dir}/${base}.aligned.bam
        sorted_bam=${out_dir}/${base}.aligned.sorted.bam

        bwa mem -t 8 $genome $fq1 $fq2 | samtools view -S -b > $bam
        samtools sort -o $sorted_bam $bam

        echo "Calling variants for sample $base"

        raw_bcf=${out_dir}/${base}_raw.bcf
        variants=${out_dir}/${base}_variants.vcf
        final_variants=${out_dir}/${base}_final_variants.vcf

        bcftools mpileup --threads 8 -O b -o $raw_bcf -f $genome $sorted_bam
        bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
        vcfutils.pl varFilter $variants > $final_variants
done
```

Save and exit the script. 

__Submitting the job using `qsub`__:

```sh
qsub run_vc.sh 
```

## Checking your job status 




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



