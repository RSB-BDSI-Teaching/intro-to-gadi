# Introduction to Gadi

## Learning Objectives 

## Logging in to Gadi

When you created an account with NCI, you will be assigned a __username__. To logging into Gadi, we can run the ssh command in our local terminal.

```sh
ssh username@gadi.nci.org.au 
```

The password would be your NCI account password. Once the password is accepted, the ssh connection will be established on one of the ten Gadi login nodes, __gadi-login-xx__, where xx=[01, 10]. Note that Gadi login nodes are separate from the Gadi compute notes on which jobs actually run. 

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

__Message of the Day:__

Under the welcome text is the Message of the Day. Please consider it a notice board and read it on every login. News and status updates relevant to Gadi users will be posted on it. 

__Usage Limit on Login Nodes__

Users can only run a limited types of commands on the login nodes, such as submitting/monitoring jobs, preparing scripts, compling small applications, and transferring small amount of data etc. All other types of workloads are encouraged to be submitted as a job. 

To ensure fair usage of login nodes, any process that use more than 30 minutes CPU time or instantly use more than 4 GiB memory on any login node will be killed immediately. 

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



