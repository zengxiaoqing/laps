/*cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis  
cdis 
cdis*/
From kbrews@cirrus.gcn.uoknor.edu Thu Apr 20 20:56:01 1995
Return-Path: <kbrews@cirrus.gcn.uoknor.edu>
Received: from cirrus.gcn.uoknor.edu by peaks.fsl.noaa.gov (4.1/SMI-4.1)
	id AA29408; Thu, 20 Apr 95 20:55:59 DST
Received: by cirrus.gcn.uoknor.edu (AIX 3.2/UCB 5.64/4.03)
          id AA14571; Thu, 20 Apr 1995 21:56:10 -0500
From: kbrews@cirrus.gcn.uoknor.edu (Keith Brewster)
Message-Id: <9504210256.AA14571@cirrus.gcn.uoknor.edu>
Subject: Re: command line FTP script
To: albers@peaks.fsl.noaa.gov (Steve Albers)
Date: Thu, 20 Apr 1995 21:56:10 +22306039 (CDT)
In-Reply-To: <9504210132.AA21203@peaks.fsl.noaa.gov> from "Steve Albers" at Apr 21, 95 01:32:22 am
X-Mailer: ELM [version 2.4 PL21]
Content-Type: text
Content-Length: 1782      
Status: RO

> 
>      Could you please resend me that command line format for LAPS_XMIT
>   that you have. I think now I might be able to modify that ftp script I
>   had to use the command line format you had - using the $1, $2, etc. as
>   I had mentioned.

Here's the comment block from the code.


/* Create the custom user command

   For this option, a command line is created from the
   LAPS environment variables and the long and short filenames
   The LAPS_XMIT variable (represented here by xmit_ptr) is the
   name of a user-programmed script or executable that accepts
   the arguments as passed here.  These arguments should give
   the script or program enough info to work with.

   The command line created consists of:
   LAPS_XMIT long_filename short_filename LAPS_USER LAPS_REMOTE LAPS_DEST

 */


I have a sample script which uses this command line to create an
rcp command - see below.   As I mentioned I tested this only once
and it didn't work.  Haven't had a chance to investigate.

-Keith

#!/bin/csh -f
#
# custom_example
#
#
# Demonstrates the interface between remap.c and
# and user-defined script for data management
#
# This one does rcp, which is supported within remap.c
# if LAPS_XMIT has been set to rcp.
#
# Keith Brewster
# CAPS, April, 1995
#

#
# FIRST: get command line args
#
if ( $#argv > 4 ) then
   set longfilnam = $argv[1]
   set shrtfilnam = $argv[2]
   set usernam = $argv[3]
   set hostnam = $argv[4]
   set destnam = $argv[5]
   echo rcp $longfilnam $usernam"@"$hostnam":"$destnam"/"$shrtfilnam
# do it
   rcp $longfilnam $usernam"@"$hostnam":"$destnam"/"$shrtfilnam

else
   echo " "
   echo Error in custom_example
   echo Received only $#argv arguments on command line from remap
   echo Check LAPS_ environment variables
   echo " "
endif
