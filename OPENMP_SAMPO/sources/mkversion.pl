#!/usr/bin/perl
# -*- Mode: cperl -*-
############################################################################
# Makes the version file Version.F90
############################################################################
#    Copyright (C) 2003 Olivier Cessenat.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
############################################################################
use Sys::Hostname; 
use File::Basename ;
use Getopt::Std;
getopts('c:d:h:i:v:t:') ;
$opt_c = '.' unless $opt_c ;
$filename = "Version.F90" ;
if (-f "$opt_c/$filename") {
  ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat("$opt_c/$filename");
} else {
  $mtime = 0 ;
}
$thetime = time() ;
(@days) = ("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat") ;
(@months) = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec") ;
if (abs($mtime-$thetime) < 120) {
  print STDOUT "\tNew request automatic: mtime=$mtime, thetime=$thetime\n" ;
  exit ;
} else {
  print STDOUT "\tRebuilding $opt_c/$filename: mtime=$mtime, thetime=$thetime\n" ;
}
unless ($opt_t) {
  $opt_t = "HEAD" ;
}
unless ($opt_c) {
  $opt_c = "." ;
}
############################################################################
# Guesses the current CVS/SVN version:
sub getVersion {
  my $version = "" ;
  my $tag = "" ;
  my $date = "" ;
  my(@versions) = () ;
  my $entryfile = "../src/CVS/Entries" ;
  unless (-e "$entryfile") {
    $entryfile = "CVS/Entries" ;
  }
  # print STDOUT "Reading CVS Entries in $entryfile\n" ;
  if (open(CVS,"$entryfile")) {
    while(<CVS>) {
      next unless (/\//) ;
      my (@a) = split(/\//, $_ ) ;
      my $fic = $a[1] ;
      my $ver = $a[2] ;
      my $d = $a[3] ;
      if ($date) {
	$date = mergeDates($d, $date) ;
      } else {
	$date = $d ;
      }
      my $t = $a[5] ;
      if ($t) {
	if ($t =~ /T(.*)/) {
	  $t = $1 ;
	  if ($tag) {
	    $tag = max($tag, $t) ;
	  } else {
	    $tag = $t ;
	  }
	}
      }
      if ($ver) {
	my (@b) = split(/\./, $ver) ;
	my $l = 0 ;
	foreach my $s (@b) {
	  $versions[$l] += $s ;
	  $l++ ;
	}
      }
    }
    close(CVS) ;
    if (@versions) {
      # Try to make it a shorter string, remove duplicates.
      # $version = join("_", @versions) ;
      $version = $versions[0] ;
      for ($i=1;$i<=$#versions;$i++) {
	if ($versions[$i] != $versions[($i-1)]) {
	  $version .= '_' . $versions[$i] ;
	}
      }
      print STDOUT "Versions = $version\n" if $debug ;
    }
  } else { 
    $entryfile = "../src/.svn/entries" ;
    unless (-e "$entryfile") {
      $entryfile = ".svn/entries" ;
    }
    if (open(SVN,"$entryfile")) {
      my $td ;
      while(<SVN>) {
	if ($. == 11) {
	  chomp ;
	  $version = $_ ;
	  $svn = $_ ;
	  print STDOUT "SVN version=$version\n" if $debug ;
	}
	# Let us get the current branch or tag:
	if ($. == 5) {
	  # indicates repository location:
	  if (/\b(branch|tag)e*s*\/(\w+)/) {
	    if ($1 eq "tag") {
	      $tag = "V_" . $2 ;
	    } else {
	      $tag = "B_" . $2 ;
	    }
	  }
	}
	if ($. == 6) {
	  if (/^file:\/\/(.*)/) {
	    $td = $1 ;
	  }
	}
	if ($. == 10) {
	  if (/(\d+)\-(\d+)\-(\d+)T([\d\:]+)[\.\d]*Z/) {
	    my (@months) = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec") ;
	    # print STDOUT "time stamp is $1,$2,$3 :: $4\n" ;
	    my $m = ( $2 - 1 ) ;
	    $date = "$months[$m] $3 $4 $1" ;
	    print STDOUT "date=$date\n" if $debug ;
	  }
	}
      }
      close(SVN) ;
      print STDOUT "tag=$tag, version=$version\n" ;
      if ($version) {
	# Now check that the SVN number is not modified:
	if ($td && ! -d "$td") {
	  print STDOUT "No SVN dir: $td\n" if $debug ;
	  if (open(SYS,"<mkobj.lst")) {
	    while(<SYS>) {
	      if (/VERSION\s*\=\s*(.*)/) {
		$version = $1 ;
		print STDOUT "Reusing version: $1\n" if $debug ;
	      }
	    }
	    close(SYS) ;
	  }
	} else {
	  # Standard check using SVN server:	
	  print STDOUT "Using SVN dir\n" if $debug ;
	  $remsh = "svn diff -r $svn" ;
	  if (open(SYS, "$remsh 2>&1 |")) {
	    while(<SYS>) {
	      print ;
	      $version .= "modif" ;
	      last ;
	    }
	    close(SYS) ;
	  }
	}
      }
      $tag = 'TRUNK' unless $tag ;
    } else {
      print STDERR "No CVS/SVN version\n" ;
    }
  }
  $tag = 'HEAD' unless $tag ;
  return($version,$date,$tag)
}
############################################################################
unless ($opt_d) {
  ($opt_v,$opt_d,$opt_t) = getVersion() ;
}
print STDOUT "\tOK: Version = '$opt_v', Library = '$opt_h', Infos='$opt_i', date='$opt_d', tag='$opt_t', write='$opt_c'\n" ;
############################################################################
my $d = "" ;
if ($opt_d) {
  $d = $opt_d ;
  if ($d =~ /^\s*(\w+)\s+(\w+)\s+(\d+)\s+([\d:]+)\s+(\d+)/) {
    $d = $1.$2.$3.",$4,$5"
  }
}
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
#    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
############################################################################
#print STDOUT "Output on $opt_c/$filename\n" ;
unless (open(VERSION, ">$opt_c/$filename")) {
  if (-d "$opt_c") {
    die "Unable to write $filename file" ;
  } else {
    die "Unable to write $filename file, $opt_c does not exist" ;
  }
}
############################################################################
# File write:
print VERSION 'module Version' . "\n" ;
print VERSION '  implicit none' . "\n" ;
#print VERSION '  include "Version.h"' . "\n" ;
print VERSION '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' . "\n" ;
print VERSION '  interface getVersion' . "\n" ;
print VERSION '     module procedure getVersionVersion' . "\n" ;
print VERSION '  end interface' . "\n" ;
print VERSION '  interface getCVS' . "\n" ;
print VERSION '     module procedure getVersionCVS' . "\n" ;
print VERSION '  end interface' . "\n" ;
print VERSION '  interface getIdent' . "\n" ;
print VERSION '     module procedure getVersionIdent' . "\n" ;
print VERSION '  end interface' . "\n" ;
print VERSION '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' . "\n" ;
print VERSION 'contains' . "\n" ;
print VERSION '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' . "\n" ;
print VERSION '  function getVersionIdent() result(myid)' . "\n" ;
print VERSION '    implicit none' . "\n" ;
print VERSION '    integer :: myid' . "\n" ;
if ("$opt_t" ne "HEAD") {
  # Provides the locked version number:
  if ($opt_t =~ /.*?(\d+)$/) {
    print VERSION "    myid = $1\n" ;
  } else {
    # Provides the today's version number for probabibly unofficial release:
    print VERSION sprintf "    myid = %2.2d%2.2d%2.2d\n", ($year-100), ($mon+1), $mday ;
  }
} else {
  # Provides the today's version number:
  print VERSION sprintf "    myid = %2.2d%2.2d%2.2d\n", ($year-100), ($mon+1), $mday ;
}
print VERSION '  end function getVersionIdent' . "\n" ;
print VERSION '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' . "\n" ;
print VERSION '  subroutine getVersionCVS(name)' . "\n" ;
print VERSION '    implicit none' . "\n" ;
print VERSION '    character(128) :: name' . "\n" ;
if ($opt_t) {
  print VERSION '    name = ' . "'$opt_t'" . "\n" ;
}
print VERSION '  end subroutine getVersionCVS' . "\n" ;
print VERSION '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' . "\n" ;
print VERSION '  function getVersionVersion() result(name)' . "\n" ;
print VERSION '    implicit none' . "\n" ;
print VERSION '    character(256) :: name' . "\n" ;
$hostname = hostname() ;
if ($hostname =~ /^(.*?)\d+/) {
  $hostname = $1 ;
}
$fci = "" ;
if ($opt_i) {
  $remsh = "$opt_i" ;
  if (open(SYS, "$remsh 2>&1 |")) {
    while(<SYS>) {
      chomp ;
      if (/^f90:\s*Forte Developer .*? Fortran ([\d\s\.]*) Patch (.*)/) { 
	$fci = "Forte $1:$2" ;
      } elsif (/^f90:(.*)/) {
	$fci = $1 ;
      } elsif (/Version\s+(.*?)\s+Build\s+(\d+).*:\s+(.*)/) {
	$fci = "${1}_${2}:$3" ;
      } elsif (/Version\s+(.*?)\s+Build\s+(\d+)/) {
	# Intel V12:
	$fci = "${1}:${2}" ;
      } elsif (/Fortran Compiler\s+(.*)/) {
	$fci = $1 ;
      } elsif (/GNU Fortran \d+ \(\s*(.*?)\s*([\(\)])\s*(.*?)\s*/) {
	$fci = $1 ;
	if ($2 eq ')') {
	  $fci = $3 ;
	}
      } elsif (/GNU Fortran \d+ \(\w+\) ([\s\d\.]+) /) {
	$fci = $1.$2 ;
      } elsif (/GNU Fortran \(\w+\) ([\s\d\.]+) /) {
	$fci = $1 ;
	$fci =~ tr/ /@/ ;
      } elsif (/GNU gcc version.*?\(\s*(.*)\s*\)/) {
	$fci = $1 ;
      }
    }
  }
}
print VERSION "    name = '" . "$opt_h($hostname)[$fci]" . "'\n" ;
print VERSION '#ifdef __OPEN64__' . "\n" ;
print VERSION '#define trim(name) name(1:len_trim(name))' . "\n" ;
print VERSION '#endif' . "\n" ;
print VERSION '#if defined __alpha' . "\n" ;
print VERSION '    name = trim(name) // "->alpha"' . "\n" ;
print VERSION '#elif defined __linux__' . "\n" ;
print VERSION '    name = trim(name) // "->linux"' . "\n" ;
#print VERSION '    name = name(1:len_trim(name)) // "->linux"' . "\n" ;
print VERSION '#elif defined __sun' . "\n" ;
print VERSION '    name = trim(name) // "->sun"' . "\n" ;
print VERSION '#elif defined __sgi' . "\n" ;
print VERSION '    name = trim(name) // "->sgi"' . "\n" ;
print VERSION '#elif defined __hpux' . "\n" ;
print VERSION '    name = trim(name) // "->hpux"' . "\n" ;
print VERSION '#elif ( __APPLE__ || macintosh )' . "\n" ;
print VERSION '#if defined __GNUC__' . "\n" ;
print VERSION '    name = trim(name) // "->macosx"' . "\n" ;
print VERSION '#else' . "\n" ;
print VERSION '    name = trim(name) // "->macos"' . "\n" ;
print VERSION '#endif' . "\n" ;
print VERSION '#elif defined AIX' . "\n" ;
print VERSION '    name = trim(name) // "->aix"' . "\n" ;
print VERSION '#elif defined WINNT' . "\n" ;
print VERSION '    name = trim(name) // "->dos"' . "\n" ;
print VERSION '#elif ( __GNUC__ || __gnu_linux__ )' . "\n" ;
print VERSION '    name = trim(name) // "->gnu"' . "\n" ;
print VERSION '#else' . "\n" ;
print VERSION '    name = trim(name) // "->auto"' . "\n" ;
print VERSION '#endif' . "\n" ;
print VERSION '#ifdef __ia64__' . "\n" ;
print VERSION '    name = trim(name) // "ia64"' . "\n" ;
print VERSION '#elif __IA64__' . "\n" ;
print VERSION '    name = trim(name) // "IA64"' . "\n" ;
print VERSION '#elif __x86_64__' . "\n" ;
print VERSION '    name = trim(name) // "x86_64"' . "\n" ;
print VERSION '#elif __arch64__' . "\n" ;
print VERSION '    name = trim(name) // "64"' . "\n" ;
print VERSION '#elif defined __arch32__' . "\n" ;
print VERSION '    name = trim(name) // "32"' . "\n" ;
print VERSION '#endif' . "\n" ;
# Major compilation flags settings:
print VERSION '#ifdef DEBUG_CODE' . "\n" ;
print VERSION '    name = trim(name) // "D"' . "\n" ;
print VERSION '#else' . "\n" ;
print VERSION '    name = trim(name) // "d"' . "\n" ;
print VERSION '#endif' . "\n" ;
print VERSION '#ifdef LIB_USER' . "\n" ;
print VERSION '    name = trim(name) // "C"' . "\n" ;
print VERSION '#else' . "\n" ;
print VERSION '    name = trim(name) // "c"' . "\n" ;
print VERSION '#endif' . "\n" ;
print VERSION '#ifdef MKL_ACTIF' . "\n" ;
print VERSION '    name = trim(name) // "M"' . "\n" ;
print VERSION '#else' . "\n" ;
print VERSION '    name = trim(name) // "m"' . "\n" ;
print VERSION '#endif' . "\n" ;
print VERSION '#ifdef PETSC_ACTIF' . "\n" ;
print VERSION '    name = trim(name) // "P"' . "\n" ;
print VERSION '#else' . "\n" ;
print VERSION '    name = trim(name) // "p"' . "\n" ;
print VERSION '#endif' . "\n" ;

if ($opt_v) {
  print VERSION '    name = trim(name) // " Version=" // ' . "'$opt_v'" . "\n" ;
} else {
  print VERSION '    name = trim(name) // " (no SVN)"' . "\n" ;
}
# Last date of commit is:
if ($opt_d) {
  print VERSION '    name = trim(name) // " SVN:" // ' . "'$d'" . "\n" ;
}
if ($opt_t) {
  print VERSION '    name = trim(name) // " Tag:" // ' . "'$opt_t'" . "\n" ;
}
$year += 1900 ;
$day = $days[$wday] ;
$month = $months[$mon] ;
$mon++ ;
if (0) {
    print VERSION '#if defined __DATE__' . "\n" ;
    print VERSION '    name = trim(name) // __DATE__' . "\n" ;
    print VERSION '#endif' . "\n" ;
    print VERSION '#if defined __TIME__' . "\n" ;
    print VERSION '    name = trim(name) // __TIME__' . "\n" ;
    print VERSION '#endif' . "\n" ;
} else {
    print VERSION '!|    name = trim(name) // ' . "'$thetime'" . "\n" ;
    print VERSION '    name = trim(name) // " Compile:" // ' . "'$day$month$mday,$hour:$min:$sec,$year'" . "\n" ;
}
print VERSION '  end function getVersionVersion' . "\n" ;
print VERSION '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' . "\n" ;
print VERSION 'end module Version' . "\n" ;
############################################################################
close(VERSION);
