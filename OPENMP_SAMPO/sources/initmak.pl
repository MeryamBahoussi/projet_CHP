#!/usr/bin/perl
# -*- Mode: cperl -*-
############################################################################
# A program to build fortran 90 module dependencies for allthe files of current
# directory.
# Guesses the best order for making build, always putting
# dependencies before so that the recompilation of a module and one
# of it's dependencies always works if no side effect.
# Should also work with C++.
# Excludes main from processing, creating a specific rule that suits the local
# makefile file in current directory.
# Usage:
#   perl $0 -d => generates all dependencies
#   perl $0    => does not generate dependencies for faster recompilation.
#   perl $0 -f => generates dependencies for the main program only, use with caution.
#   perl $0 -d -f => generates all dependencies + program's global dependencies.
############################################################################
#  -t string : link towards ftagshtml restart file (typically ./HTML/ftagshtml/ftagshtml.pm)
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
use File::Basename ;
use Cwd ;
use Getopt::Std;
getopts('A:bdfI:lL:Rs:t:V') ;
$debug = 0 ;
if (@ARGV) {
  (@dirs) = (@ARGV) ;
} else {
  (@dirs) = ('.') ;
}
$lasttime = 0 ;
# Put here the list of system uses:
(%sysuse) = ("omp_lib" => 1, "mpi" => 1, "version" => 1, "ifport" => 1, "mkl_dss" => 1, "mod_hic_f90" => 1,) ;
(@aincludes) = split(':', $opt_A) ;
(@rincludes) = split(':', $opt_I) ;
if ($opt_l) {
  $lib1 = '$(LIB)(' ;
  $lib2 = ')' ;
} else {
  $lib1 = "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/" ;
  $lib2 = '' ;
}
$mkversion = "../src/mkversion.pl" ;
unless (-e "$mkversion") {
  if (-e "./mkversion.pl") {
    $mkversion = "./mkversion.pl" ;
  #} ) {
  #  ;
  } else {
    $mkversion = "mkversion.pl" ;
  }
}
$md5 = "" ;
if (1 || -f "$mkversion") {
  $remsh = "pwd | md5sum" ;
  if (open(SYS, "$remsh 2>&1 |")) {
    while(<SYS>) {
      $md5 = "/$1" if (/^(\w+)\s*/) ;
      print ;
    }
  }
}
$reuse = 0 ;
if ($opt_t) {
  if (-s $opt_t) {
    # Load a dependencies analysis file !
    require ("$opt_t") ;
    $reuse = 1 ;
    print STDOUT "Reusing file $opt_t\n" if $debug ;
  } else {
    print STDOUT "Warning: File $opt_t does not exist, regenerate using fth.pl -code=sophie *.F90 *.f90 *.f *.F *.c\n" ;
  }
}
############################################################################
# Gets the absolute path to $2 from the $1 location: path becomes valid
# from the current working directory
sub getAbs {
  my $from = shift ;
  my $to = shift ;
  my $out ;
  my (@out) = split('/', $from) ;
  if ($to =~ /^\// ) {
    $out = $to ;
  } else {
    foreach my $l (@out) {
      if ($to =~ /^\.\.\/(.*)/ ) {
	$to = $1 ;
      }
    }
    $out = $to ;
  }
}
############################################################################
# Gets the relative location to $2 from $1
sub getRel {
  my $from = shift ;
  my $to = shift ;
  my $out ;
  if ($to =~ /^\// ) {
    $out = $to ;
  } else {
    # Tries to get ../$2 or something like this from $base...
    if ($from ne '.') {
      $to = "../$to" ;
    }
    while (dirname($from) ne '.') {
      $to = "../$to" ;
    }
    $out = $to ;
  }
}
############################################################################
# Tries to get the include file actual location according to the -I variables.
sub getInc {
  my $base = shift ;
  my $in  = shift ;
  my $out = shift ;
  if ((@aincludes)||(@rincludes)) {
    foreach $i (@aincludes) {
      if (-e "$i/$in") {
	$out = &getRel($base,"$i/$in") ;
      }
    }
    foreach $i (@rincludes) {
      my $j = $getrel{$i} ;
      if (-e "$j/$in") {
	$out = "$i/$in" ;
      }
    }
  } else {
    if (-e "$base/$in") {
      $out = $in ;
    }
  }
  return($out) ;
}
############################################################################
# Looks for include calls in CPP files
sub getDependcpp {
  my $f = shift ;
  my $s = shift ;
  my $d = shift ;
  $isamain{$s} = 0 ;
  if (open(FILE,"<$f")) {
    $lasttime = max($lasttime, (stat("$f"))[9]) ;
    while(<FILE>) {
      # eliminates some comments, not all...
      next if (/^\s*\/\//) ;
      if (/^\s*\#include\s*([\"\'])(.*)\1/) {
	$inc = &getInc($d,$2);
	push(@{$dependencies{$s}},$inc) ;
      } elsif (/^\s*(int|void)\s*main\s*\(/) {
	# A not very smart manner to guess for main...
	$isamain{$s} = 1 ;
      } elsif (/^\s*CALL_FORTRAN\s*\(\s*(\w+)\s*\)/) {
	# Says we call a fortran stuff:
	$_ = "$1(" ;
	&parseCalls($s) ;
	#$r = lc($1) ;
	#$calls{$s}->{$r} = 1 ;
      } elsif ($reuse && /(\w+)\s*\((.*)/) {
	&parseCalls($s) ;
      }
    }
    close(FILE) ;
  }
}
############################################################################
# Looks for include calls in C files
sub getDependc {
  my $f = shift ;
  my $s = shift ;
  my $d = shift ;
  $isamain{$s} = 0 ;
  if (open(FILE,"<$f")) {
    $lasttime = max($lasttime, (stat("$f"))[9]) ;
    while(<FILE>) {
      if (/^\s*\#include\s*([\"\'])(.*)\1/) {
	$inc = &getInc($d,$2);
	push(@{$dependencies{$s}},$inc) ;
      } elsif (/^\s*(int|void)\s*main\s*\(/) {
	$isamain{$s} = 1 ;
      } elsif (/^\s*CALL_FORTRAN\s*\(\s*(\w+)\s*\)/) {
	# Says we call a fortran stuff:
	$_ = "$1(" ;
	&parseCalls($s) ;
	print STDOUT "Added parse for '$1('\n" if $debug ;
      } elsif ($reuse && /(\w+)\s*\((.*)/) {
	&parseCalls($s) ;
      }
    }
    close(FILE) ;
  }
}
############################################################################
# Looks for include calls in fortran 77 files
sub getDependf {
  my $f = shift ;
  my $s = shift ;
  my $d = shift ;
  $isamain{$s} = 0 ;
  if (open(FILE,"<$f")) {
    $lasttime = max($lasttime, (stat("$f"))[9]) ;
    while(<FILE>) {
      if (/^\s*program\s*(\w+)/i) {
	$isamain{$s} = 1 ;
      } elsif (/^\s*subroutine\s+mpc_user_main\b/i) {
	$isamain{$s} = 1 ;
      } elsif (/^\s{6}\s*include\s+([\"\'])(.*)\1/i) {
	$inc = &getInc($d,$2);
	push(@{$dependencies{$s}},$inc) ;
      } elsif (/^\s{6}\s*use\s+(\w+)\b/i) {
	# Hum, not really safe, case management...
	my $r = lc($1) ;
	unless (exists $sysuse{$r}) {
	  $uses{$s}->{$r} = 1 ;
	}
      } elsif (/^\s{6}\s*module\s*(\w+)/i) {
	# Provides
	my $r = lc($1) ;
	$modules{$r} = "$s" ;
	$revmods{$s} = "$r" ;
	print STDOUT "File $s provides module $1\n" if $debug ;
      } elsif ($reuse && /(\w+)\s*\((.*)/) {
	&parseCalls($s) ;
      }
    }
    close(FILE) ;
  }
}
############################################################################
# Looks for include calls in fortran 90 files
sub getDependf90 {
  my $f = shift ;
  my $s = shift ;
  my $d = shift ;
  $isamain{$s} = 0 ;
  if (open(FILE,"<$f")) {
    $lasttime = max($lasttime, (stat("$f"))[9]) ;
    while(<FILE>) {
      next if (/^\s*\#/) ;
      if (/^\s*program\s+(\w+)/i) {
	$isamain{$s} = 1 ;
	print STDOUT "Adding a main $s\n" if $debug ;
      } elsif (/^\s*subroutine\s+mpc_user_main\b/i) {
	$isamain{$s} = 1 ;
      } elsif (/^\s*include\s+([\"\'])(.*)\1/i) {
	$inc = &getInc($d,$2);
	push(@{$dependencies{$s}},$inc) ;
      } elsif (/^\s*use iso_c_binding/i) {
      } elsif (/^\s*use petscvec/i) {
      } elsif (/^\s*use petscmat/i) {
      } elsif (/^\s*use petscksp/i) {
      } elsif (/^\s*use petscsys/i) {
      } elsif (/^\s*use MOD_HIC_F90/i) {
      } elsif (/^\s*use\s+(\w+)\b/i) {
	# Hum, not really safe, case management...
	my $r = lc($1) ;
	unless (exists $sysuse{$r}) {
	  $uses{$s}->{$r} = 1 ;
	}
	print STDOUT "\t+file $s uses module $r line $.\n" if $debug ;
      } elsif (/^\s*module procedure\s*(\w+)/i) {
      } elsif (/^\s*module\s*(\w+)/i) {
	# Provides
	my $r = lc($1) ;
	if (exists $modules{$r}) {
	  # Module  is already defined
	  print STDERR "Already defined module $r in $modules{$r} newly defined in $f: given up file $f\n" ;
	  $eliminate{$s} = 1 ;
	  last ;
	}
	$modules{$r} = "$s" ;
	$revmods{$s} = "$r" ;
	print STDOUT "File $s provides module $r\n" if $debug ;
      } elsif ($reuse && /(\w+)\s*\((.*)/) {
	&parseCalls($s) ;
      }
    }
    close(FILE) ;
  }
}
############################################################################
sub parseCalls {
  my $s = shift ;
  while(/(\w+)\s*\((.*)/) {
    $_ = $2 ;
    my $r = lc($1) ;
    print STDOUT "\tr=$r\n" if $debug ;
    if (exists $names{$r}) {
      my $p = $names{$r} ;
      print STDOUT "\t\tp=$p\n" if $debug ;
      if (exists $fthobjs{$p}) {
	my (@t) = (@{$fthobjs{$p}}) ;
	my (@f) = split(/:/, $t[0]) ;
	print STDOUT "Adding module call to file $f[0] from stuff $s:$.:$r\n" if $debug ;
	$calls{$s}->{$f[0]} = 1 ;
      } elsif (exists $fthrouts{$p}) {
	my (@f) = split(/:/, $fthrouts{$p}) ;
	print STDOUT "Adding simple call to file $f[0] from stuff $s:$.:$r\n" if $debug ;
	$calls{$s}->{$f[0]} = 1 ;
      }
    }
  }
}
############################################################################
# Tries to guess the best sorting for the compilation order in mkobj.lst
sub trySort {
  my(@oldfiles) = (@_) ;
  print STDOUT "Entering trySort with arguments = @oldfiles\n" if $debug ;
  my(@newfiles) = () ;
  my(@tmpfiles) = () ;
  # Get the files with dependencies:
  my $remains = 0 ;
  # my $debug = 1 ;

  foreach (@oldfiles) {
    if (exists $uses{$_}) {
      push(@tmpfiles, $_) ;
      $remains++ ;
    } else {
      push(@newfiles, $_) ;
      $dones{$_} = 1 ;
    }
  }
  my $ier = 0 ;
  my(%errs) = () ;
  my $level = 1 ;
  $noinfty = 1 ;
  while($remains && $noinfty) {
    print STDOUT "trySort: Level $level tmpfiles=@tmpfiles\n" if $debug ;
    print STDOUT "trySort: Level $level newfiles=@newfiles\n" if $debug ;
    $noinfty = 0 ;
    (@oldfiles) = (@tmpfiles) ;
    (@tmpfiles) = () ;
    $remains = 0 ;
    $level++ ;
    print STDOUT "level = $level @oldfiles\n" if $debug ;
    if ($ier) {
      print STDERR "At least one of the following modules are not set:\n" ;
      foreach my $module (keys %errs) {
	print STDOUT "\t$module, required by at least $errs{$module}\n" ;
      }
      print STDERR "initmak.pl: KO" ;
    }
    foreach (@oldfiles) {
      my $ko = 0 ;
      print STDOUT "Analysis for $_:\n" if $debug ;
      foreach my $module (keys %{$uses{$_}}) {
	if (exists $modules{$module}) {
	  my $file = $modules{$module} ;
	  print STDOUT "\tUsed module for $_ is $module defined in $file.$types{$file}\n"  if $debug;
	  my $dep = $modules{$module} ; # This is a file
	  if (exists $dones{$dep}) {
	    print STDOUT "\tDepend for $module is $dep: already done\n" if $debug ;
	  } elsif ($dep eq $modules{$_} || $dep eq $_) {
	    print STDOUT "\tDepend for $module is $dep: equals self\n" if $debug ;
	  } else {
	    print STDOUT "\tDepend for $module is $dep: not already done (from $module, not equal $_, ie $modules{$_})\n" if $debug ;
	    $ko = 1 ;
	    last ;
	  }
	} elsif (!exists $errs{$module}) {
	  print STDOUT "\tModule $module required by $_ not known\n" if $debug ;
	  $errs{$module} = "$_" ;
	  $ier++ ;
	}
      }
      if ($ko) {
	push(@tmpfiles, $_) ;
	if ($debug) {
	  my(@totos) = (keys %{$uses{$_}}) ;
	  print STDOUT "\t$_: uses list (@totos) not found in dones\n" if $debug;
	}
	$remains++ ;
     } else {
	push(@newfiles, $_) ;
	$dones{$_} = $level ;
	print STDOUT "\tDone for $_ at level $level is performed\n" if $debug ;
	$noinfty = 1 ;
      }
    }
  }
  if (!$noinfty) {
    print STDERR "trySort(*): Error on level = $level, $#oldfiles: @oldfiles\n";
    print STDERR "tmpfiles $#tmpfiles: @tmpfiles\n";
    print STDERR "remains is $remains\n";
    die "initmak.pl: KO" ;
  } elsif ($debug) {
    print STDOUT "trySort(W): Error on level = $level, $#oldfiles: @oldfiles\n";
    print STDOUT "tmpfiles $#tmpfiles: @tmpfiles\n";
    print STDOUT "remains is $remains\n";
  }
  # Now returns files sorted by level then by priority if any then by name:
  return(@newfiles) ;
}
############################################################################
sub bylevel {
  # First lower level functions then by alphabetical order
  $dones{$a} <=> $dones{$b} || $a cmp $b ;
}
############################################################################
# Flat by file sort:
sub byname {
  my $a1 = $a ;
  my $b1 = $b ;
  my $a2 = $a ;
  my $b2 = $b ;
  if ($a =~ /(.*)\/(.*)/) {
    $a1 = $2 ; $a2 = $1 ;
  }
  if ($b =~ /(.*)\/(.*)/) {
    $b1 = $2 ; $b2 = $1 ;
  }
  $a1 cmp $b1 || $a2 cmp $b2 ;
}

############################################################################
# A more complicated routine that sorts the modules by level then by
# one level dependency, this is not very safe but could be sufficient, it
# depends on your code hierarchy:
sub glolevel {
  my $xa = 0 ;
  my $xb = 0 ;
  if ($dones{$a} == $dones{$b}) {
    # The following unsafe stuff is costly, run it only if same level:
    my(%depsa)=(%{$uses{$a}}) ;
    my(%depsb)=(%{$uses{$b}}) ;
    foreach my $moda (keys %depsa) {
      if ($modules{$moda} eq $b) {
	$xa = 1 ; 
      }
    }
    foreach my $modb (keys %depsb) {
      if ($modules{$modb} eq $a) {
	$xb = 1 ; 
      }
    }
  }
  $xa <=> $xb || $dones{$a} <=> $dones{$b} || $a cmp $b ;
}
############################################################################
# A conservative one level only sorting
sub constlevel {
  my $xa = 0 ;
  my $xb = 0 ;
  my $ldebug = 0 ;
  if (1 == 0 && ($a eq 'sFDTD' || $b eq 'sFDTD')) {
    $ldebug = 1 ;
    print STDOUT "Comparing $a with $b\n" if $ldebug ;
    if (exists $uses{$a}) {
      my(%depsa)=(%{$uses{$a}}) ;
      foreach my $moda (keys %depsa) {
	print STDOUT " $moda" ;
      }
    }
    print STDOUT "\n" ;
  }
  # my $sa = $a ;
  # my $sb = $b ;
  # $sa =~ s/\.\w+$// ;
  # $sb =~ s/\.\w+$// ;
  if (exists $uses{$a}) {
    # The following unsafe stuff is costly, run it only if same level:
    my(%depsa)=(%{$uses{$a}}) ;
    foreach my $moda (keys %depsa) {
      if ($modules{$moda} eq $b) {
	# Means a uses module b, so b is first
	print STDOUT "Resorting(1), put $a after $b, since $b used by $a\n" if $ldebug ;
	$xa = 1 ; 
	last ;
      }
    }
  }
  if (exists $uses{$b}) {
    my(%depsb)=(%{$uses{$b}}) ;
    foreach my $modb (keys %depsb) {
      if ($modules{$modb} eq $a) {
	# Means b uses module a, so a is first
	print STDOUT "Resorting(2), put $b after $a, since $a used by $b\n" if $ldebug ;
	$xb = 1 ; 
	last ;
      }
    }
  }
  if ($xa == 0 && $xb == 0) {
    # Make a second level search:
    if (exists $uses{$a}) {
      # The following unsafe stuff is costly, run it only if same level:
      my(%depsa1)=(%{$uses{$a}}) ;
      foreach my $moda1 (keys %depsa1) {
	my $depa1 = $modules{$moda1} ;
	if (exists $uses{$depa1}) {
	  my(%depsa2)=(%{$uses{$depa1}}) ;
	  foreach my $moda2 (keys %depsa2) {
	    if ($modules{$moda2} eq $b) {
	      # Means a uses module b, so b is first
	      print STDOUT "Resorting(3), put $a after $b, since $b used by $a\n" if $ldebug ;
	      $xa = 1 ;
	      last ;
	    }
	  }
	  last if $xa ;
	}
      }
    }
    if (exists $uses{$b}) {
      # The following unsafe stuff is costly, run it only if same level:
      my(%depsb1)=(%{$uses{$b}}) ;
      foreach my $modb1 (keys %depsb1) {
	my $depb1 = $modules{$modb1} ;
	if (exists $uses{$depb1}) {
	  my(%depsb2)=(%{$uses{$depb1}}) ;
	  foreach my $modb2 (keys %depsb2) {
	    if ($modules{$modb2} eq $a) {
	      # Means b uses module a, so a is first
	      print STDOUT "Resorting(4), put $b after $a, since $a used by $b\n" if $ldebug ;
	      $xb = 1 ;
	      last ;
	    }
	  }
	  last if $xb ;
	}
      }
    }
  }
  if ($xa == 0 && $xb == 0) {
    # Look for b, stuffs being called by a
    if (exists $calls{$a}) {
      my(%depsa)=(%{$calls{$a}}) ;
      foreach my $f (keys %depsa) {
	my $moda = $f ;
	$moda =~ s/\.\w+$// ;
	if ($moda eq $b) {
	  $moda = lc($moda) ;
	  unless (exists $modules{$moda}) {
	    # Means a calls module b, so b is first
	    print STDOUT "Resorting(5), put $b before $a, since called 1 ($moda not a module)\n" if $ldebug ;
	    $xa = 1 ; 
	    last ;
	  }
	}
      }
    }
    if (exists $calls{$b}) {
      my(%depsb)=(%{$calls{$b}}) ;
      foreach my $f (keys %depsb) {
	my $modb = $f ;
	$modb =~ s/\.\w+$// ;
	if ($modb eq $a) {
 	  $modb = lc($modb) ;
	  unless (exists $modules{$modb}) {
	    # Means b calls module a, so a is first
	    print STDOUT "Resorting(6), put $a before $b, since called 1 ($modb not a module)\n" if $ldebug ;
	    $xb = 1 ; 
	    last ;
	  }
	}
      }
    }
  }
  if (1 == 0 && $xa == 0 && $xb == 0) {
    # Considered as probably too dangerous:
    if (!exists $uses{$a} && exists $uses{$b}) {
      # Probably that a is before b
      print STDOUT "Resorting(7), put $a before $b, since $a does not use\n" if $ldebug ;
      $xa = 1 ; 
    }
    if (exists $uses{$a} && !exists $uses{$b}) {
      # Probably that b is before a
      print STDOUT "Resorting(8), put $b before $a, since $b does not use\n" if $ldebug ;
      $xb = 1 ; 
    }
  }
  if ($xa == 0 && $xb == 0) {
    print STDOUT "Unable to sort, leave $a before $b\n" if $ldebug ;
  }
  $xa <=> $xb || 0 ;
  # $xb <=> $xa || 0 ;
}
############################################################################
sub basiclevel {
  my $xa = 0 ;
  my $xb = 0 ;
  # The following unsafe stuff is costly, run it only if same level:
  my(%depsa)=(%{$uses{$a}}) ;
  my(%depsb)=(%{$uses{$b}}) ;
  foreach my $moda (keys %depsa) {
    if ($modules{$moda} eq $b) {
      $xa = 1 ; 
    }
  }
  foreach my $modb (keys %depsb) {
    if ($modules{$modb} eq $a) {
      $xb = 1 ; 
    }
  }
  $xa <=> $xb || $a ;
}
############################################################################
# Prints the dependency rule for an f90 main program
# Returns the no dependency stuffs, then the ordered modules then the mains:
#  (@lnouse, @orders, @lmains90) ;
sub sortallf90 {
  my(@olds) = (@_) ;
  foreach my $msub (@olds) {
    if ($isamain{$msub}) {
      push(@lmains90, $msub) ;
    }
  }
  local(%writens) = () ;
  local($outline) = " " ;
  my $space = "  " ;
  foreach my $msub (@lmains90) {
    # last if ($msub eq "fdtd") ;
    print STDOUT "Analysis for $msub\n" if $debug ;
    &setUses($msub, 0, "$space") ;
  }
  foreach my $msub (@olds) {
    unless ($isamain{$msub}) {
      unless (exists $writens{$msub}) {
	unless (exists $revmods{$msub}) {
	  # Not a module, a global stuff
	  push(@lnouse90, $msub) ;
	}
      }
    }
  }
}
############################################################################
# Prints the dependency rule for an f90 main program
sub printMain {
  my $msub = shift ;
  my $typ = shift ;
  my $mode = shift ;
  my $term = "" ;
  my $code = $msub ;
  # Remove directory:
  $code = $2 if ($code =~ /(.*)\/(.*)/) ;
  if ($mode == 1) {
    print LSTDEP "\n\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o: " ;
  } else {
    # Making the library:
    print LSTDEP "\n\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/lib$msub.so: " ;
    $term = ".\$(DYLIB)" ;
  }
  # Remove the dependencies from this list ?
  if (!$reuse || !$opt_f) {
    if (%{$calls{$msub}} || %{$uses{$msub}}) {
      print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.F90 " if $opt_V ;
      print LSTDEP " \$(OBJF) \$(OBJf) \$(OBJs90) \$(OBJC) \$(OBJCPP) " ;
    }
  }
  if ($opt_f) {
    # Individual dependency is written: force this for the C file...
    if ($opt_l) {
      $lib1 = '$(LIB)(' ;
      $lib2 = ')' ;
    } else {
      $lib1 = "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/" ;
      $lib2 = '' ;
    }
    local(%writens) = () ;
    local($outline) = " " ;
    local(@orders) = () ;
    if ($reuse) {
      print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.F90 " if $opt_V ;
    }
    #print LSTDEP " AA($msub $typ) " ;
    &printrUses($msub, $typ) ;
    #print LSTDEP " BB " ;
    # OK, write the OBJF list but with some files being removed:
    if ($reuse) {
      foreach (@listef,@listeF,@listec,@listepp,@lnouse90) {
	next if ($isamain{$_}) ;
	next if ($writens{$_}) ;
	next if (%{$uses{$_}}) ;
	next unless $_ ;
	if ($opt_l) {
	  print LSTDEP "\$(LIB)($_.o) " ;
	} else {
	  print LSTDEP "$lib1$_.o$lib2 " ;
	}
      }
    }
  } elsif ($opt_d) {
    # Recursive analysis shall be performed, independent file are stored in OBJs90
    # 24/03/2010: adding all the f90 files in doubt...
    print LSTDEP " \$(OBJf90) \$(OBJF90) " ;
    &printUses($msub, 1) ;
  } else {
    # Simple compilation must be performed, adding all the f90 file dependencies
    print LSTDEP " \$(OBJf90) \$(OBJF90) " ;
    &printUses($msub, 0) ;
  }
  if ($mode == 1) {
    print LSTDEP " $msub.$typ @{$dependencies{$msub}} \n" ;
  } else {
    print LSTDEP "\n" ;
  }
  # Compilation rule:
  print LSTDEP "\t\@echo \"Compiled $msub library\"\n" ;
  # Release update:
  print LSTDEP "\tperl $mkversion -c '\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)' -i '\$(HELP)' -h '\$(SECTION)\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)' -v '\$(VERSION)' -d '\$(DATE)' -t '\$(TAG)'\n" ;
  if (0) {
    print LSTDEP "\t\$(MAKE) theversion TYPEREEL=\$(TYPEREEL) OUT=\$(OUT)\n" ;
  } else {
    if (exists $ldcmpl{"F90"}) {
      print LSTDEP "\t\$($ldcmpl{'F90'}) $ldopts{'F90'} -c \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.F90 -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o\n" ;
    } else {
      print LSTDEP "\t\$(F90) -DTYPEREEL=\$(TYPEREEL) \$(FOPTS) \$(F90FLAGS) \$(OPT) \$(INCLUDES) -c \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.F90 -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o\n" ;
    }
  }
  # Main program compilation:
  if ($mode == 1) {
    if ($typ eq 'c' || $typ eq 'cpp') {
      # print STDOUT "$typ has ldcmpl=$ldcmpl{$typ}, opts=$ldopts{$typ}\n" ;
      if (exists $ldcmpl{$typ}) {
	print LSTDEP "\t\$($ldcmpl{$typ}) $ldopts{$typ} -c $msub.$typ -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \n" ;
      } else {
	print LSTDEP "\t\$(CC) \$(CFLAGS) \$(INCLUDES) -c $msub.$typ -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o\n" ;
      }
    } elsif (exists $ldcmpl{$typ}) {
      print LSTDEP "\t\$($ldcmpl{$typ}) $ldopts{$typ} -c $msub.$typ -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \n" ;
    } else {
      print LSTDEP "\t\$(F90) -DTYPEREEL=\$(TYPEREEL) \$(FOPTS) \$(F90FLAGS) \$(OPT) \$(INCLUDES) -c $msub.$typ -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o\n" ;
    }
  }
  my $target ;
  my $shared = "" ;
  if ($mode == 1) {
    $target = "a." ;
  } else {
    $target = "lib" ;
    $shared = " \$(SHARED)" ;
  }
  # Linking rule:
  # Executable (or shared library) creation rule:
  print LSTDEP "\$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term: \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o\n" ;
  if ($typ eq 'c' || $typ eq 'cpp') {
    # When C calls fortran, set the compilation rule as followed
    my $lds = 'LD' . uc($typ) ;
    if ($typ eq 'c') {
	$lds = 'LDC' ; # 'CC'
    }
    # print STDOUT "lds=$lds\n" ;
    if (%{$calls{$msub}}) {
      if ($opt_l) {
	print LSTDEP "\t\$($lds)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO)" ;
	if ($mode == 1) {
	  print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o \$(LIB) \n" ;
	} else {
	  print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o \$(LIB) \n" ;
	}
      } else {
	if ($opt_f && $reuse) {
	  print LSTDEP "\t\$($lds)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO)" ;
	  if ($mode == 1) {
	    print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o " ;
	  } else {
	    print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o " ;
	  }
	  $lib1 = "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/" ;
	  $lib2 = '' ;
	  local(%writens) = () ;
	  local($outline) = " " ;
	  local(@orders) = () ;
	  &printrUses($msub, $typ) ;
	  foreach (@listef,@listeF,@listec,@listepp,@lnouse90) {
	    next if ($isamain{$_}) ;
	    next if ($writens{$_}) ;
	    next if (%{$uses{$_}}) ;
	    if ($opt_l) {
	      print LSTDEP "\$(LIB)($_.o) " ;
	    } else {
	      print LSTDEP "$lib1$_.o$lib2 " ;
	    }
	  }
	  #print LSTDEP " \$(myOBJF) \$(myOBJs90) \$(myOBJC) \$(myOBJCPP)" ;
	  print LSTDEP " \n" ;
	  # What the hell is that for ?
	  if ($opt_l) {
	    $lib1 = '$(LIB)(' ;
	    $lib2 = ')' ;
	  } else {
	    $lib1 = "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/" ;
	    $lib2 = '' ;
	  }
	} else {
	  print LSTDEP "\t\$($lds)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO)" ;
	  if ($mode == 1) {
	    print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o \$(OBJC) \$(OBJF) \$(OBJf) \$(OBJf90) \$(OBJF90) \n" ;
	  } else {
	    print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o \$(OBJC) \$(OBJF) \$(OBJf) \$(OBJf90) \$(OBJF90) \n" ;
	  }
	}
      }
    } else {
      # No longer add  -nofor_main for Sun...
      if ($mode == 1) {
	print LSTDEP "\t\$($lds)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO) \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \n" ;
      } else {
	print LSTDEP "\t\$($lds)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO) \n" ;
      }
    }
  } elsif (%{$uses{$msub}}) {
    if ($opt_f) {
      print LSTDEP "\t\$(LD)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term " ;
      if ($mode == 1) {
	print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o" ;
      } else {
	print LSTDEP " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o" ;
      }
      if (!$reuse) {
	print LSTDEP " \$(myOBJF) \$(myOBJf) \$(myOBJs90) \$(myOBJC) \$(myOBJCPP)" ;
      }
      $lib1 = "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/" ;
      $lib2 = '' ;
      local(%writens) = () ;
      local($outline) = " " ;
      local(@orders) = () ;
      &printrUses($msub, $typ) ;
      if ($reuse) {
	foreach (@listef,@listeF,@listec,@listepp,@lnouse90) {
	  next if ($isamain{$_}) ;
	  next if ($writens{$_}) ;
	  next if (%{$uses{$_}}) ;
	  if ($opt_l) {
	    print LSTDEP "\$(LIB)($_.o) " ;
	  } else {
	    print LSTDEP "$lib1$_.o$lib2 " ;
	  }
	}
      }
      print LSTDEP "  \$(BIBLIO)\n" ;
      if ($opt_l) {
	$lib1 = '$(LIB)(' ;
	$lib2 = ')' ;
      } else {
	$lib1 = "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/" ;
	$lib2 = '' ;
      }
    } else {
      my $wl = "" ;
      if ($mode == 1) {
	$wl = " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o" ;
      } else {
	$wl = " \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o" ;
      }
      if ($opt_l) {
	print LSTDEP "\t\$(LD)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO)$wl \$(LIB) \n" ;
      } else {
	print LSTDEP "\t\$(LD)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO)$wl \$(myOBJF) \$(myOBJf) \$(myOBJs90) \$(myOBJC) \$(myOBJf90) \$(myOBJF90) \n" ;
      }
    }
  } else {
    if ($mode == 1) {
      print LSTDEP "\t\$(LD)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO) \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$msub.o \n" ;
    } else {
      print LSTDEP "\t\$(LD)$shared \$(OPL) -o \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term \$(BIBLIO) \n" ;
    }
  }
  print LSTDEP "\t\@echo \"OK, Generated \$(BIN)$target$code\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)$term\"\n" ;
}
############################################################################
# Tries to guess how to compile main program:
sub getRules {
  my $make = shift ;
  my $section = "" ;
  my (%compilers) = ( 'F77' => 'f', 'F90' => 'F90', 'f90' => 'f90', 'CC' => 'c', 'CPP' => 'cpp', 'cpp' => 'cpp',) ;
  if (open(MK, "<$make")) {
    while(<MK>) {
      if ($section) {
	if (/\s+\$\((\w+)\)\s+(.*)/) {
	  if (exists $compilers{$1}) {
	    # Another rule may be built:
	    # next if exists $ldopts{$section} ;
	    print STDOUT "\tNew rule($.) $section compiler $1 is $2\n" if $debug ;
	    $ldcmpl{$section} = $1 ;
	    $ldopts{$section} = $2 ;
	    if ($ldopts{$section} =~ /(.*)\s+\-c\s+(.*?)\s+\-o\s+.*\/(.*?)\s+(.*)$/) {
	      $ldopts{$section} = $1 . $4 ;
	    } elsif ($ldopts{$section} =~ /(.*)\s+\-c\s+(.*?)\s+\-o\s+.*\/(.*)$/) {
	      $ldopts{$section} = $1 ;
	    } elsif ($ldopts{$section} =~ /(.*)\s+\-c\s+(.*?)\s+\-o\s+(\S+)\s*$/) {
	      # New definition style:
	      $ldopts{$section} = $1 ;
	    } elsif ($ldopts{$section} =~ /(.*)\s+\-c\s+(.*?)\s+\-o\s+\$\@\s+\$\<$/) {
	      $ldopts{$section} = "$1 $2" ;
	    }
	    print STDOUT "\t\tNew rule($.) $section compiler $ldopts{$section}\n" if $debug ;
	  }
	}
      }
      if (/^\.(\w+)\.o\s*:\s*$/) {
	# Old section style:
	print STDOUT "New section $1\n" if $debug ;
	$section = $1 ;
      } elsif (/^.*\%\.o\s*:\s*\%\.(\w+)$/) {
	# New section style:
	print STDOUT "New section $1\n" if $debug ;
	$section = $1 ;
      } elsif (/^\.\w+\.\w+\s*:\s*$/) {
	$section = "" ;
      } elsif (/^\w+\s*:\s*$/) {
	$section = "" ;
      }
    }
    close(MK) ;
  }
}
############################################################################
# Prints the dependency rule for an f90 subroutine
sub printRoutine {
  my $msub = shift ;
  my $typ = shift ;
  if ($opt_l) {
    if ($msub =~ /mpi/i) {
      print LSTDEP "\$(LIB)($msub.o):" ;
    } else {
      print LSTDEP "\$(LIB)($msub.o):" ;
    }
  } else {
    if ($msub =~ /mpi/i) {
      print LSTDEP "$lib1$msub.o$lib2:" ;
    } else {
      print LSTDEP "$lib1$msub.o$lib2:" ;
    }
  }
  if ($opt_d) {
    &printUses($msub, 1) ;
  } else {
    &printUses($msub, 0) ;
  }
  print LSTDEP " $msub.$typ @{$dependencies{$msub}}\n" ;
}
############################################################################
# Prints dependencies for an f90 routine or program
sub printUses {
  my $msub = shift ;
  my $depend = shift ;
  if ($depend) {
    (%writens) = () ;
    foreach $module (sort bylevel keys %{$uses{$msub}}) {
      $dep = $modules{$module} ;
      print "File for module='$module' is $dep.\n" if ($depend == 2);
      unless (exists $writens{$dep}) {
	$writens{$dep} = 1 ;
	next if ($dep eq "Version") ;
	if ($dep =~ /mpi/i) {
	  print LSTDEP ' ' . $lib1 . $dep . '.o' . $lib2;
	} else {
	  print LSTDEP ' ' . $lib1 . $dep . '.o' . $lib2;
	}
      }
    }
  }
}
############################################################################
# Recursively builds the stuff:
sub setCallUses {
  my $msub = shift ;
  my $level = shift ;
  my $space = shift ;
  $maxlev = max($maxlev, $level) ;
  print "$space Analyzing $msub dependencies:\n" if $debug;
  my(@a) = () ;
  my(%b) = () ;
  my(@c) = () ;
  foreach my $f (keys %{$calls{$msub}}) {
    print STDOUT "\tf=$f\n" if $debug ;
    my $dep = $f ;
    $dep = $1 if ($dep =~ /(.*)\./) ;
    next if ($f =~ /\.\w*h/i) ;
    # Strange, has no effect:
    next if $isamain{$f} ;
    next if $isamain{$dep} ;
    # Refuse the modules as being called directly without any interface.
    next if exists $writens{$dep} ;
    next if exists $revmods{$dep} ;
    next unless (-s $f) ;
    unless (exists $writens{$dep}) {
      print "$space File for calls of $msub include $dep.\n" if $debug ;
      unless (exists $b{$dep}) {
	if (%{$uses{$dep}} || %{$calls{$dep}}) {
	  push(@a,$dep) ;
	} else {
	  push(@c,$dep) ;
	}
	$b{$dep} = 1 ;
      }
    }
  }
  foreach my $module (keys %{$uses{$msub}}) {
    my $dep = $modules{$module} ;
    #next if $isamain{$dep} ;
    unless (exists $writens{$dep}) {
      print "$space File for module='$module'($msub) is $dep.\n" if $debug;
      unless (exists $b{$dep}) {
	if (%{$uses{$dep}} || %{$calls{$dep}}) {
	  push(@a,$dep) ;
	} else {
	  push(@c,$dep) ;
	}
	$b{$dep} = 1 ;
      }
    }
  }
  foreach my $dep (@a,@c) {
    next if exists $writens{$dep} ;
    # $outline .= ' $(LIB)(' . $dep . '.o) ' . "@{$dependencies{$dep}}" ;
    $outline .= ' ' . $lib1 . $dep . '.o' . $lib2 . ' ' ;
    $writens{$dep} = ' ' . $lib1 . $dep . '.o' . $lib2 . ' ' ;
    if (%{$uses{$dep}} || %{$calls{$dep}}) {
      print STDOUT "$space further use analysis for $dep from $msub\n" if $debug ;
      &setCallUses($dep, ($level+1), "$space  ")  ;
      push(@orders, $dep) ;
    } else {
      # Already performed stuff:
      print STDOUT "$space Storing use $dep\n" if $debug ;
      push(@orders, $dep) ;
    }
  }
  print "$space Performed $msub dependencies:\n" if $debug;
}
############################################################################
# Recursively builds the stuff:
sub setCalls {
  my $msub = shift ;
  my $level = shift ;
  my $space = shift ;
  # $debug = 1 if ($msub eq 'main') ;
  $maxlev = max($maxlev, $level) ;
  print "$space Analyzing $msub dependencies:$calls{$msub}\n" if $debug;
  foreach my $f (keys %{$calls{$msub}}) {
    print STDOUT "\tf=$f\n" if $debug ;
    my $dep = $f ;
    $dep = $1 if ($dep =~ /(.*)\./) ;
    next if ($f =~ /\.\w*h/i) ;
    # Strange, has no effect:
    next if $isamain{$f} ;
    next if $isamain{$dep} ;
    # Refuse the modules as being called directly without any interface.
    next if exists $writens{$dep} ;
    next if exists $revmods{$dep} ;
    next unless (-s $f) ;
    unless (exists $writens{$dep}) {
      print "$space File for calls of $msub include $dep.\n" if $debug ;
      # $outline .= ' $(LIB)(' . $dep . '.o) ' . "@{$dependencies{$dep}}" ;
      $outline .= ' ' . $lib1 . $dep . '.o' . $lib2 . ' ' ;
      $writens{$dep} = ' ' . $lib1 . $dep . '.o' . $lib2 . ' ' ;
      if (%{$calls{$dep}}) {
	print STDOUT "$space further call analysis for $dep from $msub\n" if $debug ;
	&setCalls($dep, ($level+1), "$space  ")  ;
	push(@orders, $dep) ;
      } else {
	# Already performed stuff:
	print STDOUT "$space Storing call $dep\n" if $debug ;
	push(@orders, $dep) ;
      }
    }
  }
  #$debug = 0 ;
}
############################################################################
# Recursively builds the stuff:
sub setUses {
  my $msub = shift ;
  my $level = shift ;
  my $space = shift ;
  $maxlev = max($maxlev, $level) ;
  print "$space Analyzing $msub dependencies:\n" if $debug;
  foreach my $module (keys %{$uses{$msub}}) {
    my $dep = $modules{$module} ;
    unless (exists $writens{$dep}) {
      print "$space File for module='$module'($msub) is $dep.\n" if $debug;
      # $outline .= ' $(LIB)(' . $dep . '.o) ' . "@{$dependencies{$dep}}" ;
      $outline .= ' ' . $lib1 . $dep . '.o' . $lib2 . ' ' ;
      $writens{$dep} = ' ' . $lib1 . $dep . '.o' . $lib2 . ' ' ;
      if (%{$uses{$dep}}) {
	print STDOUT "$space further use analysis for $dep from $msub\n" if $debug ;
	&setUses($dep, ($level+1), "$space  ")  ;
	push(@orders, $dep) ;
      } else {
	# Already performed stuff:
	print STDOUT "$space Storing use $dep\n" if $debug ;
	push(@orders, $dep) ;
      }
    }
  }
  print "$space Performed $msub dependencies:\n" if $debug;
}
############################################################################
# Prints recursive dependencies for an f90 routine or program:
sub printrUses {
  my $msub = shift ;
  my $typ = shift ;
  my $space = "  " ;
  # Case of C call... n,ot very clean manner:
  if (1 == 0 && exists $calls{$msub}) {
    foreach my $lsub (keys %{$calls{$msub}}) {
      print STDOUT "Parsing calls for $lsub from $msub\n" if $debug ;
      &setUses($lsub, 1, "$space") ;
    }
    foreach my $dep (@orders) {
      next if ($dep eq "Version") ;
      print LSTDEP " $writens{$dep}" ;
      print STDOUT " $writens{$dep}" if $debug ;
    }
    # Not very safe stuff, based on the fact the called routine MUST
    # be defined in a file with same name !!!!
    # foreach my $lsub (keys %{$calls{$msub}}) {
    #   print LSTDEP " $lib1$lsub.o$lib2" ;
    # }
  } else {
    # Old technique:
    if (1 == 0) {
      $maxlev = 0 ;
      &setUses($msub, 0, "$space") ;
      # Only now add the call dependencies:
      $imaxlev = $maxlev ;
      &setCalls($msub, $imaxlev, "$space") ;
      # Make a first sort:
      (@porders) = (sort constlevel @orders) ;
      #die "orders=@orders\n" ;
      # Add the call dependencies for all routines:
      foreach my $dep (@porders) {
	&setCalls($dep, $imaxlev+1, "$space$space") ;
      }
      # Now add the uses from all the calls again, this may not be in perfect order...
      foreach my $dep (@porders) {
	&setUses($dep, $imaxlev, "$space") ;
      }
      push(@orders,@porders) ;
    } else {
      &setCallUses($msub, 0, "$space") ;
    }
    print STDOUT "Final dep $msub: @orders\n" if $debug ;
    for (my $i=0; $i < $#orders; $i++) {
      $rang{$orders[$i]} = $i ;
    }
    if (0) {
      print STDOUT "Working for $msub\n" ;
      my(@out) = (sort constlevel @orders) ;
      for (my $i=0; $i < $#orders; $i++) {
	if ($orders[$i] ne $out[$i]) {
	  print STDOUT "Index $i, file $orders[$i] replaced by $out[$i]\n" if $debug;
	}
      }
      print STDOUT "Worked for $msub\n" ;
    }
    # Here another last chance to geta good sort, or a chance to get a bad one from a proper one !
    # print LSTDEP "orders=@orders\n" ;
    foreach my $dep (sort constlevel @orders) {
    # foreach my $dep (@orders) {
      next if ($dep eq "Version") ;
      next unless $dep ;
      print LSTDEP " $writens{$dep}" ;
      print STDOUT " $writens{$dep}" if $debug;
    }
  }
  print STDOUT "\nRecursive dep is $outline\n" if $debug ;
}
############################################################################
# A spcial max format
sub max {
  my $max = shift(@_);
  foreach my $foo (@_) {
    my $b1 = $max ;
    my $a1 = $foo ;
    my $b2 = 0 ;
    my $a2 = 0 ;
    if ($max =~ /V(\d+)\_(\d+)/) {
      $b1 = $1 ; $b2 = $2 ;
    }
    if ($foo =~ /V(\d+)\_(\d+)/) {
      $a1 = $1 ; $a2 = $2 ;
    }
    if ($a1 > $b1) {
      $max = $foo ;
    } elsif ($a1 == $b1) {
      if ($a2 > $b2) {
	$max = $foo ;
      }
    }
  }
  return $max;
}
############################################################################
# Takes the most recent date of change
sub mergeDates {
  my $a = shift ;
  my $b = shift ;
  my $r = "" ;
  my (@months) = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec") ;
  # Date format is Tue Jun 20 12:14:06 2006
  #                0   1   2   3        4
  my (@as) = split(/\s+/, $a) ;
  my (@bs) = split(/\s+/, $b) ;
  # Year compare:
  if ($bs[4] > $as[4]) {
    $r = $b ;
  } elsif ($bs[4] < $as[4]) {
    $r = $a ;
  } else {
    for (my $i=11; $i>=0; $i--) {
      if ($bs[1] eq $months[$i] || $as[1] eq $months[$i]) {
	if ($bs[1] eq $months[$i] && $as[1] eq $months[$i]) {
	  # Both months equal, look for day:
	  if ($bs[2] > $as[2]) {
	    $r = $b ;
	  } elsif ($bs[2] < $as[2]) {
	    $r = $a ;
	  } else {
	    # Look for time, hours:
	    my (@ad) = split(/:/, $as[3]) ;
	    my (@bd) = split(/:/, $bs[3]) ;
	    if ($bd[0] > $ad[0]) {
	      $r = $b ;
	    } elsif ($bd[0] < $ad[0]) {
	      $r = $a ;
	    } else {
	      # Look for minutes:
	      if ($bd[1] > $ad[1]) {
		$r = $b ;
	      } elsif ($bd[1] < $ad[1]) {
		$r = $a ;
	      } else {
		# Look for seconds:
		if ($bd[2] > $ad[2]) {
		  $r = $b ;
		} elsif ($bd[2] < $ad[2]) {
		  $r = $a ;
		} else {
		  # Times are absolutely equal
		  $r = $b ;
		}
	      }
	    }
	  }
	} elsif ($bs[1] eq $months[$i]) {
	  $r = $b ;
	} elsif ($as[1] eq $months[$i]) {
	  $r = $a ;
	}
	last ;
      }
    }
  }
#  print STDOUT "cmp $a with $b is $r\n";
  return($r) ;
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
    $entryfile = "../../sophie/src/.svn/entries" ;
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
	  print STDOUT "SVN version=$version file $entryfile\n" if 1 || $debug ;
	}
	# Let us get the current branch or tag:
	if ($. == 5) {
	  # indicates repository location:
	  if (/\b(branch|tag)e*s*\/([\w\.]+)/) {
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
	    #die "data=$date" ;
	    print STDOUT "date=$date\n" if $debug ;
	  }
	}
      }
      close(SVN) ;
      print STDOUT "tag=$tag, version=$version\n" if $debug ;
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
	  print STDOUT "Using SVN dir '$svn'\n" if 1 || $debug ;
	  $remsh = "svn diff -r $svn" ;
	  if (open(SYS, "$remsh 2>&1 |")) {
	    while(<SYS>) {
	      print ;
	      if (/^svn:/) {
		$version .= "unknown" ;
	      } else {
		$version .= "modif" ;
	      }
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
# Recursive add files to @liste:
sub pushFiles {
  my $dir = shift ;
  my $target = $dir ;
  $target =~ s/^\.\/*// ;
  $target .= '/' if $target ;
  if (opendir(MONDIRECTOIRE,"$dir")) {
    my (@dirs) = grep { $_ if (-d "$dir/$_" && $_ ne '.' && $_ ne '..' && !($_ =~ /\.(svn|cvs)/i)) } readdir(MONDIRECTOIRE);
    # Open/Close required:
    closedir MONDIRECTOIRE;
    opendir(MONDIRECTOIRE,"$dir") ;
    my (@files) = map { "$target$_" } grep { $_ if (-e "$dir/$_") } grep /\.(f90|f|c|cc|cu|cpp|c\+\+)$/i, readdir(MONDIRECTOIRE);
    # print STDOUT "List($dir [@dirs]) is: @files\n" if (@files && $opt_v) ;
    closedir MONDIRECTOIRE;
    push(@liste, @files) ;
    if ($opt_L) {
      foreach my $d (@dirs) {
	#print STDOUT "Adding $dir/$d\n" ;
	print STDOUT "Testing $dir/$d\n" if $debug ;
	next unless $okdirs{"$dir/$d"} ;
	print STDOUT "Accepting $dir/$d\n" if $debug ;
	&pushFiles("$dir/$d") ;
      }
    } else {
      foreach my $d (@dirs) {
	#print STDOUT "Adding $dir/$d\n" ;
	&pushFiles("$dir/$d") ;
      }
    }
  }
}
############################################################################
sub getList {
  if ($opt_s) {
    (@liste) = ("$opt_s") ;
  } elsif ($opt_R) {
    # Recursively find the fortran files:
    my (@dirs) = grep { $_ if (-d $_ && $_ ne '.' && $_ ne '..' && !($_ =~ /\.(svn|cvs)/i)) } readdir(DIR) ;
    if ($opt_L) {
      $opt_L =~ s/\s+/:/g;
      my @a = split(/:/, $opt_L) ;
      (%okdirs) = () ;
      foreach my $b (@a) {
	# We ought to accept all the previous dirs:
	$okdirs{$b} = 1 ;
	$okdirs{"./$b"} = 1 ;
	# Accept all the father directories:
	my(@c) = split(/\//, $b) ;
	my $f = $c[0] ;
	for (my $i=1; $i <= $#c; $i++) {
	  print STDOUT "Validating directory $f\n" if $debug ;
	  $okdirs{$f} = 1 ;
	  $okdirs{"./$f"} = 1 ;
	  $f .= '/' . $c[$i] ; 
	}
      }
      print STDOUT "Used directory list is @a\n" if $debug ;
      foreach my $d (@dirs) {
	next unless $okdirs{$d} ;
	&pushFiles("./$d") ;
      }
    } else {
      foreach my $d (@dirs) {
	&pushFiles("./$d") ;
      }
    }
    print STDOUT "Recursively added from @dirs following files\n\t@liste\n" if $debug ;
  } else {
    (@liste) = grep { !/#\./ && -f $_ } readdir(DIR) ;
  }
}
############################################################################
sub printAll {
  open(LSTOBJ,">$dir/mkobj.lst") or die "Write permission $dir/mkobj.lst denied" ;
  open(LSTDEP,">$dir/mkdep.lst") or die "Write permission $dir/mkdep.lst denied" ;
  print LSTOBJ "VERSION = $version\n" ;
  print LSTOBJ "DATE = $date\n" ;
  print LSTOBJ "TAG = $tag\n" ;
  print LSTOBJ "MD5 ?= $md5\n" if $md5 ;
  print LSTOBJ "SVN = $svn\n" if $svn ;
  if (@listecpp) {
    print LSTOBJ "OBJCPP = " ;
    foreach (@listecpp) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTOBJ "\$(LIB)($_.o) " ;
      } else {
	print LSTOBJ "$lib1$_.o$lib2 " ;
      }
    }
    print LSTOBJ "\n" ;
    print LSTOBJ "myOBJCPP = " ;
    foreach (@listecpp) {
      next if ($isamain{$_}) ;
      print LSTOBJ "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$_.o " ;
    }
    print LSTOBJ "\n" ;
    foreach (@listecpp) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTDEP "\$(LIB)($_.o): $_.cpp @{$dependencies{$_}} \n" ;
      } else {
	print LSTDEP "$lib1$_.o$lib2: $_.cpp @{$dependencies{$_}} \n" ;
      }
    }
  }
  if (@listec) {
    print LSTOBJ "OBJC = " ;
    foreach (@listec) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTOBJ "\$(LIB)($_.o) " ;
      } else {
	print LSTOBJ "$lib1$_.o$lib2 " ;
      }
    }
    print LSTOBJ "\n" ;
    print LSTOBJ "myOBJC = " ;
    foreach (@listec) {
      next if ($isamain{$_}) ;
      print LSTOBJ "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$_.o " ;
    }
    print LSTOBJ "\n" ;
    foreach (@listec) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTDEP "\$(LIB)($_.o): $_.c @{$dependencies{$_}} \n" ;
      } else {
	print LSTDEP "$lib1$_.o$lib2: $_.c @{$dependencies{$_}} \n" ;
      }
    }
  }
    else {
      print LSTOBJ "OBJC = \n" ;
    }
  if (@listef) {
    print LSTOBJ "OBJf = " ;
    foreach (@listef) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTOBJ "\$(LIB)($_.o) " ;
      } else {
	print LSTOBJ "$lib1$_.o$lib2 " ;
      }
    }
    print LSTOBJ "\n" ;
    print LSTOBJ "myOBJf = " ;
    foreach (@listef) {
      next if ($isamain{$_}) ;
      print LSTOBJ "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$_.o " ;
    }
    print LSTOBJ "\n" ;
    foreach (@listef) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTDEP "\$(LIB)($_.o): $_.f @{$dependencies{$_}} \n" ;
      } else {
	print LSTDEP "$lib1$_.o$lib2: $_.f @{$dependencies{$_}} \n" ;
      }
    }
  }
  if (@listeF) {
    print LSTOBJ "OBJF = " ;
    foreach (@listeF) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTOBJ "\$(LIB)($_.o) " ;
      } else {
	print LSTOBJ "$lib1$_.o$lib2 " ;
      }
    }
    print LSTOBJ "\n" ;
    print LSTOBJ "myOBJF = " ;
    foreach (@listeF) {
      next if ($isamain{$_}) ;
      print LSTOBJ "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$_.o " ;
    }
    print LSTOBJ "\n" ;
    foreach (@listeF) {
      next if ($isamain{$_}) ;
      if ($opt_l) {
	print LSTDEP "\$(LIB)($_.o): $_.F @{$dependencies{$_}} \n" ;
      } else {
	print LSTDEP "$lib1$_.o$lib2: $_.F @{$dependencies{$_}} \n" ;
      }
    }
  }
  if (@listef90) {
    print LSTOBJ "\nOBJf90 = " ;
    # From sortallf90 and 
    foreach (@lnouse90,@orders,@othf90) {
      next if ($types{$_} ne 'f90') ;
      if ($isamain{$_}) {
	print STDOUT "Main $_\n" if ($debug) ;
      } else {
	if ($opt_l) {
	  print LSTOBJ "\$(LIB)($_.o) " ;
	} else {
	  print LSTOBJ "$lib1$_.o$lib2 " ;
	}
      }
    }
    print LSTOBJ "\n" ;
    print LSTOBJ "myOBJf90 = " ;
    foreach (@lnouse90,@orders,@othf90) {
      next if ($types{$_} ne 'f90') ;
      unless ($isamain{$_}) {
	print LSTOBJ "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$_.o " ;
      }
    }
    print LSTOBJ "\n" ;
  }
  if (@listeF90) {
    print LSTOBJ "\nOBJF90 = " ;
    foreach (@lnouse90,@orders,@othf90) {
      next if ($types{$_} ne 'F90') ;
      next if ($_ eq "Version") ;
      if ($isamain{$_}) {
	print STDOUT "Main $_\n" if ($debug) ;
      } else {
	if ($opt_l) {
	  print LSTOBJ "\$(LIB)($_.o) " ;
	} else {
	  print LSTOBJ "$lib1$_.o$lib2 " ;
	}
      }
    }
    print LSTOBJ "\n" ;
    print LSTOBJ "myOBJF90 = " ;
    foreach (@lnouse90,@orders,@othf90) {
      next if ($types{$_} ne 'F90') ;
      next if ($_ eq "Version") ;
      unless ($isamain{$_}) {
	print LSTOBJ "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$_.o " ;
      }
    }
    print LSTOBJ "\n" ;
  }
  if (@listeall90) {
    print LSTOBJ "\nOBJs90 = " ;
    foreach (@lnouse90) {
      next if ($isamain{$_}) ;
      next if (%{$uses{$_}}) ;
      if ($opt_l) {
	print LSTOBJ "\$(LIB)($_.o) " ;
      } else {
	print LSTOBJ "$lib1$_.o$lib2 " ;
      }
    }
    print LSTOBJ "\n" ;
    print LSTOBJ "myOBJs90 = " ;
    foreach (@lnouse90) {
      next if ($isamain{$_}) ;
      next if (%{$uses{$_}}) ;
      next if ($eliminate{$_}) ;
      if ($opt_l) {
	print LSTOBJ "\$(LIB)($_.o) " ;
      } else {
	print LSTOBJ "$lib1$_.o$lib2 " ;
      }
      # print LSTOBJ "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/$_.o " ;
    }
    print LSTOBJ "\n" ;
  }
  if (@listeall90) {
    print LSTDEP "\n" ;
    # First prints routine dependencies,
    foreach (@lnouse90, @orders,@othf90) {
      next if ($_ eq "Version") ;
      unless ($isamain{$_}) {
	&printRoutine($_, $types{$_}) ;
      }
    }
    # then main dependencies
    foreach (@lmains90) {
      if ($isamain{$_}) {
	&printMain($_, $types{$_}, 1) ;
	&printMain($_, $types{$_}, 2) ;
      }
    }
  }
  if (@listec || @listecpp) {
    # main dependencies for C and C++ programs:
    foreach (@listec,@listecpp) {
      next unless $_ ;
      if ($isamain{$_}) {
	&printMain($_, $types{$_}, 1) ;
	&printMain($_, $types{$_}, 2) ;
      }
    }
  }
  # Prints theversion rule in case of -d:
  # if ($opt_d && !$opt_f) {
  #   print LSTDEP "theversion:\n" ;
  # Prints the version rule unless -V:
  if ($opt_V) {
    print LSTDEP "\n" ;
    print LSTDEP "\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.F90: \n" ;
    print LSTDEP "\t\@echo \"Depend VERSION\"\n" ;
    print LSTDEP "\tperl $mkversion -c '\$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)' -i '\$(HELP)' -h '\$(SECTION)\$(OUT)\$(TYPEREEL)_\$(HOSTARCH)' -v '\$(VERSION)' -d '\$(DATE)' -t '\$(TAG)'\n" ;
    if (exists $ldcmpl{"F90"}) {
      print LSTDEP "\t\$($ldcmpl{'F90'}) $ldopts{'F90'} -c \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.F90 -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o\n" ;
    } else {
      print LSTDEP "\t\$(F90) -DTYPEREEL=\$(TYPEREEL) \$(FOPTS) \$(F90FLAGS) \$(OPT) \$(INCLUDES) -c \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.F90 -o \$(TMPLIB)/objets\$(OUT)\$(TYPEREEL)/\$(HOSTARCH)/Version.o\n" ;
    } 
  }
  # Closes the target dependencies files
  close(LSTOBJ) ;
  close(LSTDEP) ;
}
############################################################################
# Main loop on all directories
sub mainRout {
  foreach $dir (@dirs) {
    unless (opendir(DIR,".")) {
      print STDERR "Read permission denied on directory $dir\n" ;
      next ;
    }
    # Gets the actual disk location for includes from current local directory
    (%getrel) = () ;
    if (@rincludes) {
      foreach my $i (@rincludes) {
	# Gets a location that exists on disk from current directory.
	$getrel{$i} = &getAbs($dir,$i) ;
      }
    }
    # Gets the directory's file list
    &getList() ;
    closedir(DIR) ;
    # Gets the list of possible makefiles:
    (@listemkf) = sort grep /makefile/i, @liste ;
    # Looking for makefile rules:
    foreach my $file (@listemkf) {
      next if ($file =~ /[\%\~]$/) ;
      print STDOUT "Read processing rules for $file\n" if $debug ;
      &getRules("$dir/$file") ;
    }
    ($version,$date,$tag) = &getVersion() ;
    # print "tag=$tag\n" ;
    (%isamain) = () ;
    (%dones) = () ;
    (%types) = () ;

    (@listecpp) = sort grep /\.cpp$/, @liste ;
    foreach my $file (@listecpp) {
      $file =~ s/\.cpp$// ;
      $types{$file} = 'cpp' ;
      &getDependc("$dir/$file.cpp", $file, $dir) ;
    }
    (@listec) = sort grep /\.c$/, @liste ;
    foreach my $file (@listec) {
      $file =~ s/\.c$// ;
      $types{$file} = 'c' ;
      &getDependc("$dir/$file.c", $file, $dir) ;
    }
    (@listef) = grep /\.f$/, @liste ;
    foreach my $file (@listef) {
      if ($file =~ /^(.*)\.(f)$/) {
	$types{$1} = $2 ;
      }
      $file =~ s/\.f$// ;
      &getDependf("$dir/$file.f", $file, $dir) ;
    }
    (@listeF) = grep /\.F$/, @liste ;
    foreach my $file (@listeF) {
      if ($file =~ /^(.*)\.(F)$/) {
	$types{$1} = $2 ;
      }
      $file =~ s/\.F$// ;
      &getDependf("$dir/$file.F", $file, $dir) ;
    }
    (@listef90) = sort grep /\.f90$/, @liste ;
    (@listeF90) = sort grep /\.F90$/, @liste ;
    foreach my $file (@listef90) {
      $file =~ s/\.f90$// ;
      $types{$file} = 'f90' ;
      &getDependf90("$dir/$file.f90", $file, $dir) ;
    }
    foreach my $file (@listeF90) {
      $file =~ s/\.F90$// ;
      $types{$file} = 'F90' ;
      &getDependf90("$dir/$file.F90", $file, $dir) ;
    }
    (@listeall90) = sort glolevel &trySort(@listef90, @listeF90); 
    # List of main programs in fortran 90:
    (@lmains90) = () ;
    # List of not used programs in fortran 90
    (@lnouse90) = () ;
    # List of dependencies:
    (@orders) = () ;
    &sortallf90(@listef90, @listeF90) ;
    # Now sets the modules of subroutines:
    my(%oks) = () ;
    foreach (@lnouse90,@orders) {
      $oks{$_} = 1 ;
    }
    foreach (@listef90, @listeF90) {
      push(@othf90, $_) unless exists $oks{$_} ;
    }
    (%oks) = () ;
    print STDOUT "lnouse90=@lnouse90\n" if $debug ;
    print STDOUT "lmains90=@lmains90\n" if $debug ;
    print STDOUT "orders=@orders\n" if $debug ;
    print STDOUT "othf90=@othf90\n" if $debug ;
    # Open the target files for writing
    unless ($opt_b) {
      &printAll() ;
    }
  }
  if ($reuse) {
    if ($lasttime > $currenttime) {
      print STDOUT "WARNING(Z): $opt_t time $currenttime from ".localtime($currenttime)." older than last modified $lasttime on ".localtime($lasttime).", mkdep.lst may be incomplete (rerun fth.pl if unsure)\n" ;
    }
  }
  # print STDOUT "reuse=$reuse, lasttime=$lasttime, currenttime=$currenttime\n" ;
}
############################################################################
# Main loop on all directories
sub mainDirs {
  unless (opendir(DIR,".")) {
    print STDERR "Read permission denied on directory $dir\n" ;
    return ;
  }
  # Gets the actual disk location for includes from current local directory
  (%getrel) = () ;
  if (@rincludes) {
    foreach my $i (@rincludes) {
      # Gets a location that exists on disk from current directory.
      $getrel{$i} = &getAbs($dir,$i) ;
    }
  }
  # Gets the directory's file list
  &getList() ;
  closedir(DIR) ;

  $dir = "." ;
  # Gets the list of possible makefiles:
  if ($opt_R) {
     if (opendir(DIR,".")) {
       (@listemkf) = sort grep /makefile/i, readdir DIR ;
       closedir(DIR) ;
     }
  } else {
    (@listemkf) = sort grep /makefile/i, @liste ;
  }
  print STDOUT "listemkf=@listemkf on liste=@liste\n" if $debug ;
  # Looking for makefile rules:
  foreach my $file (@listemkf) {
    next if ($file =~ /[\%\~]$/) ;
    print STDOUT "Read processing rules for $file\n" if $debug ;
    &getRules("$file") ;
  }
  ($version,$date,$tag) = &getVersion() ;
  # print "tag=$tag\n" ;
  (%isamain) = () ;
  (%dones) = () ;
  (%types) = () ;
  
  (@listecpp) = sort grep /\.cpp$/, @liste ;
  foreach my $file (@listecpp) {
    $file =~ s/\.cpp$// ;
    $types{$file} = 'cpp' ;
    &getDependc("$file.cpp", $file, $dir) ;
  }
  (@listec) = sort grep /\.c$/, @liste ;
  foreach my $file (@listec) {
    $file =~ s/\.c$// ;
    $types{$file} = 'c' ;
    &getDependc("$file.c", $file, $dir) ;
  }
  (@listef) = grep /\.f$/, @liste ;
  foreach my $file (@listef) {
    if ($file =~ /^(.*)\.(f)$/) {
      $types{$1} = $2 ;
    }
    $file =~ s/\.f$// ;
    &getDependf("$file.f", $file, $dir) ;
  }
  (@listeF) = grep /\.F$/, @liste ;
  foreach my $file (@listeF) {
    if ($file =~ /^(.*)\.(F)$/) {
      $types{$1} = $2 ;
    }
    $file =~ s/\.F$// ;
    &getDependf("$file.F", $file, $dir) ;
  }
  (@listef90) = sort byname grep /\.f90$/, @liste ;
  (@listeF90) = sort byname grep /\.F90$/, @liste ;
  foreach my $file (@listef90) {
    $file =~ s/\.f90$// ;
    $types{$file} = 'f90' ;
    &getDependf90("$file.f90", $file, $dir) ;
  }
  foreach my $file (@listeF90) {
    $file =~ s/\.F90$// ;
    $types{$file} = 'F90' ;
    &getDependf90("$file.F90", $file, $dir) ;
  }
  (@listeall90) = sort glolevel &trySort(@listef90, @listeF90); 
  # List of main programs in fortran 90:
  (@lmains90) = () ;
  # List of not used programs in fortran 90
  (@lnouse90) = () ;
  # List of dependencies:
  (@orders) = () ;
  &sortallf90(@listef90, @listeF90) ;
  # Now sets the modules of subroutines:
  my(%oks) = () ;
  foreach (@lnouse90,@orders) {
    $oks{$_} = 1 ;
  }
  foreach (@listef90, @listeF90) {
    push(@othf90, $_) unless exists $oks{$_} ;
  }
  (%oks) = () ;
  print STDOUT "lnouse90=@lnouse90\n" if $debug ;
  print STDOUT "lmains90=@lmains90\n" if $debug ;
  print STDOUT "orders=@orders\n" if $debug ;
  print STDOUT "othf90=@othf90\n" if $debug ;
  # Open the target files for writing
  unless ($opt_b) {
    &printAll() ;
  }
  if ($reuse) {
    if ($lasttime > $currenttime) {
      print STDOUT "WARNING(Z): $opt_t time $currenttime from ".localtime($currenttime)." older than last modified $lasttime on ".localtime($lasttime).", mkdep.lst may be incomplete (rerun fth.pl if unsure)\n" ;
    }
  }
  # print STDOUT "reuse=$reuse, lasttime=$lasttime, currenttime=$currenttime\n" ;
}
############################################################################
if ($opt_R) {
  &mainDirs() ;
} else {
  &mainRout() ;
}
############################################################################
1;
