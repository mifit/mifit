#!/usr/bin/perl

open(FOO,$ARGV[0]);
$level=0;
$seg=0;
$count=0;
$name[2]= $ARGV[2] ."_pref";
$name[3]= $ARGV[2] . "_allowed";
$printed=0;

while(<FOO>)
{
  if ($_ =~ m/vectorlist/) {
    $level++;
    if ($count>0) {
      $count=-1;
    }
  }
  if ($level < 1) {
    next;
  }

  @w=split(/ /,$_);
  if ($w[0] =~ m/\{level/) {
    if ($printed==1) {
      # close up last segment
      printf("%9.3ff\n};\n",0);
      $count++;
      $printed=0;
    }
    $printed=1;
    printf("static float %s_%s%d[] =  {\n",$name[$level],$ARGV[1],$count+1);
    if ($ARGV[1] =~ m/x/) {
      printf("%9.3ff, ",$w[2]);
    } else {
      printf("%9.3ff, ",$w[3]);
    }
    $seg=1;
  } elsif ($w[0] =~ m:\{\"\}L:) {
    if ($ARGV[1] =~ m/x/) {
      printf("%9.3ff, ",$w[1]);
    } else {
      printf("%9.3ff, ",$w[2]);
    }
    $seg++;
  }
  if ($seg>6) {
    printf("\n");
    $seg=0;
  }
}
printf("%9.3ff\n};\n",0);
