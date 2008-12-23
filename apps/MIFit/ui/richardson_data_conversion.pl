#!/usr/bin/perl

# convert filename to array name
$name=$ARGV[0];
$name =~ s!\.data!_data!g;
$name =~ s!\-!_!g;
$name =~ s:\.\/::g;
$name =~ s!rama500!richardson!g;
printf("static const char *%s[180] = {\n",$name);


# get map of phi/psi data
open(INFILE,$ARGV[0]);
while(<INFILE>)
{
  if ($_ =~ m/\#/) {
    next;
  }
  @words=split(/ /,$_);
  $phi=$words[0];
  $psi=$words[1];
  $percent=$words[2];

  $i=int(($phi+180.0)/2.0);
  $j=int(($psi+180.0)/2.0);
  $phipsi{$i}{$j}=$percent;
}

# print it out
for ($j=180-1;$j>=0;$j--)
{
  printf "\t\"";
  for ($i=0;$i<180;$i++) {
    if ($phipsi{$i}{$j} >= 0.02) {
      printf "X";
    } elsif ($name == "richardson_general_data" && $phipsi{$i}{$j} >= 0.0005) {
      printf "x"
    } elsif ($phipsi{$i}{$j} >= 0.002) # special cases
      {
        printf "x"
      } else {
        printf " ";
      }
  }
  printf "\",\n";
}
printf "};\n"
