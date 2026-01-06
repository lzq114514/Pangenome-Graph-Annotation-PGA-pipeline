use strict;
use warnings;

open(IN,$ARGV[0]);
my $begin=1;
my $end=0;

my %hash;

while(<IN>){
    chomp;
    if(/^S/){
        my @arr = split(/\t/);
        $hash{$arr[1]} = length($arr[2]);
    }if(/^P/){
        my @arr = split(/\t/);
        $arr[2]=~s/-//g;
        $arr[2]=~s/\+//g;
        my @asd = split(/,/,$arr[2]);
        foreach my $node(@asd){
            $end+=$hash{$node};
            print $arr[1]."\t".$node."\t".$begin."\t".$end."\n";
            $begin+=$hash{$node};
        }
    }
	$begin=1;
	$end=0;
}
