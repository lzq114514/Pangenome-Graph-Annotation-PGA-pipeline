use strict;
use warnings;

# 打开输入和输出文件
open(IN, $ARGV[0]) or die "Cannot open input file: $!";
open(REF, ">$ARGV[1]") or die "Cannot open ref output file: $!";
open(ALT, ">$ARGV[2]") or die "Cannot open alt output file: $!";

# 创建哈希来存储样本信息
my %hash;
my $begin = 0;
my $end = 0;
my @AT;
my %printed_gt; # 用于记录已经输出的基因型

while (<IN>) {
    chomp;
    
    # 处理标题行，记录样本名称
    if (/^#CHROM/) {
        my @arr = split(/\t/);
        for (my $i = 9; $i < @arr; $i++) {
            $hash{$i} = $arr[$i];
        }
    } elsif (/^#/) { 
        next; 
    }

    my @arr = split(/\t/);

    # 检查 $arr[1] 是否为数字
    next unless $arr[1] =~ /^\d+$/;

    # 输出 REF 区域
    if (defined($arr[3])) {
        $begin = $arr[1];
        $end = $arr[1] + length($arr[3]);
        print REF $arr[0] . "\t" . $begin . "\t" . $end . "\n";
    }

    # 输出 ALT 区域
    if (defined($arr[4])) {
        # 从 INFO 列中提取 AT 信息
        $arr[7] =~ s/AC.*AT=//;
        $arr[7] =~ s/;NS.*//;
        @AT = split(/,/, $arr[7]);
        
        # 清空已输出基因型记录
        %printed_gt = ();
        
        # 遍历样本基因型信息
        for (my $j = 9; $j <= $#arr; $j++) {
            # 跳过缺失数据
            next if $arr[$j] eq ".";
            
            # 只处理单个数字的基因型
            if ($arr[$j] =~ /^(\d+)$/) {
                my $gt = $1;
                # 跳过0
                next if $gt == 0;
                # 如果这个基因型已经输出过，跳过
                next if exists $printed_gt{$gt};
                
                # 检查AT是否存在
                if (defined $AT[$gt - 1]) {
                    print ALT $hash{$j} . "\t" . $AT[$gt - 1] . "\n";
                    $printed_gt{$gt} = 1; # 标记这个基因型已经输出
                }
            }
        }
    }

    # 重置开始和结束位置
    $begin = 0;
    $end = 0;
}

# 关闭文件句柄
close(IN);
close(REF);
close(ALT);

print "筛选完成！\n";