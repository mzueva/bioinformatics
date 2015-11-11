
use libraries::Stat::Basic;
use Time::HiRes qw(gettimeofday tv_interval);

my $stat = new Stat::Basic;

my $file;

$before = gettimeofday;

open $file, '<', $ARGV[0];
my @lr = ();
$ind = 0;
while(<$file>)
{
	chomp;
	$lr[$ind] = $_;
	$ind++;
}

print "Number of entries: ";
print $ind;
print  "\n";
print  "\n";
$elapsed = gettimeofday( ) - $before;
print "Time for preparation: ";
print sprintf("%.2f", $elapsed*1000);
print  "\n";

$before = gettimeofday;
$res = $stat -> mean(\@lr);
$elapsed = gettimeofday( ) - $before;
print $stat -> mean(\@lr);
print  "\n";
print "Time for mean: ";
print sprintf("%.2f", $elapsed*1000);
print " ms";
print  "\n";

$before = gettimeofday;
print $stat -> median(\@lr);
print  "\n";
$elapsed = gettimeofday( ) - $before;
print "Time for median: ";
print sprintf("%.2f", $elapsed*1000);
print " ms";
print  "\n";

$before = gettimeofday;

print $stat -> max(\@lr);
print  "\n";
$elapsed = gettimeofday( ) - $before;
print "Time for max: ";
print sprintf("%.2f", $elapsed*1000);
print " ms";
print  "\n";

$before = gettimeofday;
print $stat -> min(\@lr);
print  "\n";
$elapsed = gettimeofday( ) - $before;
print "Time for min: ";
print sprintf("%.2f", $elapsed*1000);
print " ms";
print  "\n";

$before = gettimeofday;
print $stat -> sum(\@lr);
print  "\n";
$elapsed = gettimeofday( ) - $before;
print "Time for sum: ";
print sprintf("%.2f", $elapsed*1000);
print " ms";
print  "\n";

$before = gettimeofday;
print $stat -> prctile(\@lr, 10, 0);
print  "\n";
$elapsed = gettimeofday( ) - $before;
print "Time for prctile: ";
print sprintf("%.2f", $elapsed*1000);
print " ms";
print  "\n";

