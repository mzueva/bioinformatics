
use libraries::Stat::Basic;
use Time::HiRes qw(gettimeofday);

my $stat = new Stat::Basic;

my $file;
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

$before = gettimeofday;
print "Result mean(): ";
print $stat -> mean(\@lr);
print  "\n";
$elapsed = gettimeofday() - $before;
print "Time mean(): ";
print $elapsed;
print " ms";
print  "\n";
print  "\n";


$before = gettimeofday;
print "Result median() ";
print $stat -> median(\@lr);
print  "\n";
$elapsed = gettimeofday() - $before;
print "Time median(): ";
print $elapsed;
print " ms";
print  "\n";
print  "\n";

$before = gettimeofday;
print "Result prctile() ";
print $stat -> prctile(\@lr, 0.1, 0);
print  "\n";
$elapsed = gettimeofday() - $before;
print "Time prctile(): ";
print $elapsed;
print " ms";
print  "\n";
print  "\n";


$before = gettimeofday;
print "Result max() ";
print $stat -> max(\@lr);
print  "\n";
$elapsed = gettimeofday() - $before;
print "Time max(): ";
print $elapsed;
print " ms";
print  "\n";
print  "\n";

$before = gettimeofday;
print "Result min() ";
print $stat -> min(\@lr);
print  "\n";
$elapsed = gettimeofday() - $before;
print "Time min(): ";
print $elapsed;
print " ms";
print  "\n";
print  "\n";


$before = gettimeofday;
print "Result sum() ";
print $stat -> sum(\@lr);
print  "\n";
$elapsed = gettimeofday() - $before;
print "Time sum(): ";
print $elapsed;
print " ms";
print  "\n";
print  "\n";

