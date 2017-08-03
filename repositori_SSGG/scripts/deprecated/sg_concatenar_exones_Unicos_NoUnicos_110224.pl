#!/usr/bin/perl

sub usage 
{
    print "\nPARA QUE SIRVE: Sirve para concatenar secuencias en formato fasta (por ejemplo varios exones de un mismo transcrito), los transcritos con un unico exon tambien son reportados. No necesita orientacion (positiva o negativa) funciona para ambas  \n COMO SE USA: El input es un fichero fasta\n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

my %hash=();

my $line="";
my $entry="";
my $keys="";

open(FASTA,"<",@ARGV[0]);
while($line=<FASTA>)
{
        chomp($line);
        if($line=~/^>/)
	{
		$entry=$line;
	}
	else
	{
		if(!exists($hash{$entry}))
		{
			$hash{$entry}=$line;
		}
		else
		{
			$hash{$entry}=$hash{$entry}."".$line;
		}
	}
}
foreach $keys (sort  keys %hash)
{
	print "$keys\n$hash{$keys}\n";

}

