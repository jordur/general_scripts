#!/usr/bin/perl -w

open (LISTADO,"<",$ARGV[0]);
open (MIRBASEFASTA,"<",$ARGV[1]);

while ($lines=<MIRBASEFASTA>)
{
	chomp ($lines);
	push (@mirbasetmp,$lines);
}
while ($linea=<LISTADO>)
{
	chomp ($linea);
	@split_lalinea=  split ("-",$linea);
	for ($i=0;$i<=$#mirbasetmp;$i++)
	{
		$primer_char = substr ($mirbasetmp[$i],0,1);
		@split_linea= split ("-",$mirbasetmp[$i]);
		if ($primer_char eq ">")
		{
			$string = $split_linea[0] . "-" . $split_linea[1] . "-" . $split_linea[2];
			$string =~s/miR/mir/;
			$string =~s/>//;
			#print $strings,"\t",$string,"\n";
			if ($strings eq $string)
			{
			
				print $mirbasetmp[$i],"\n";
				print $mirbasetmp[$i+1],"\n";
			}
		}
	}
    
}
