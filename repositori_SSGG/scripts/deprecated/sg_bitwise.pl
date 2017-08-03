#!/usr/bin/perl -w

$num=$ARGV[0];

$resultado=dec2bin($num);
print $resultado,"\n\n\n";
$largo = length($resultado);
$reverso = reverse $resultado;
sub dec2bin  
{
    $numero =  unpack("B32", pack("N", shift));
    $numero =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $numero;
}

#$contador=0;
for ($i=0; $i< $largo; $i++)
{
	$caracter = substr($reverso,$i,1);
	if ($i==0)
	{
		if ($caracter ==1)
		{
			 print $caracter,"\tpaired\n";
		}
		else
		{
			print $caracter,"\tunpaired\n";
		}
	}
	elsif ($i==1)
	{
		if ($caracter ==1)
                {
                         print $caracter,"\tproper_pair\n";
                }
                else
                {
                        print $caracter,"\tno proper_pair\n";
 
                }
	}
	elsif ($i==2)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tunmapped\n";
                }
                else
                {
                        print $caracter,"\tmapped\n";
                }
        }
	elsif ($i==3)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tunmapped_mate\n";
                }
                else
                {
                        print $caracter,"\tmapped_mate\n";
                }
        }
	elsif ($i==4)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\treverse\n";
                }
                else
                {
                        print $caracter,"\tno reverse\n";
                }
        }
	elsif ($i ==5)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tmate reverse\n";
                }
                else
                {
                        print $caracter,"\tmate forward\n";
                }
        }
	elsif ($i==6)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tfirst in pair\n";
                }
                else
                {
                        print $caracter,"\tnot first in pair\n";
                }
        }
	elsif ($i==7)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tsecond in pair\n";
                }
                else
                {
                        print $caracter,"\tnot second in pair\n";
                }
        }
	elsif ($i==8)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tsplitted\n";
                }
                else
                {
                        print $caracter,"\tunsplitted\n";
                }
        }
	elsif ($i==9)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tfail_QV\n";
                }
                else
                {
                        print $caracter,"\tOK_QV\n";
                }
        }
	elsif ($i==11)
        {
                if ($caracter ==1)
                {
                         print $caracter,"\tPCR_duplicate\n";
                }
                else
                {
                        print $caracter,"\tnon_duplicate\n";
                }
        }


}



#if ($largo ==1)
#{
#	if ($resultado ==1)
#	{
#		print $resultado,"\tpaired\n";
#	}

#	else
#	{
#		print $resultado,"\tunpaired\n";
#	}

#}

#elsif ($largo==2)
#{
	

#}


